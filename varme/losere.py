"""Løsere for varmelikningen  du/dt = alfa · ∇²u  med isolerte kanter.

* `simuler_eksplisitt`: numerisk løsning med endelige differanser (FTCS),
  fungerer i 1D, 2D og 3D.
* `FourierCosinusrekke`: analytisk løsning som Fourier-cosinusrekke (1D).
"""

import math
from functools import lru_cache

import numpy as np
from scipy import sparse

from .rutenett import Rutenett


def laplace(u: np.ndarray, dx: tuple[float, ...]) -> np.ndarray:
    """Diskret Laplace-operator ∇²u med isolerte kanter (du/dn = 0).

    Regnes ut som varmestrøm mellom naboceller: strømmen over flaten mellom
    celle i og i+1 er proporsjonal med (u[i+1] - u[i]). Den legges til celle i
    og trekkes fra celle i+1, så varme verken skapes eller forsvinner. Over
    ytterkantene er det ingen flate og dermed ingen strøm, altså isolert.
    Til sammen blir dette (u[i+1] - 2·u[i] + u[i-1]) / h² i det indre.
    """
    resultat = np.zeros_like(u)
    for akse, h in enumerate(dx):
        strom = np.diff(u, axis=akse) / h**2
        r = np.moveaxis(resultat, akse, 0)  # visning: endringer går rett i `resultat`
        s = np.moveaxis(strom, akse, 0)
        r[:-1] += s
        r[1:] -= s
    return resultat


@lru_cache(maxsize=8)  # samme rutenett gir samme matrise
def laplace_matrise(rutenett: Rutenett) -> sparse.csr_matrix:
    """Samme operator som `laplace`, men som en (glissen) matrise M, slik at
    ∇²u = M @ u.ravel(). Matrisen er mye raskere å bruke i tidsløkka.

    I 1D er M tridiagonal: 1, -2, 1 på hver rad, men -1 på diagonalen i
    endecellene, siden de bare har én nabo (isolert kant). I flere dimensjoner
    summeres 1D-matrisene langs hver akse med Kronecker-produkt.
    """
    M = sparse.csr_matrix((np.prod(rutenett.form),) * 2)
    for akse, (n, h) in enumerate(zip(rutenett.form, rutenett.dx)):
        diagonal = np.full(n, -2.0)
        diagonal[[0, -1]] = -1.0
        D = sparse.diags([np.ones(n - 1), diagonal, np.ones(n - 1)], [-1, 0, 1]) / h**2
        faktorer = [sparse.identity(m) for m in rutenett.form]
        faktorer[akse] = D
        ledd = faktorer[0]
        for f in faktorer[1:]:
            ledd = sparse.kron(ledd, f)
        M = M + ledd
    return sparse.csr_matrix(M)


def stabilt_tidssteg(rutenett: Rutenett, alfa: float) -> float:
    """Største tidssteg der den eksplisitte metoden er stabil.

    Kravet er  alfa · dt · Σ(1/dx²) ≤ 1/2,  som i 1D blir  alfa·dt/dx² ≤ 1/2.
    """
    return 0.5 / (alfa * sum(1.0 / h**2 for h in rutenett.dx))


def simuler_eksplisitt(
    u0: np.ndarray,
    rutenett: Rutenett,
    alfa: float,
    tider,
    sikkerhet: float = 0.9,
) -> np.ndarray:
    """Løser varmelikningen med eksplisitt Euler i tid og sentrale differanser i rom.

        u(t + dt) = u(t) + dt · alfa · ∇²u(t)

    `u0` er temperaturen ved t = 0, og `tider` er stigende tidspunkt (≥ 0) der
    temperaturen skal lagres. Avstanden mellom tidspunktene kan variere: mellom
    hvert par tas så mange like steg som trengs for å holde metoden stabil.

    Returnerer bildene, med tid langs første akse.
    """
    dt_maks = sikkerhet * stabilt_tidssteg(rutenett, alfa)
    tider = np.asarray(tider, dtype=float)
    if np.any(np.diff(tider) < 0) or tider[0] < 0:
        raise ValueError("Tidspunktene må være stigende og ≥ 0.")

    M = laplace_matrise(rutenett)
    bilder = np.empty((len(tider), *u0.shape))
    u = np.array(u0, dtype=float).ravel()
    t = 0.0
    for i, t_neste in enumerate(tider):
        varighet = t_neste - t
        if varighet > 0:
            antall_steg = math.ceil(varighet / dt_maks)
            c = varighet / antall_steg * alfa
            for _ in range(antall_steg):
                u = u + c * (M @ u)
        t = t_neste
        bilder[i] = u.reshape(u0.shape)
    return bilder


class FourierCosinusrekke:
    """Analytisk løsning for en isolert stav av lengde L:

        u(x, t) = a_0 + Σ_{n≥1} a_n · cos(nπx/L) · exp(-alfa·(nπ/L)²·t)

    Isolerte ender (du/dx = 0 i x = 0 og x = L) gir en cosinusrekke. Hvert ledd
    dør ut eksponentielt, og ledd med høyere n (raskere svingninger) dør ut
    raskere. Til slutt står bare a_0 igjen, som er snittemperaturen.

    Koeffisientene regnes ut eksakt fra celleverdiene i u0, ved å se på u0 som
    konstant innenfor hver celle:

        a_0 = (1/L) ∫ u0 dx
        a_n = (2/L) ∫ u0 · cos(nπx/L) dx
    """

    def __init__(self, u0: np.ndarray, rutenett: Rutenett, alfa: float, antall_ledd: int = 200):
        if rutenett.dimensjon != 1:
            raise NotImplementedError("Fourier-løsningen er foreløpig bare laget for 1D.")
        self.L = rutenett.lengder[0]
        self.alfa = alfa
        self.n = np.arange(antall_ledd)

        kanter = rutenett.kanter[0]
        k = self.n[1:, np.newaxis] * np.pi / self.L
        # ∫ over celle i av cos(kx) dx = (sin(k·x_{i+1/2}) - sin(k·x_{i-1/2})) / k
        integraler = np.diff(np.sin(k * kanter), axis=1) / k

        self.a = np.empty(antall_ledd)
        self.a[0] = np.sum(u0 * np.diff(kanter)) / self.L
        self.a[1:] = (2.0 / self.L) * integraler @ u0

    def ledd(self, n: int, x: np.ndarray, t: float) -> np.ndarray:
        """Ledd nummer n alene: a_n · cos(nπx/L) · exp(-alfa·(nπ/L)²·t)."""
        k = n * np.pi / self.L
        return self.a[n] * np.cos(k * x) * np.exp(-self.alfa * k**2 * t)

    def evaluer(self, x: np.ndarray, t: float, antall_ledd: int | None = None) -> np.ndarray:
        """Summen av de første `antall_ledd` leddene (n = 0, 1, …, antall_ledd - 1)."""
        n = self.n[:antall_ledd]
        k = n * np.pi / self.L
        demping = np.exp(-self.alfa * k**2 * t)
        return (self.a[: len(n)] * demping) @ np.cos(np.outer(k, x))
