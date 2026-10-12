"""Løsere for varmelikningen  du/dt = alfa · ∇²u  med isolerte kanter.

* `simuler_eksplisitt`: numerisk løsning med endelige differanser (FTCS),
  for en stav (1D) eller en plate (2D).
* `FourierCosinusrekke`: analytisk løsning som Fourier-cosinusrekke for en stav.
* `FourierCosinusrekke2D`: tilsvarende dobbel cosinusrekke for en plate.
"""

import math
from functools import lru_cache

import numpy as np
from scipy import sparse

from .rutenett import Rutenett


@lru_cache(maxsize=8)  # samme rutenett gir samme matrise
def laplace_matrise(rutenett: Rutenett) -> sparse.csr_matrix:
    """Diskret Laplace-operator ∇²u med isolerte kanter, som en glissen
    matrise M slik at ∇²u = M @ u.ravel().

    Operatoren beskriver varmestrøm mellom naboceller: strømmen over flaten
    mellom celle i og i+1 er proporsjonal med (u[i+1] - u[i]). Over ytterkantene
    er det ingen nabo og dermed ingen strøm, altså isolert.

    I 1D er M tridiagonal: 1, -2, 1 på hver rad (delt på dx²), men -1 på
    diagonalen i endecellene, siden de bare har én nabo. I 2D summeres
    1D-matrisene for x- og y-retningen med Kronecker-produkt.
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

    Kravet er  alfa · dt · Σ(1/dx²) ≤ 1/2.  I 1D blir det  alfa·dt/dx² ≤ 1/2,
    og i 2D med like celler  alfa·dt/dx² ≤ 1/4.
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


def cosinus_operator(kanter: np.ndarray, antall_ledd: int) -> np.ndarray:
    """Matrise C slik at a = C @ u gir cosinus-koeffisientene a_0 … a_{N-1}
    langs én akse, når u er konstant innenfor hver celle:

        a_0 = (1/L) ∫ u dx
        a_n = (2/L) ∫ u · cos(nπx/L) dx

    `kanter` er cellekantene fra 0 til L. Integralet over celle i regnes ut
    eksakt: ∫ cos(kx) dx = (sin(k·x_{i+1/2}) - sin(k·x_{i-1/2})) / k.
    """
    L = kanter[-1]
    C = np.empty((antall_ledd, len(kanter) - 1))
    C[0] = np.diff(kanter) / L
    k = np.arange(1, antall_ledd)[:, np.newaxis] * np.pi / L
    C[1:] = (2.0 / L) * np.diff(np.sin(k * kanter), axis=1) / k
    return C


class FourierCosinusrekke:
    """Analytisk løsning for en isolert stav av lengde L:

        u(x, t) = a_0 + Σ_{n≥1} a_n · cos(nπx/L) · exp(-alfa·(nπ/L)²·t)

    Isolerte ender (du/dx = 0 i x = 0 og x = L) gir en cosinusrekke. Hvert ledd
    dør ut eksponentielt, og ledd med høyere n (raskere svingninger) dør ut
    raskere. Til slutt står bare a_0 igjen, som er snittemperaturen.

    Koeffisientene regnes ut eksakt fra celleverdiene i u0, se `cosinus_operator`.
    """

    def __init__(self, u0: np.ndarray, rutenett: Rutenett, alfa: float, antall_ledd: int = 200):
        self.L = rutenett.lengder[0]
        self.alfa = alfa
        self.n = np.arange(antall_ledd)
        self.a = cosinus_operator(rutenett.kanter[0], antall_ledd) @ u0

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


class FourierCosinusrekke2D:
    """Analytisk løsning for en isolert plate med sidekanter Lx og Ly:

        u(x, y, t) = Σ_m Σ_n a_mn · cos(mπx/Lx) · cos(nπy/Ly)
                               · exp(-alfa·π²·(m²/Lx² + n²/Ly²)·t)

    Dette er produktet av to cosinusrekker, én i hver retning. Leddet a_00 er
    snittemperaturen. `u0` har x langs første akse og y langs andre.
    """

    def __init__(self, u0: np.ndarray, rutenett: Rutenett, alfa: float, antall_ledd: int = 60):
        self.alfa = alfa
        kanter_x, kanter_y = rutenett.kanter
        self.kx = np.arange(antall_ledd) * np.pi / rutenett.lengder[0]
        self.ky = np.arange(antall_ledd) * np.pi / rutenett.lengder[1]
        # Koeffisientene: cosinus-operatoren brukt langs x, og så langs y
        self.a = cosinus_operator(kanter_x, antall_ledd) @ u0 @ cosinus_operator(kanter_y, antall_ledd).T

    def evaluer(self, x: np.ndarray, y: np.ndarray, t: float, antall_ledd: int | None = None) -> np.ndarray:
        """Summen av de første `antall_ledd` × `antall_ledd` leddene, i punktene (x_i, y_j)."""
        kx, ky = self.kx[:antall_ledd], self.ky[:antall_ledd]
        demping = np.exp(-self.alfa * t * (kx[:, np.newaxis] ** 2 + ky[np.newaxis, :] ** 2))
        a = self.a[: len(kx), : len(ky)] * demping
        return np.cos(np.outer(x, kx)) @ a @ np.cos(np.outer(ky, y))
