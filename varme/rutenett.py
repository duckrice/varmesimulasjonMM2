"""Cellesentrert rutenett for en stav (1D) eller en plate (2D).

Domenet deles i like store celler, og temperaturen lagres i midten av hver
celle. Det gjør isolerte ender enkle å håndtere: kanten av domenet ligger
nøyaktig mellom to celler.
"""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Rutenett:
    lengder: tuple[float, ...]  # domenets utstrekning langs hver akse [m]
    antall_celler: tuple[int, ...]  # antall celler langs hver akse

    @classmethod
    def stav(cls, lengde: float, antall_celler: int) -> "Rutenett":
        """Et 1D-rutenett for en stav."""
        return cls((lengde,), (antall_celler,))

    @property
    def dimensjon(self) -> int:
        return len(self.lengder)

    @property
    def form(self) -> tuple[int, ...]:
        return self.antall_celler

    @property
    def dx(self) -> tuple[float, ...]:
        """Cellestørrelse langs hver akse [m]."""
        return tuple(L / n for L, n in zip(self.lengder, self.antall_celler))

    @property
    def sentre(self) -> tuple[np.ndarray, ...]:
        """Koordinatene til cellesentrene langs hver akse [m]."""
        return tuple((np.arange(n) + 0.5) * h for n, h in zip(self.antall_celler, self.dx))

    @property
    def kanter(self) -> tuple[np.ndarray, ...]:
        """Koordinatene til cellekantene langs hver akse (n + 1 verdier) [m]."""
        return tuple(np.linspace(0.0, L, n + 1) for L, n in zip(self.lengder, self.antall_celler))
