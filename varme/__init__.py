"""Simulering av varmelikningen i en stav (1D) eller en plate (2D)."""

from .losere import (
    FourierCosinusrekke,
    FourierCosinusrekke2D,
    FourierSinusrekke,
    FourierSinusrekke2D,
    laplace_matrise,
    simuler_eksplisitt,
    stabilt_tidssteg,
)
from .materialer import MATERIALER, Materiale
from .rutenett import Rutenett
from .startbetingelser import varmt_segment
from .visualisering import Varmesimulering

__all__ = [
    "MATERIALER",
    "FourierCosinusrekke",
    "FourierCosinusrekke2D",
    "FourierSinusrekke",
    "FourierSinusrekke2D",
    "Materiale",
    "Rutenett",
    "Varmesimulering",
    "laplace_matrise",
    "simuler_eksplisitt",
    "stabilt_tidssteg",
    "varmt_segment",
]
