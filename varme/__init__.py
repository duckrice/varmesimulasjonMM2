"""Simulering av varmelikningen (1D nå, laget for å kunne utvides til 2D/3D)."""

from .losere import FourierCosinusrekke, laplace, simuler_eksplisitt, stabilt_tidssteg
from .materialer import MATERIALER, Materiale
from .rutenett import Rutenett
from .startbetingelser import varmt_segment
from .visualisering import InteraktivStav

__all__ = [
    "MATERIALER",
    "FourierCosinusrekke",
    "InteraktivStav",
    "Materiale",
    "Rutenett",
    "laplace",
    "simuler_eksplisitt",
    "stabilt_tidssteg",
    "varmt_segment",
]
