"""Materialdata for metaller (omtrentlige verdier ved romtemperatur)."""

from dataclasses import dataclass


@dataclass(frozen=True)
class Materiale:
    navn: str
    varmeledningsevne: float  # k   [W/(m·K)]
    tetthet: float  # rho [kg/m³]
    spesifikk_varme: float  # c_p [J/(kg·K)]

    @property
    def alfa(self) -> float:
        """Termisk diffusivitet alfa = k / (rho · c_p)  [m²/s]."""
        return self.varmeledningsevne / (self.tetthet * self.spesifikk_varme)


MATERIALER = {
    "kobber": Materiale("Kobber", 401.0, 8960.0, 385.0),
    "aluminium": Materiale("Aluminium", 237.0, 2700.0, 897.0),
    "messing": Materiale("Messing", 109.0, 8530.0, 380.0),
    "jern": Materiale("Jern", 80.0, 7874.0, 449.0),
    "rustfritt stål": Materiale("Rustfritt stål", 16.0, 8000.0, 500.0),
}
