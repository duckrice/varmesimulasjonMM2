"""Åpner varmesimuleringen i et eget vindu.

Kjør med:  .venv\\Scripts\\python.exe simulering.py
"""

import matplotlib

matplotlib.rcParams["toolbar"] = "None"

import matplotlib.pyplot as plt  # noqa: E402

from varme import MATERIALER, Varmesimulering  # noqa: E402

if __name__ == "__main__":
    # modus: "1D" (stav) eller "2D" (plate). lengde: stavens lengde og platens sidekant i meter.
    simulering = Varmesimulering(MATERIALER["kobber"], modus="1D", lengde=1.0)
    plt.show()
