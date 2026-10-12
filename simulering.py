"""Åpner varmesimuleringen i et eget vindu.

Kjør med:  .venv\\Scripts\\python.exe simulering.py
"""

import matplotlib

matplotlib.rcParams["toolbar"] = "None"

import matplotlib.pyplot as plt  # noqa: E402

from varme import MATERIALER, InteraktivStav  # noqa: E402

if __name__ == "__main__":
    stav = InteraktivStav(MATERIALER["kobber"])
    plt.show()
