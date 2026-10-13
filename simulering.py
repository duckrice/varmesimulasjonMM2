"""Åpner varmesimuleringen i et eget vindu.

Kjør med:  .venv\\Scripts\\python.exe simulering.py
"""

import os

# Matrisene i dette prosjektet er små. numpys matriseprodukt bytter likevel til
# mange tråder når de passerer en viss størrelse (f.eks. 200 Fourier-ledd), og da
# går det rundt 50 ganger saktere enn med én tråd. Må settes før numpy lastes.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import matplotlib  # noqa: E402

matplotlib.rcParams["toolbar"] = "None"

import matplotlib.pyplot as plt  # noqa: E402

from varme import MATERIALER, Varmesimulering  # noqa: E402

if __name__ == "__main__":
    # modus: "1D" (stav) eller "2D" (plate). lengde: stavens lengde og platens sidekant i meter.
    simulering = Varmesimulering(MATERIALER["kobber"], modus="1D", lengde=1.0)
    plt.show()
