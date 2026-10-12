"""Starttemperaturer, dvs. u(x, 0)."""

import numpy as np

from .rutenett import Rutenett


def varmt_segment(
    rutenett: Rutenett,
    sentrum,
    bredde,
    T_omgivelse: float,
    T_varm: float,
) -> np.ndarray:
    """Hele domenet holder T_omgivelse, bortsett fra et varmt område med T_varm.

    `sentrum` og `bredde` er tall i 1D, eller tupler med én verdi per akse i 2D/3D.
    Cellene der cellesenteret ligger innenfor området blir varme.
    """
    sentrum = np.atleast_1d(sentrum)
    bredde = np.broadcast_to(np.atleast_1d(bredde), sentrum.shape)

    inni = np.ones(rutenett.form, dtype=bool)
    for akse, x in enumerate(rutenett.sentre):
        innenfor = np.abs(x - sentrum[akse]) <= bredde[akse] / 2
        # Sørg for at minst én celle blir varm selv om segmentet er smalere enn en celle
        if not innenfor.any():
            innenfor[np.argmin(np.abs(x - sentrum[akse]))] = True
        form = [1] * rutenett.dimensjon
        form[akse] = -1
        inni &= innenfor.reshape(form)

    return np.where(inni, T_varm, T_omgivelse).astype(float)
