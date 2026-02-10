from __future__ import annotations

import numpy as np


def gauge_fix_parallel_transport(u: np.ndarray) -> np.ndarray:
    """Band-wise phase smoothing over 2D mesh.

    u shape: (nb, nb_sel, Nkx, Nky)
    """
    ug = np.array(u, copy=True)
    _, nb_sel, nkx, nky = ug.shape

    # Sweep along kx
    for iy in range(nky):
        for ix in range(1, nkx):
            uprev = ug[:, :, ix - 1, iy]
            ucur = ug[:, :, ix, iy].copy()
            for n in range(nb_sel):
                ov = np.vdot(uprev[:, n], ucur[:, n])
                ph = ov / max(abs(ov), 1e-30)
                ucur[:, n] /= ph
            ug[:, :, ix, iy] = ucur

    # Sweep along ky
    for ix in range(nkx):
        for iy in range(1, nky):
            uprev = ug[:, :, ix, iy - 1]
            ucur = ug[:, :, ix, iy].copy()
            for n in range(nb_sel):
                ov = np.vdot(uprev[:, n], ucur[:, n])
                ph = ov / max(abs(ov), 1e-30)
                ucur[:, n] /= ph
            ug[:, :, ix, iy] = ucur
    return ug

