from __future__ import annotations

import numpy as np

from .utils import get_opt


def calculate_ef(enk: np.ndarray, u: float) -> float:
    """Estimate Fermi energy by filling ratio u in [0,1]."""
    if not (0.0 <= u <= 1.0):
        raise ValueError(f"u must be within [0,1], got {u}")
    flat = np.ravel(enk)
    n_occ = int(np.ceil(flat.size * u))
    if n_occ < 1:
        n_occ = 1
    # partition is O(N) and avoids full sort
    idx = n_occ - 1
    part = np.partition(flat, idx)
    return float(part[idx])


def extract_spinvalley_blocks(h: np.ndarray, n_layer: int) -> dict[str, np.ndarray]:
    """Extract K_up/Kp_up/K_dn/Kp_dn blocks from H(nkx,nky,8N,8N)."""
    nkx, nky, dim1, dim2 = h.shape
    if dim1 != dim2:
        raise ValueError("h must be (...,Nh,Nh)")
    dim_layer = 2 * int(n_layer)
    idx_k_up = slice(0, dim_layer)
    idx_kp_up = slice(dim_layer, 2 * dim_layer)
    idx_k_dn = slice(2 * dim_layer, 3 * dim_layer)
    idx_kp_dn = slice(3 * dim_layer, 4 * dim_layer)

    return {
        "K_up": h[:, :, idx_k_up, idx_k_up],
        "Kp_up": h[:, :, idx_kp_up, idx_kp_up],
        "K_dn": h[:, :, idx_k_dn, idx_k_dn],
        "Kp_dn": h[:, :, idx_kp_dn, idx_kp_dn],
    }


def diag_mesh_from_hk(hk: np.ndarray, opts: dict | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Diagonalize H(k) on mesh.

    hk: (nb, nb, Nkx, Nky)
    returns:
      U: (nb, nb_sel, Nkx, Nky)
      E: (Nkx, Nky, nb_sel)
    """
    if hk.ndim != 4 or hk.shape[0] != hk.shape[1]:
        raise ValueError("hk must have shape (nb, nb, Nkx, Nky)")

    nb, _, nkx, nky = hk.shape
    band_list = get_opt(opts, "band_list", list(range(1, nb + 1)))
    band_idx = np.asarray(band_list, dtype=np.int64) - 1
    if np.any(band_idx < 0) or np.any(band_idx >= nb):
        raise ValueError("band_list out of range for hk")
    nb_sel = band_idx.size

    u = np.zeros((nb, nb_sel, nkx, nky), dtype=np.complex128)
    e = np.zeros((nkx, nky, nb_sel), dtype=np.float64)

    for ix in range(nkx):
        for iy in range(nky):
            h = hk[:, :, ix, iy]
            h = 0.5 * (h + h.conj().T)
            evals, vecs = np.linalg.eigh(h)
            e[ix, iy, :] = np.real(evals[band_idx])
            u[:, :, ix, iy] = vecs[:, band_idx]
    return u, e

