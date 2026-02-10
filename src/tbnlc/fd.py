from __future__ import annotations

import numpy as np


def grid_step_vectors(kx: np.ndarray, ky: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Infer dk1 (ix direction) and dk2 (iy direction) from mesh arrays."""
    nkx, nky = kx.shape

    dkx1 = np.diff(kx[:, 0])
    dky1 = np.diff(ky[:, 0])
    nz1 = np.where((np.abs(dkx1) + np.abs(dky1)) > 0)[0]
    if nkx < 2 or nz1.size == 0:
        raise ValueError("Cannot infer dk1 from k-grid.")
    i1 = int(nz1[0])
    dk1 = np.array([dkx1[i1], dky1[i1]], dtype=np.float64)

    dkx2 = np.diff(kx[0, :])
    dky2 = np.diff(ky[0, :])
    nz2 = np.where((np.abs(dkx2) + np.abs(dky2)) > 0)[0]
    if nky < 2 or nz2.size == 0:
        raise ValueError("Cannot infer dk2 from k-grid.")
    i2 = int(nz2[0])
    dk2 = np.array([dkx2[i2], dky2[i2]], dtype=np.float64)
    return dk1, dk2


def skewgrid_ops_from_k(
    kx: np.ndarray,
    ky: np.ndarray,
    periodic_fd: bool = False,
    trim_boundary: bool = True,
) -> dict:
    nkx, nky = kx.shape
    if ky.shape != (nkx, nky):
        raise ValueError("kx and ky must have same shape")

    dk1, dk2 = grid_step_vectors(kx, ky)
    j = np.column_stack([dk1, dk2])
    inv_jt = np.linalg.inv(j.T)
    det_j = float(np.linalg.det(j))
    d2k = abs(det_j)
    w_k = d2k / (2.0 * np.pi) ** 2

    if trim_boundary and not periodic_fd:
        ix_list = np.arange(1, nkx - 1, dtype=np.int64)
        iy_list = np.arange(1, nky - 1, dtype=np.int64)
        mask_k = np.zeros((nkx, nky), dtype=np.float64)
        if ix_list.size > 0 and iy_list.size > 0:
            mask_k[np.ix_(ix_list, iy_list)] = 1.0
    else:
        ix_list = np.arange(nkx, dtype=np.int64)
        iy_list = np.arange(nky, dtype=np.int64)
        mask_k = np.ones((nkx, nky), dtype=np.float64)

    return {
        "nkx": nkx,
        "nky": nky,
        "dk1": dk1,
        "dk2": dk2,
        "J": j,
        "invJT": inv_jt,
        "detJ": det_j,
        "d2k": d2k,
        "w_k": w_k,
        "ix_list": ix_list,
        "iy_list": iy_list,
        "mask_k": mask_k,
        "mask_k4": mask_k.reshape(nkx, nky, 1, 1),
    }


def fd_du_skew(u: np.ndarray, inv_jt: np.ndarray, periodic_fd: bool) -> tuple[np.ndarray, np.ndarray]:
    """Cartesian derivatives for eigenvectors on skew mesh.

    u shape: (nb, nb_sel, Nkx, Nky)
    """
    _, _, nkx, nky = u.shape
    if periodic_fd:
        idxp1 = np.r_[1:nkx, 0]
        idxm1 = np.r_[nkx - 1, 0 : nkx - 1]
        idxp2 = np.r_[1:nky, 0]
        idxm2 = np.r_[nky - 1, 0 : nky - 1]
        du_k1 = (u[:, :, idxp1, :] - u[:, :, idxm1, :]) * 0.5
        du_k2 = (u[:, :, :, idxp2] - u[:, :, :, idxm2]) * 0.5
    else:
        du_k1 = np.zeros_like(u)
        du_k2 = np.zeros_like(u)
        if nkx > 2:
            du_k1[:, :, 1:-1, :] = (u[:, :, 2:, :] - u[:, :, :-2, :]) * 0.5
        if nky > 2:
            du_k2[:, :, :, 1:-1] = (u[:, :, :, 2:] - u[:, :, :, :-2]) * 0.5

    du_x = inv_jt[0, 0] * du_k1 + inv_jt[0, 1] * du_k2
    du_y = inv_jt[1, 0] * du_k1 + inv_jt[1, 1] * du_k2
    return du_x, du_y


def fd_tensor_skew(f: np.ndarray, inv_jt: np.ndarray, periodic_fd: bool) -> tuple[np.ndarray, np.ndarray]:
    """Central differences for tensor field f(ix,iy,...) on skew mesh."""
    sz = f.shape
    nkx, nky = sz[0], sz[1]
    rest = int(np.prod(sz[2:])) if len(sz) > 2 else 1
    f3 = f.reshape(nkx, nky, rest)

    if periodic_fd:
        idxp1 = np.r_[1:nkx, 0]
        idxm1 = np.r_[nkx - 1, 0 : nkx - 1]
        idxp2 = np.r_[1:nky, 0]
        idxm2 = np.r_[nky - 1, 0 : nky - 1]
        d_k1 = (f3[idxp1, :, :] - f3[idxm1, :, :]) * 0.5
        d_k2 = (f3[:, idxp2, :] - f3[:, idxm2, :]) * 0.5
    else:
        d_k1 = np.zeros_like(f3)
        d_k2 = np.zeros_like(f3)
        if nkx > 2:
            d_k1[1:-1, :, :] = (f3[2:, :, :] - f3[:-2, :, :]) * 0.5
        if nky > 2:
            d_k2[:, 1:-1, :] = (f3[:, 2:, :] - f3[:, :-2, :]) * 0.5

    dfx3 = inv_jt[0, 0] * d_k1 + inv_jt[0, 1] * d_k2
    dfy3 = inv_jt[1, 0] * d_k1 + inv_jt[1, 1] * d_k2
    return dfx3.reshape(sz), dfy3.reshape(sz)


def fd_rmn_skew(r_mn: np.ndarray, inv_jt: np.ndarray, periodic_fd: bool) -> np.ndarray:
    """r_mn: (Nkx,Nky,nb,nb,2) -> dr: (Nkx,Nky,nb,nb,2,2)."""
    if r_mn.ndim != 5 or r_mn.shape[4] != 2:
        raise ValueError("r_mn must have shape (Nkx,Nky,nb,nb,2)")
    nkx, nky, nb, _, nc = r_mn.shape
    r3 = r_mn.reshape(nkx, nky, nb * nb * nc)
    dx3, dy3 = fd_tensor_skew(r3, inv_jt, periodic_fd)
    dx = dx3.reshape(nkx, nky, nb, nb, nc)
    dy = dy3.reshape(nkx, nky, nb, nb, nc)
    dr = np.zeros((nkx, nky, nb, nb, nc, 2), dtype=r_mn.dtype)
    dr[:, :, :, :, :, 0] = dx
    dr[:, :, :, :, :, 1] = dy
    return dr

