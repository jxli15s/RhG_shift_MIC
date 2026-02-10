from __future__ import annotations

import numba as nb
import numpy as np


@nb.njit(cache=True, parallel=True)
def shift_integrate_kernel(
    e: np.ndarray,  # (Nkx,Nky,nb)
    f_n: np.ndarray,  # (Nkx,Nky,nb)
    r_mn: np.ndarray,  # (Nkx,Nky,nb,nb,2) complex
    r_cov: np.ndarray,  # (Nkx,Nky,nb,nb,2,2) complex
    eph_list: np.ndarray,  # (Nw,)
    ix_list: np.ndarray,
    iy_list: np.ndarray,
    pref: float,
    w_k: float,
    eta: float,
    sym_bc: int,
) -> np.ndarray:
    nb_sel = e.shape[2]
    n_w = eph_list.size
    sigma = np.zeros((2, 2, 2, n_w), dtype=np.float64)

    for comp in nb.prange(8):
        a = comp // 4
        b = (comp // 2) % 2
        c = comp % 2
        sig_comp = np.zeros(n_w, dtype=np.float64)
        for ixi in range(ix_list.size):
            ix = ix_list[ixi]
            for iyi in range(iy_list.size):
                iy = iy_list[iyi]
                for n in range(nb_sel):
                    fn = f_n[ix, iy, n]
                    en = e[ix, iy, n]
                    for m in range(nb_sel):
                        if m == n:
                            continue
                        d_e = e[ix, iy, m] - en
                        f_nm = fn - f_n[ix, iy, m]
                        if abs(f_nm) < 1e-14:
                            continue

                        z1 = r_mn[ix, iy, m, n, b] * r_cov[ix, iy, n, m, c, a]
                        i1 = z1.imag
                        if sym_bc == 1:
                            z2 = r_mn[ix, iy, m, n, c] * r_cov[ix, iy, n, m, b, a]
                            i_val = 0.5 * (i1 + z2.imag)
                        else:
                            i_val = i1

                        fac = pref * f_nm * i_val * w_k / np.pi
                        for iw in range(n_w):
                            de = d_e - eph_list[iw]
                            lor = eta / (de * de + eta * eta)
                            sig_comp[iw] += fac * lor
        sigma[a, b, c, :] = sig_comp
    return sigma


@nb.njit(cache=True, parallel=True)
def mic_integrate_kernel(
    d_e_mn: np.ndarray,  # (Nkx,Nky,nb,nb)
    f_nm: np.ndarray,  # (Nkx,Nky,nb,nb)
    dv_mn: np.ndarray,  # (Nkx,Nky,nb,nb,2)
    g2: np.ndarray,  # (Nkx,Nky,nb,nb,2,2)
    mask_total: np.ndarray,  # (Nkx,Nky,nb,nb)
    eph_list: np.ndarray,
    pref: float,
    w_k: float,
    eta: float,
) -> np.ndarray:
    n_w = eph_list.size
    nb_sel = d_e_mn.shape[2]
    out = np.zeros((2, 2, 2, n_w), dtype=np.float64)

    for comp in nb.prange(8):
        a = comp // 4
        b = (comp // 2) % 2
        c = comp % 2
        eta_comp = np.zeros(n_w, dtype=np.float64)

        for ix in range(d_e_mn.shape[0]):
            for iy in range(d_e_mn.shape[1]):
                for m in range(nb_sel):
                    for n in range(nb_sel):
                        w = mask_total[ix, iy, m, n]
                        if w == 0.0:
                            continue
                        kval = f_nm[ix, iy, m, n] * dv_mn[ix, iy, m, n, a] * g2[ix, iy, m, n, b, c] * w
                        if kval == 0.0:
                            continue
                        d_e = d_e_mn[ix, iy, m, n]
                        fac = pref * w_k * kval / np.pi
                        for iw in range(n_w):
                            de = d_e - eph_list[iw]
                            lor = eta / (de * de + eta * eta)
                            eta_comp[iw] += fac * lor
        out[a, b, c, :] = eta_comp

    return out

