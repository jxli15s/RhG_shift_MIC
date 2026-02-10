from __future__ import annotations

from typing import Any

import numpy as np

from .fd import fd_du_skew, fd_rmn_skew, fd_tensor_skew, skewgrid_ops_from_k
from .gauge import gauge_fix_parallel_transport
from .kernels import mic_integrate_kernel, shift_integrate_kernel
from .utils import get_opt


def _reshape_u_e(kx: np.ndarray, ky: np.ndarray, u: np.ndarray, e: np.ndarray) -> tuple[np.ndarray, np.ndarray, int, int, int, int]:
    nkx, nky = kx.shape
    if ky.shape != (nkx, nky):
        raise ValueError("kx and ky must have same shape")

    if u.ndim == 3:
        nb, nb_sel, nflat = u.shape
        if nflat != nkx * nky:
            raise ValueError("If u is 3D, u.shape[2] must be Nkx*Nky")
        u = u.reshape(nb, nb_sel, nkx, nky)
    elif u.ndim == 4:
        nb, nb_sel, ukx, uky = u.shape
        if ukx != nkx or uky != nky:
            raise ValueError("u must have shape (nb,nb_sel,Nkx,Nky)")
    else:
        raise ValueError("u must be 3D or 4D")

    if e.shape != (nkx, nky, nb_sel):
        if e.shape == (nkx * nky, nb_sel):
            e = e.reshape(nkx, nky, nb_sel)
        else:
            raise ValueError("e must have shape (Nkx,Nky,nb_sel) or (Nkx*Nky,nb_sel)")

    return u, e, nkx, nky, nb, nb_sel


def fermi_dirac(e: np.ndarray, ef: float, k_t: float) -> np.ndarray:
    if k_t <= 0:
        return (e < ef).astype(np.float64)
    x = np.clip((e - ef) / k_t, -700.0, 700.0)
    return 1.0 / (1.0 + np.exp(x))


def hermitianize_rmn(r_mn: np.ndarray) -> np.ndarray:
    return 0.5 * (r_mn + np.transpose(np.conjugate(r_mn), (0, 1, 3, 2, 4)))


def shift_current_plane_fd_energy_skew_from_ue(
    kx: np.ndarray,
    ky: np.ndarray,
    u: np.ndarray,
    e: np.ndarray,
    eph_list: np.ndarray,
    ef: float,
    k_t: float,
    eta: float,
    opts: Any | None = None,
    return_details: bool = False,
):
    periodic_fd = bool(get_opt(opts, "periodicFD", False))
    trim_boundary = bool(get_opt(opts, "trimBoundary", True))
    sym_bc = bool(get_opt(opts, "symBC", True))
    verbose = bool(get_opt(opts, "verbose", True))
    g_s = float(get_opt(opts, "g_s", 1.0))
    do_gauge_fix = bool(get_opt(opts, "doGaugeFix", True))
    save_intermediates = bool(get_opt(opts, "saveIntermediates", True))

    e_charge = 1.602176634e-19
    hbar_js = 1.054571817e-34
    pref = np.pi * g_s * (e_charge**2) / hbar_js

    u, e, nkx, nky, nb, nb_sel = _reshape_u_e(kx, ky, u, e)

    band_list = get_opt(opts, "band_list", None)
    if band_list is not None and len(band_list) > 0:
        bidx = np.asarray(band_list, dtype=np.int64) - 1
        if np.any(bidx < 0) or np.any(bidx >= nb_sel):
            raise ValueError(f"band_list out of range, nb_sel={nb_sel}")
        u = u[:, bidx, :, :]
        e = e[:, :, bidx]
        nb_sel = bidx.size

    eph_list = np.asarray(eph_list, dtype=np.float64).reshape(-1)
    ops = skewgrid_ops_from_k(kx, ky, periodic_fd=periodic_fd, trim_boundary=trim_boundary)

    if verbose:
        print(
            f"[tbnlc.shift] Nkx={nkx} Nky={nky} nb={nb} nb_sel={nb_sel} "
            f"periodicFD={int(periodic_fd)} trimBoundary={int(trim_boundary)} symBC={int(sym_bc)}"
        )

    ug = gauge_fix_parallel_transport(u) if do_gauge_fix else u
    du_x, du_y = fd_du_skew(ug, ops["invJT"], periodic_fd=periodic_fd)

    a_diag = np.zeros((nkx, nky, nb_sel, 2), dtype=np.float64)
    r_mn = np.zeros((nkx, nky, nb_sel, nb_sel, 2), dtype=np.complex128)
    for ix in ops["ix_list"]:
        for iy in ops["iy_list"]:
            u0 = ug[:, :, ix, iy]
            dux = du_x[:, :, ix, iy]
            duy = du_y[:, :, ix, iy]
            ax = 1j * np.diag(u0.conj().T @ dux)
            ay = 1j * np.diag(u0.conj().T @ duy)
            a_diag[ix, iy, :, 0] = np.real(ax)
            a_diag[ix, iy, :, 1] = np.real(ay)
            r_mn[ix, iy, :, :, 0] = 1j * (u0.conj().T @ dux)
            r_mn[ix, iy, :, :, 1] = 1j * (u0.conj().T @ duy)

    off = np.ones((nb_sel, nb_sel), dtype=np.float64) - np.eye(nb_sel, dtype=np.float64)
    r_mn *= off.reshape(1, 1, nb_sel, nb_sel, 1)
    r_mn = hermitianize_rmn(r_mn)
    r_mn *= off.reshape(1, 1, nb_sel, nb_sel, 1)

    dr = fd_rmn_skew(r_mn, ops["invJT"], periodic_fd=periodic_fd)

    r_cov = np.zeros_like(dr)
    for a in range(2):
        amm = a_diag[:, :, :, a].reshape(nkx, nky, nb_sel, 1)
        ann = a_diag[:, :, :, a].reshape(nkx, nky, 1, nb_sel)
        d_a = amm - ann
        for c in range(2):
            r_cov[:, :, :, :, c, a] = dr[:, :, :, :, c, a] - 1j * d_a * r_mn[:, :, :, :, c]

    f_n = fermi_dirac(e, ef, k_t).astype(np.float64)
    sigma = shift_integrate_kernel(
        e.astype(np.float64),
        f_n,
        r_mn.astype(np.complex128),
        r_cov.astype(np.complex128),
        eph_list.astype(np.float64),
        ops["ix_list"].astype(np.int64),
        ops["iy_list"].astype(np.int64),
        float(pref),
        float(ops["w_k"]),
        float(eta),
        1 if sym_bc else 0,
    )

    if not return_details:
        return sigma

    out = {
        "dk1": ops["dk1"],
        "dk2": ops["dk2"],
        "J": ops["J"],
        "invJT": ops["invJT"],
        "w_k": ops["w_k"],
        "pref": pref,
        "ix_list": ops["ix_list"],
        "iy_list": ops["iy_list"],
        "band_list": band_list,
    }
    if save_intermediates:
        out.update({"E": e, "Ug": ug, "A_diag": a_diag, "r_mn": r_mn, "dr": dr, "r_cov": r_cov})
    return sigma, out


def mic_metric_plane_from_ue_mn(
    kx: np.ndarray,
    ky: np.ndarray,
    u: np.ndarray,
    e: np.ndarray,
    eph_list: np.ndarray,
    ef: float,
    k_t: float,
    eta: float,
    opts: Any | None = None,
    return_details: bool = False,
):
    periodic_fd = bool(get_opt(opts, "periodicFD", False))
    trim_boundary = bool(get_opt(opts, "trimBoundary", True))
    do_gauge_fix = bool(get_opt(opts, "doGaugeFix", True))
    g_s = float(get_opt(opts, "g_s", 1.0))
    verbose = bool(get_opt(opts, "verbose", True))
    positive_de = bool(get_opt(opts, "positiveDE", True))
    save_full_mn = bool(get_opt(opts, "saveFullMN", True))
    save_intermediates = bool(get_opt(opts, "saveIntermediates", True))

    e_charge = 1.602176634e-19
    hbar_js = 1.054571817e-34
    pref = -np.pi * g_s * (e_charge**2) / (2.0 * hbar_js)
    vel_pref = e_charge / hbar_js

    u, e, nkx, nky, nb, nb_sel = _reshape_u_e(kx, ky, u, e)

    band_list = get_opt(opts, "band_list", None)
    if band_list is not None and len(band_list) > 0:
        bidx = np.asarray(band_list, dtype=np.int64) - 1
        if np.any(bidx < 0) or np.any(bidx >= nb_sel):
            raise ValueError(f"band_list out of range, nb_sel={nb_sel}")
        u = u[:, bidx, :, :]
        e = e[:, :, bidx]
        nb_sel = bidx.size

    eph_list = np.asarray(eph_list, dtype=np.float64).reshape(-1)
    ops = skewgrid_ops_from_k(kx, ky, periodic_fd=periodic_fd, trim_boundary=trim_boundary)

    if verbose:
        print(
            f"[tbnlc.mic] Nkx={nkx} Nky={nky} nb={nb} nb_sel={nb_sel} "
            f"periodicFD={int(periodic_fd)} trimBoundary={int(trim_boundary)} positiveDE={int(positive_de)}"
        )

    ug = gauge_fix_parallel_transport(u) if do_gauge_fix else u
    du_x, du_y = fd_du_skew(ug, ops["invJT"], periodic_fd=periodic_fd)

    r_mn = np.zeros((nkx, nky, nb_sel, nb_sel, 2), dtype=np.complex128)
    for ix in ops["ix_list"]:
        for iy in ops["iy_list"]:
            u0 = ug[:, :, ix, iy]
            dux = du_x[:, :, ix, iy]
            duy = du_y[:, :, ix, iy]
            r_mn[ix, iy, :, :, 0] = 1j * (u0.conj().T @ dux)
            r_mn[ix, iy, :, :, 1] = 1j * (u0.conj().T @ duy)

    off = np.ones((nb_sel, nb_sel), dtype=np.float64) - np.eye(nb_sel, dtype=np.float64)
    r_mn_h = hermitianize_rmn(r_mn)
    r_mn_h *= off.reshape(1, 1, nb_sel, nb_sel, 1)
    r_nm_mn = np.transpose(r_mn_h, (0, 1, 3, 2, 4))

    g2 = np.zeros((nkx, nky, nb_sel, nb_sel, 2, 2), dtype=np.float64)
    for b in range(2):
        for c in range(2):
            g2[:, :, :, :, b, c] = np.real(r_mn_h[:, :, :, :, b] * r_nm_mn[:, :, :, :, c] + r_mn_h[:, :, :, :, c] * r_nm_mn[:, :, :, :, b])

    em = e.reshape(nkx, nky, nb_sel, 1)
    en = e.reshape(nkx, nky, 1, nb_sel)
    d_e_mn = em - en
    mask_pos = (d_e_mn > 0).astype(np.float64) if positive_de else np.ones_like(d_e_mn, dtype=np.float64)

    f = fermi_dirac(e, ef, k_t).astype(np.float64)
    fm = f.reshape(nkx, nky, nb_sel, 1)
    fn = f.reshape(nkx, nky, 1, nb_sel)
    f_nm = fn - fm

    d_e_dx, d_e_dy = fd_tensor_skew(e.astype(np.float64), ops["invJT"], periodic_fd=periodic_fd)
    v_band = np.zeros((nkx, nky, nb_sel, 2), dtype=np.float64)
    v_band[:, :, :, 0] = vel_pref * d_e_dx
    v_band[:, :, :, 1] = vel_pref * d_e_dy

    dv_mn = np.zeros((nkx, nky, nb_sel, nb_sel, 2), dtype=np.float64)
    for a in range(2):
        vm = v_band[:, :, :, a].reshape(nkx, nky, nb_sel, 1)
        vn = v_band[:, :, :, a].reshape(nkx, nky, 1, nb_sel)
        dv_mn[:, :, :, :, a] = vm - vn

    mask_total = ops["mask_k4"] * off.reshape(1, 1, nb_sel, nb_sel) * mask_pos

    eta_abc = mic_integrate_kernel(
        d_e_mn.astype(np.float64),
        f_nm.astype(np.float64),
        dv_mn.astype(np.float64),
        g2.astype(np.float64),
        mask_total.astype(np.float64),
        eph_list.astype(np.float64),
        float(pref),
        float(ops["w_k"]),
        float(eta),
    )

    if not return_details:
        return eta_abc

    out = {
        "pref": pref,
        "vel_pref": vel_pref,
        "w_k": ops["w_k"],
        "J": ops["J"],
        "invJT": ops["invJT"],
        "dk1": ops["dk1"],
        "dk2": ops["dk2"],
        "mask_k": ops["mask_k"],
        "mask_off": off,
        "positiveDE": positive_de,
        "band_list": band_list,
    }
    if save_intermediates:
        out["Ug"] = ug
    if save_full_mn:
        out.update({"r_mn": r_mn_h, "g2": g2, "dE_mn": d_e_mn, "f_nm": f_nm, "v_band": v_band, "dv_mn": dv_mn})
    return eta_abc, out


# MATLAB-style aliases
shift_current_plane_fd_energy_skew_fromUE = shift_current_plane_fd_energy_skew_from_ue
mic_metric_plane_fromUE_mn = mic_metric_plane_from_ue_mn
