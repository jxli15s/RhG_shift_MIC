from __future__ import annotations

from typing import Any

import numpy as np

from .responses import mic_metric_plane_from_ue_mn, shift_current_plane_fd_energy_skew_from_ue
from .utils import get_opt


def compute_nonlinear_conductivity_from_ue(
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
    """Unified entry point for nonlinear conductivity from eigenpairs.

    Matches the MATLAB +tbNLC wrapper behavior:
      - `do_sigma` (default True)
      - `do_mic` (default True)
      - `sigma_opts` (dict-like options for shift current)
      - `mic_opts` (dict-like options for MIC)
    """
    do_sigma = bool(get_opt(opts, "do_sigma", True))
    do_mic = bool(get_opt(opts, "do_mic", True))
    sigma_opts = get_opt(opts, "sigma_opts", {})
    mic_opts = get_opt(opts, "mic_opts", {})

    if not (do_sigma or do_mic):
        raise ValueError("At least one of do_sigma/do_mic must be True.")

    resp: dict[str, np.ndarray] = {}
    out: dict[str, Any] = {}

    if do_sigma:
        if return_details:
            sigma_abc, out_shift = shift_current_plane_fd_energy_skew_from_ue(
                kx, ky, u, e, eph_list, ef, k_t, eta, sigma_opts, return_details=True
            )
            out["shift"] = out_shift
        else:
            sigma_abc = shift_current_plane_fd_energy_skew_from_ue(
                kx, ky, u, e, eph_list, ef, k_t, eta, sigma_opts, return_details=False
            )
        resp["sigma_abc"] = sigma_abc

    if do_mic:
        if return_details:
            eta_abc, out_mic = mic_metric_plane_from_ue_mn(
                kx, ky, u, e, eph_list, ef, k_t, eta, mic_opts, return_details=True
            )
            out["mic"] = out_mic
        else:
            eta_abc = mic_metric_plane_from_ue_mn(
                kx, ky, u, e, eph_list, ef, k_t, eta, mic_opts, return_details=False
            )
        resp["eta_abc"] = eta_abc

    if not return_details:
        return resp

    out["meta"] = {
        "do_sigma": do_sigma,
        "do_mic": do_mic,
        "num_omega": int(np.asarray(eph_list).size),
        "kgrid_size": tuple(np.asarray(kx).shape),
    }
    return resp, out


def compute_nonlinear_conductivity_fromUE(*args, **kwargs):
    """MATLAB-style alias."""
    return compute_nonlinear_conductivity_from_ue(*args, **kwargs)

