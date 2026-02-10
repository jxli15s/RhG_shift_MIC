from __future__ import annotations

import numpy as np

from tbnlc import (
    compute_nonlinear_conductivity_from_ue,
    get_bulk2d_kmesh,
    mic_metric_plane_from_ue_mn,
    reciprocal_lattice,
    rhg_lattice,
    shift_current_plane_fd_energy_skew_from_ue,
)


def _random_unitary_field(nb: int, nkx: int, nky: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    out = np.zeros((nb, nb, nkx, nky), dtype=np.complex128)
    for ix in range(nkx):
        for iy in range(nky):
            z = rng.normal(size=(nb, nb)) + 1j * rng.normal(size=(nb, nb))
            q, _ = np.linalg.qr(z)
            out[:, :, ix, iy] = q
    return out


def _toy_inputs():
    nkx = 5
    nky = 5
    nb = 2
    eph = np.linspace(0.0, 0.08, 7)

    a_mat = rhg_lattice(2.46)
    b_mat = reciprocal_lattice(a_mat)
    kx, ky, _ = get_bulk2d_kmesh((-0.02, 0.02), (-0.02, 0.02), nkx - 1, b_mat)

    u = _random_unitary_field(nb, nkx, nky, seed=7)
    x = np.linspace(-1.0, 1.0, nkx)[:, None]
    y = np.linspace(-1.0, 1.0, nky)[None, :]
    e = np.zeros((nkx, nky, nb), dtype=np.float64)
    e[:, :, 0] = -0.04 + 0.005 * (x + y)
    e[:, :, 1] = 0.04 + 0.008 * (x - y)
    return kx, ky, u, e, eph


def test_compute_api_smoke():
    kx, ky, u, e, eph = _toy_inputs()
    opts = {
        "do_sigma": True,
        "do_mic": True,
        "sigma_opts": {"band_list": [1, 2], "verbose": False, "doGaugeFix": False},
        "mic_opts": {"band_list": [1, 2], "verbose": False, "doGaugeFix": False},
    }
    resp = compute_nonlinear_conductivity_from_ue(
        kx=kx, ky=ky, u=u, e=e, eph_list=eph, ef=0.0, k_t=0.0, eta=1e-3, opts=opts
    )
    assert set(resp.keys()) == {"sigma_abc", "eta_abc"}
    assert resp["sigma_abc"].shape == (2, 2, 2, eph.size)
    assert resp["eta_abc"].shape == (2, 2, 2, eph.size)
    assert np.all(np.isfinite(resp["sigma_abc"]))
    assert np.all(np.isfinite(resp["eta_abc"]))


def test_band_list_single_band_zero_response():
    kx, ky, u, e, eph = _toy_inputs()
    sigma = shift_current_plane_fd_energy_skew_from_ue(
        kx, ky, u, e, eph, ef=0.0, k_t=0.0, eta=1e-3, opts={"band_list": [1], "verbose": False}
    )
    eta = mic_metric_plane_from_ue_mn(
        kx, ky, u, e, eph, ef=0.0, k_t=0.0, eta=1e-3, opts={"band_list": [1], "verbose": False}
    )
    assert np.allclose(sigma, 0.0)
    assert np.allclose(eta, 0.0)
