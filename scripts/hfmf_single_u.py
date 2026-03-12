#!/usr/bin/env python3
from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import numpy as np

from tbnlc import (
    calculate_ef,
    compute_nonlinear_conductivity_from_ue,
    diag_mesh_from_hk,
    extract_spinvalley_blocks,
    get_bulk2d_kmesh,
    load_hfmf_mat,
    reciprocal_lattice,
    rhg_lattice,
)


def _parse_band_list(text: str) -> list[int]:
    vals = [x.strip() for x in text.split(",") if x.strip()]
    if not vals:
        raise ValueError("band list is empty")
    out = [int(v) for v in vals]
    if any(v <= 0 for v in out):
        raise ValueError("band_list indices must be positive integers")
    return out


def _parse_blocks(text: str) -> list[str]:
    vals = [x.strip() for x in text.split(",") if x.strip()]
    if not vals:
        raise ValueError("block list is empty")
    return vals


def _normalize_h_int(h_int: np.ndarray) -> np.ndarray:
    h = np.asarray(h_int)
    if h.ndim != 4:
        raise ValueError(f"H_int must be 4D, got {h.shape}")
    if h.shape[2] == h.shape[3]:
        return h
    if h.shape[0] == h.shape[1]:
        return np.transpose(h, (2, 3, 0, 1))
    raise ValueError(f"Cannot infer H_int axis order from {h.shape}")


def _build_parser() -> argparse.ArgumentParser:
    frac = 0.15 / (2.0 * np.pi)
    p = argparse.ArgumentParser(
        description="Compute shift current + MIC for a single HFMF file (single U)."
    )
    p.add_argument("--hfmf-file", type=Path, required=True, help="Path to one MAT file with H_int/E_int")
    p.add_argument("--u-value", type=float, default=np.nan, help="Metadata only (optional)")
    p.add_argument("--band-list", type=str, default="1,2")
    p.add_argument("--blocks", type=str, default="K_up,Kp_up,K_dn,Kp_dn")
    p.add_argument("--n-layer", type=int, default=1)
    p.add_argument("--kT", type=float, default=0.0)
    p.add_argument("--eta", type=float, default=5e-4)
    p.add_argument("--eph-min", type=float, default=0.0)
    p.add_argument("--eph-max", type=float, default=0.1)
    p.add_argument("--eph-num", type=int, default=5000)
    p.add_argument("--periodic-fd", action="store_true")
    p.add_argument("--no-trim-boundary", action="store_true")
    p.add_argument("--disable-gauge-fix", action="store_true")
    p.add_argument("--disable-symbc", action="store_true")
    p.add_argument("--disable-mic", action="store_true")
    p.add_argument("--verbose-kernel", action="store_true")
    p.add_argument("--kx-min-frac", type=float, default=-frac)
    p.add_argument("--kx-max-frac", type=float, default=frac)
    p.add_argument("--ky-min-frac", type=float, default=-frac)
    p.add_argument("--ky-max-frac", type=float, default=frac)
    p.add_argument("--output", type=Path, default=None, help="Output .npz path")
    return p


def main() -> int:
    args = _build_parser().parse_args()

    if args.eph_num < 2:
        raise ValueError("--eph-num must be >= 2")

    band_list = _parse_band_list(args.band_list)
    block_names = _parse_blocks(args.blocks)
    eph_list = np.linspace(args.eph_min, args.eph_max, args.eph_num, dtype=np.float64)

    periodic_fd = bool(args.periodic_fd)
    trim_boundary = not bool(args.no_trim_boundary)
    do_gauge_fix = not bool(args.disable_gauge_fix)
    sym_bc = not bool(args.disable_symbc)
    do_mic = not bool(args.disable_mic)

    sigma_opts = {
        "band_list": band_list,
        "periodicFD": periodic_fd,
        "trimBoundary": trim_boundary,
        "symBC": sym_bc,
        "verbose": bool(args.verbose_kernel),
        "doGaugeFix": do_gauge_fix,
        "g_s": 1.0,
        "saveIntermediates": False,
    }
    mic_opts = {
        "band_list": band_list,
        "periodicFD": periodic_fd,
        "trimBoundary": trim_boundary,
        "verbose": bool(args.verbose_kernel),
        "doGaugeFix": do_gauge_fix,
        "g_s": 1.0,
        "positiveDE": True,
        "saveIntermediates": False,
        "saveFullMN": False,
    }
    run_opts = {
        "do_sigma": True,
        "do_mic": do_mic,
        "sigma_opts": sigma_opts,
        "mic_opts": mic_opts,
    }
    diag_opts = {"band_list": band_list}

    hfmf_file = args.hfmf_file.expanduser().resolve()
    if not hfmf_file.is_file():
        raise FileNotFoundError(f"HFMF file not found: {hfmf_file}")

    print(f"[single-u] load: {hfmf_file}")
    mat = load_hfmf_mat(hfmf_file)
    h_int = _normalize_h_int(mat["H_int"]).astype(np.complex128) / 1000.0
    e_int = np.asarray(mat["E_int"], dtype=np.float64) / 1000.0
    ef = calculate_ef(np.ravel(e_int), 0.5)

    nkx_data, nky_data, dim_h1, dim_h2 = h_int.shape
    if dim_h1 != dim_h2:
        raise ValueError(f"invalid H_int shape: {h_int.shape}")

    knum = nkx_data - 1
    a_mat = rhg_lattice(2.46)
    b_mat = reciprocal_lattice(a_mat)
    kxline = (args.kx_min_frac, args.kx_max_frac)
    kyline = (args.ky_min_frac, args.ky_max_frac)
    kx, ky, _ = get_bulk2d_kmesh(kxline, kyline, knum, b_mat)
    if kx.shape != (nkx_data, nky_data):
        raise ValueError(f"k-grid mismatch: Kx={kx.shape} vs H_int={(nkx_data, nky_data)}")

    h_blocks = extract_spinvalley_blocks(h_int, args.n_layer)
    sigma_sum = np.zeros((2, 2, 2, eph_list.size), dtype=np.float64)
    eta_sum = np.zeros((2, 2, 2, eph_list.size), dtype=np.float64)

    print(
        f"[single-u] blocks={block_names}, bands={band_list}, "
        f"Eph[{eph_list.size}] in [{eph_list[0]:.4f}, {eph_list[-1]:.4f}] eV"
    )
    for bname in block_names:
        if bname not in h_blocks:
            raise KeyError(f"unknown block: {bname}")
        print(f"[single-u] compute block: {bname}")
        hk = np.transpose(h_blocks[bname], (2, 3, 0, 1))
        u_blk, e_blk = diag_mesh_from_hk(hk, diag_opts)
        resp = compute_nonlinear_conductivity_from_ue(
            kx, ky, u_blk, e_blk, eph_list, ef, args.kT, args.eta, run_opts
        )
        sigma_sum += resp["sigma_abc"]
        if do_mic:
            eta_sum += resp["eta_abc"]

    out_path = args.output
    if out_path is None:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_dir = Path.cwd() / "outputs_py"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"single_u_{ts}.npz"
    else:
        out_path = out_path.expanduser().resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)

    payload = {
        "hfmf_file": str(hfmf_file),
        "U_value": float(args.u_value),
        "Kx": kx,
        "Ky": ky,
        "Eph_list": eph_list,
        "Ef": np.array([ef], dtype=np.float64),
        "sigma_abc": sigma_sum,
    }
    if do_mic:
        payload["eta_abc"] = eta_sum

    np.savez_compressed(out_path, **payload)
    print(f"[single-u] saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
