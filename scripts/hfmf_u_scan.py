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


def _parse_block_list(text: str) -> list[str]:
    out = [x.strip() for x in text.split(",") if x.strip()]
    if not out:
        raise ValueError("block list is empty")
    return out


def _normalize_h_int(h_int: np.ndarray) -> np.ndarray:
    h = np.asarray(h_int)
    if h.ndim != 4:
        raise ValueError(f"H_int must be 4D, got shape={h.shape}")
    # Expected shape for this package: (Nkx, Nky, Nh, Nh)
    if h.shape[2] == h.shape[3]:
        return h
    # Common alternative from MAT loaders: (Nh, Nh, Nkx, Nky)
    if h.shape[0] == h.shape[1]:
        return np.transpose(h, (2, 3, 0, 1))
    raise ValueError(f"Cannot infer H_int axis order from shape={h.shape}")


def _render_file_name(pattern: str, u_value: float) -> str:
    if "{" in pattern:
        return pattern.format(u=u_value)
    return pattern % u_value


def _build_parser() -> argparse.ArgumentParser:
    frac = 0.15 / (2.0 * np.pi)
    p = argparse.ArgumentParser(
        description="Scan U and compute shift current / MIC from HFMF Hamiltonian MAT files."
    )
    p.add_argument("--hfmf-data-dir", type=Path, required=True)
    p.add_argument("--file-pattern", type=str, default="ne=0.0000e12_U={u:.3f}data.mat")
    p.add_argument("--u-start", type=float, default=1.0)
    p.add_argument("--u-stop", type=float, default=31.0)
    p.add_argument("--u-step", type=float, default=1.0)
    p.add_argument("--band-list", type=str, default="1,2")
    p.add_argument("--blocks", type=str, default="K_up,Kp_up,K_dn,Kp_dn")
    p.add_argument("--n-layer", type=int, default=1)
    p.add_argument("--kT", type=float, default=0.0)
    p.add_argument("--eta", type=float, default=5e-4)
    p.add_argument("--eph-min", type=float, default=0.0)
    p.add_argument("--eph-max", type=float, default=0.1)
    p.add_argument("--eph-num", type=int, default=5000)
    p.add_argument("--kx-min-frac", type=float, default=-frac)
    p.add_argument("--kx-max-frac", type=float, default=frac)
    p.add_argument("--ky-min-frac", type=float, default=-frac)
    p.add_argument("--ky-max-frac", type=float, default=frac)
    p.add_argument("--periodic-fd", action="store_true")
    p.add_argument("--no-trim-boundary", action="store_true")
    p.add_argument("--disable-gauge-fix", action="store_true")
    p.add_argument("--disable-symbc", action="store_true")
    p.add_argument("--disable-mic", action="store_true")
    p.add_argument("--verbose-kernel", action="store_true")
    p.add_argument("--output", type=Path, default=None)
    return p


def main() -> int:
    args = _build_parser().parse_args()

    band_list = _parse_band_list(args.band_list)
    block_names = _parse_block_list(args.blocks)
    if args.eph_num < 2:
        raise ValueError("--eph-num must be >= 2")
    if args.u_step <= 0:
        raise ValueError("--u-step must be > 0")

    u_list = np.arange(args.u_start, args.u_stop + 0.5 * args.u_step, args.u_step, dtype=np.float64)
    eph_list = np.linspace(args.eph_min, args.eph_max, args.eph_num, dtype=np.float64)

    do_mic = not args.disable_mic
    periodic_fd = bool(args.periodic_fd)
    trim_boundary = not bool(args.no_trim_boundary)
    do_gauge_fix = not bool(args.disable_gauge_fix)
    sym_bc = not bool(args.disable_symbc)

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

    print(f"[scan] U points: {u_list.size}, bands: {band_list}, blocks: {block_names}")
    print(f"[scan] do_mic={int(do_mic)} periodicFD={int(periodic_fd)} trimBoundary={int(trim_boundary)}")

    hfmf_data_dir = args.hfmf_data_dir.expanduser().resolve()
    if not hfmf_data_dir.is_dir():
        raise FileNotFoundError(f"HFMF data dir not found: {hfmf_data_dir}")

    sigma_shift_u: np.ndarray | None = None
    eta_mic_u: np.ndarray | None = None
    ef_list = np.full(u_list.size, np.nan, dtype=np.float64)
    success_mask = np.zeros(u_list.size, dtype=np.bool_)
    kx = ky = None

    for iu, u_val in enumerate(u_list):
        file_name = _render_file_name(args.file_pattern, float(u_val))
        mat_path = hfmf_data_dir / file_name
        if not mat_path.is_file():
            print(f"[scan] skip U={u_val:g}: file not found: {mat_path}")
            continue

        try:
            mat = load_hfmf_mat(mat_path)
            h_int = _normalize_h_int(mat["H_int"]).astype(np.complex128) / 1000.0
            e_int = np.asarray(mat["E_int"], dtype=np.float64) / 1000.0

            nkx_data, nky_data, dim_h1, dim_h2 = h_int.shape
            if dim_h1 != dim_h2:
                raise ValueError(f"invalid H_int shape={h_int.shape}")

            if kx is None:
                knum = nkx_data - 1
                a_mat = rhg_lattice(2.46)
                b_mat = reciprocal_lattice(a_mat)
                kxline = (args.kx_min_frac, args.kx_max_frac)
                kyline = (args.ky_min_frac, args.ky_max_frac)
                kx, ky, _ = get_bulk2d_kmesh(kxline, kyline, knum, b_mat)
                if kx.shape != (nkx_data, nky_data):
                    raise ValueError(
                        f"k-grid mismatch: Kx={kx.shape}, H_int grid={(nkx_data, nky_data)}"
                    )

                sigma_shift_u = np.full((2, 2, 2, eph_list.size, u_list.size), np.nan, dtype=np.float64)
                if do_mic:
                    eta_mic_u = np.full((2, 2, 2, eph_list.size, u_list.size), np.nan, dtype=np.float64)
            else:
                if kx.shape != (nkx_data, nky_data):
                    print(f"[scan] skip U={u_val:g}: K-grid changed to {(nkx_data, nky_data)}")
                    continue

            ef = calculate_ef(np.ravel(e_int), 0.5)
            ef_list[iu] = ef

            h_blocks = extract_spinvalley_blocks(h_int, args.n_layer)
            sigma_sum = np.zeros((2, 2, 2, eph_list.size), dtype=np.float64)
            eta_sum = np.zeros((2, 2, 2, eph_list.size), dtype=np.float64)

            print(f"[scan] U={u_val:g} -> {file_name}")
            for bname in block_names:
                if bname not in h_blocks:
                    raise KeyError(f"unknown block '{bname}'")

                hk = np.transpose(h_blocks[bname], (2, 3, 0, 1))
                u_blk, e_blk = diag_mesh_from_hk(hk, diag_opts)
                resp = compute_nonlinear_conductivity_from_ue(
                    kx, ky, u_blk, e_blk, eph_list, ef, args.kT, args.eta, run_opts
                )
                sigma_sum += resp["sigma_abc"]
                if do_mic:
                    eta_sum += resp["eta_abc"]

            sigma_shift_u[:, :, :, :, iu] = sigma_sum
            if do_mic and eta_mic_u is not None:
                eta_mic_u[:, :, :, :, iu] = eta_sum
            success_mask[iu] = True

        except Exception as exc:
            print(f"[scan] skip U={u_val:g}: {exc}")
            continue

    if sigma_shift_u is None or kx is None or ky is None:
        raise RuntimeError("No valid U file was processed.")

    out_path = args.output
    if out_path is None:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_dir = Path.cwd() / "outputs_py"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"scan_U_sigma_mic_{ts}.npz"
    else:
        out_path = out_path.expanduser().resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)

    payload = {
        "sigma_shift_u": sigma_shift_u,
        "Kx": kx,
        "Ky": ky,
        "Eph_list": eph_list,
        "U_list": u_list,
        "success_mask": success_mask,
        "Ef_list": ef_list,
    }
    if do_mic and eta_mic_u is not None:
        payload["eta_mic_u"] = eta_mic_u

    np.savez_compressed(out_path, **payload)
    print(f"[scan] saved: {out_path}")
    print(f"[scan] valid U count: {int(np.count_nonzero(success_mask))}/{u_list.size}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
