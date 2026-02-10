from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np


def _maybe_fix_v73_shape(arr: np.ndarray) -> np.ndarray:
    """Best-effort shape fix for MATLAB v7.3 (HDF5) arrays."""
    if arr.ndim <= 1:
        return arr
    # MATLAB stores column-major; h5py typically reads reversed dimension order.
    return np.transpose(arr, axes=tuple(reversed(range(arr.ndim))))


def load_hfmf_mat(path: str | Path) -> dict[str, Any]:
    """Load HFMF .mat containing at least H_int and E_int.

    Supports both legacy MAT (scipy.io.loadmat) and v7.3 MAT (h5py).
    """
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise FileNotFoundError(f"HFMF mat file not found: {p}")

    try:
        from scipy.io import loadmat

        data = loadmat(p, squeeze_me=True, struct_as_record=False)
        out: dict[str, Any] = {}
        for k in ("H_int", "E_int"):
            if k in data:
                out[k] = np.array(data[k])
        if "H_int" in out and "E_int" in out:
            return out
    except Exception:
        pass

    try:
        import h5py

        out = {}
        with h5py.File(p, "r") as f:
            for k in ("H_int", "E_int"):
                if k in f:
                    out[k] = _maybe_fix_v73_shape(np.array(f[k]))
        if "H_int" in out and "E_int" in out:
            return out
    except Exception as e:
        raise RuntimeError(f"Failed to parse MAT file {p}: {e}") from e

    raise KeyError(f"File {p} does not contain required keys H_int and E_int.")

