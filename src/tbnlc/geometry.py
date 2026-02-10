from __future__ import annotations

import numpy as np


def rhg_lattice(a: float = 2.46) -> np.ndarray:
    """Return RhG real-space lattice matrix (3x3), same convention as MATLAB code."""
    return np.array(
        [
            [0.5, -np.sqrt(3) / 2, 0.0],
            [0.5, np.sqrt(3) / 2, 0.0],
            [0.0, 0.0, 1.0 / a],
        ],
        dtype=np.float64,
    ) * a


def reciprocal_lattice(a_mat: np.ndarray) -> np.ndarray:
    """Return reciprocal lattice matrix b = 2*pi*inv(a^T)."""
    return 2.0 * np.pi * np.linalg.inv(a_mat.T)


def get_bulk2d_kmesh(
    kxline: tuple[float, float] | list[float],
    kyline: tuple[float, float] | list[float],
    knum: int,
    b_mat: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build bulk 2D mesh, matching MATLAB get_Bulk2Dkmesh behavior.

    Returns shape (knum+1, knum+1).
    """
    x_frac = np.linspace(float(kxline[0]), float(kxline[1]), int(knum) + 1)
    y_frac = np.linspace(float(kyline[0]), float(kyline[1]), int(knum) + 1)
    xg, yg = np.meshgrid(x_frac, y_frac, indexing="xy")
    kx = xg * b_mat[0, 0] + yg * b_mat[1, 0]
    ky = xg * b_mat[0, 1] + yg * b_mat[1, 1]
    kz = xg * b_mat[0, 2] + yg * b_mat[1, 2]
    return kx, ky, kz

