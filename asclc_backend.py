"""Reusable A-SCLC functions, independent of the notebook frontend."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class Measurements:
    """Raw arrays and their declared units; no unit conversion is applied."""

    U: np.ndarray
    I: np.ndarray
    T: np.ndarray
    units: dict[str, str]
    finite: np.ndarray
    path: Path


def load_measurements(path):
    """Read the standard CSV: header U(V),I(A),T(C), then three numeric columns.

    Preserve raw Celsius temperatures, signs, row order and nonfinite values.
    Blank lines are ignored; missing or extra columns are rejected.
    """
    path = Path(path)
    with path.open(encoding="utf-8-sig", newline="") as stream:
        header = stream.readline().strip()
        if header != "U(V),I(A),T(C)":
            raise ValueError("Expected CSV header U(V),I(A),T(C) in that order.")
        data = np.loadtxt(stream, delimiter=",", ndmin=2)
    if data.size == 0 or data.shape[1] != 3:
        raise ValueError("Expected measurement rows with exactly three numeric columns.")
    U, I, T = data.T.copy()
    finite = np.isfinite(U) & np.isfinite(I) & np.isfinite(T)
    return Measurements(U, I, T, {"U": "V", "I": "A", "T": "C"}, finite, path)


def current_density(current, area_m2, *, current_unit="A"):
    """Return J in A/m^2, preserving current signs and measurement order."""
    scales = {"A": 1.0, "mA": 1e-3, "uA": 1e-6, "nA": 1e-9, "pA": 1e-12}
    if current_unit not in scales:
        raise ValueError(f"Unsupported current unit: {current_unit}. Use {tuple(scales)}.")
    if not np.isfinite(area_m2) or area_m2 <= 0:
        raise ValueError("area_m2 must be finite and positive.")
    return np.asarray(current, dtype=float) * scales[current_unit] / area_m2
