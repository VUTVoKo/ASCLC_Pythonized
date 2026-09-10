"""Reusable A-SCLC functions, independent of the notebook frontend."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class Measurements:
    """Raw arrays and their declared units; no unit conversion is applied."""

    U: np.ndarray
    I: np.ndarray
    units: dict[str, str]
    finite: np.ndarray
    path: Path


def load_measurements(path):
    """Read the standard CSV: header U(V),I(A), then two numeric columns.

    Preserve signs, row order and nonfinite values.
    Blank lines are ignored; missing or extra columns are rejected.
    """
    path = Path(path)
    with path.open(encoding="utf-8-sig", newline="") as stream:
        header = stream.readline().strip()
        if header != "U(V),I(A)":
            raise ValueError("Expected CSV header U(V),I(A) in that order.")
        data = np.loadtxt(stream, delimiter=",", ndmin=2)
    if data.size == 0 or data.shape[1] != 2:
        raise ValueError("Expected measurement rows with exactly two numeric columns.")
    U, I = data.T.copy()
    finite = np.isfinite(U) & np.isfinite(I)
    return Measurements(U, I, {"U": "V", "I": "A"}, finite, path)


def current_density(current, area_m2, *, current_unit="A"):
    """Return J in A/m^2, preserving current signs and measurement order."""
    scales = {"A": 1.0, "mA": 1e-3, "uA": 1e-6, "nA": 1e-9, "pA": 1e-12}
    if current_unit not in scales:
        raise ValueError(f"Unsupported current unit: {current_unit}. Use {tuple(scales)}.")
    if not np.isfinite(area_m2) or area_m2 <= 0:
        raise ValueError("area_m2 must be finite and positive.")
    return np.asarray(current, dtype=float) * scales[current_unit] / area_m2


def trailing_mean(values, *, window):
    """Trailing moving average: point ``i`` is the mean of ``values[i : i+window]``.

    This is the averaging in the workbook's ``Data-calculations`` columns E, G
    and N (``=AVERAGE(<col><row>:<col><row+I1>)``), with ``window`` the sheet's
    ``I1 + 1``. The smoothed current density of column G is
    ``current_density(trailing_mean(I, window=w), area)``; the sheet also adds
    the voltage offset to column E and the 273.15 K offset to column N.

    The last ``window - 1`` points, and any window covering a nonfinite value,
    return ``nan``; the sheet instead shrinks the window as it runs off the end
    of the data. Signs and row order are preserved. ``window = 1`` returns the
    input unchanged (the sheet's shipped ``I1 = 0``).
    """
    x = np.asarray(values, dtype=float)
    if x.ndim != 1:
        raise ValueError("values must be a 1-D array.")
    window = int(window)
    if window < 1:
        raise ValueError("window must be at least 1 point.")
    out = np.full(x.shape, np.nan)
    for start in range(x.size - window + 1):
        block = x[start:start + window]
        if np.isfinite(block).all():
            out[start] = block.mean()
    return out


def centered_mean(values, *, window):
    """Centred moving average: point ``i`` is the mean of the window around ``i``.

    Unlike ``trailing_mean`` (the workbook's columns E/G/N), this keeps the
    smoothed series aligned with the raw points on the voltage axis, so the two
    can be plotted against each other. For an even ``window`` the extra point is
    taken from the leading (lower-index) side. The points at each end without a
    full centred window, and any window covering a nonfinite value, return
    ``nan``. ``window = 1`` returns the input unchanged.
    """
    x = np.asarray(values, dtype=float)
    if x.ndim != 1:
        raise ValueError("values must be a 1-D array.")
    window = int(window)
    if window < 1:
        raise ValueError("window must be at least 1 point.")
    lead = window // 2
    out = np.full(x.shape, np.nan)
    for i in range(lead, x.size - (window - 1 - lead)):
        block = x[i - lead:i - lead + window]
        if np.isfinite(block).all():
            out[i] = block.mean()
    return out


def local_gamma(voltage, current, *, window, centered=False):
    """Local gamma = d ln|U| / d ln|j| from a moving least-squares fit.

    This ports the workbook's ``Data-calculations`` column K (headed ``g (-)``):
    for each point it fits a straight line to (ln|j|, ln|U|) over a ``window``
    of points and returns the slope, which is gamma = 1/m with m the logarithmic
    J-V slope. ``window`` is the workbook's ``I1 + 1`` (see docs/ASCLC_guide.md,
    "Local-slope smoothing window").

    Like the sheet, it takes ln of the absolute values, so noise-floor points
    with sign flips still contribute. ``centered=False`` (default) matches the
    sheet's trailing window (anchored at ``i``); ``centered=True`` centres the
    window on ``i`` so the result stays aligned with ``voltage`` — the extra
    point of an even window goes to the leading side, as in ``centered_mean``.
    Points without a full window, and any window spanning a nonfinite ln (a zero
    or nonfinite reading), return ``nan`` rather than the sheet's silently shrunk
    window and ``0``.

    ``current`` may be raw current or current density; a constant positive
    scale factor shifts ln|j| by a constant and does not change the slope.
    Returns an array the length of ``voltage``.
    """
    v = np.asarray(voltage, dtype=float)
    i = np.asarray(current, dtype=float)
    if v.ndim != 1 or v.shape != i.shape:
        raise ValueError("voltage and current must be 1-D arrays of equal length.")
    window = int(window)
    if window < 2:
        raise ValueError("window must be at least 2 points.")
    with np.errstate(divide="ignore"):    # zeros become -inf, handled per window below
        x = np.log(np.abs(i))     # workbook column J: ln|j|
        y = np.log(np.abs(v))     # workbook column I: ln|U|
    lead = window // 2 if centered else 0
    gamma = np.full(v.shape, np.nan)
    for anchor in range(lead, v.size - (window - 1 - lead)):
        block = slice(anchor - lead, anchor - lead + window)
        xs, ys = x[block], y[block]
        if not (np.isfinite(xs).all() and np.isfinite(ys).all()):
            continue
        xc = xs - xs.mean()
        denominator = np.dot(xc, xc)
        if denominator == 0.0:
            continue
        gamma[anchor] = np.dot(xc, ys - ys.mean()) / denominator
    return gamma
