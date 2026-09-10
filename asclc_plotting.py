"""Matplotlib figure for measured A-SCLC results.

Formatting is plain and publication-style (default colour cycle, light grid,
power-of-ten tick labels on logarithmic axes, frameless legends). Axis titles
and series names use the notation of the published equations: Voltage, U (V) /
Current density, j (A/m2); series "j". Each function takes the ``Calculation``
returned by ``run_calculation`` and returns ``(figure, axes)``.
"""
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterSciNotation
import numpy as np

from asclc_backend import centered_mean

_FIGSIZE = (7.2, 4.5)

_U_AXIS = "Voltage, $U$ (V)"
_J_AXIS = "Current density, $j$ (A/m$^2$)"


def _style(ax, *, xlabel, ylabel):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, which="major", alpha=0.3)
    for name in ("x", "y"):
        if getattr(ax, f"get_{name}scale")() == "log":
            getattr(ax, f"{name}axis").set_major_formatter(LogFormatterSciNotation())
    ax.legend(frameon=False)
    ax.figure.tight_layout()


def _finite_positive(*arrays):
    """Boolean mask where every array is finite and strictly positive."""
    mask = None
    for values in arrays:
        values = np.asarray(values, dtype=float)
        ok = np.isfinite(values) & (values > 0)
        mask = ok if mask is None else mask & ok
    return mask


def plot_measured_jv(result):
    """Measured current density against voltage on logarithmic axes.

    Nonpositive readings are kept in the data but cannot appear on a
    logarithmic axis.
    """
    v, j = np.asarray(result.V, float), np.asarray(result.J, float)
    ok = _finite_positive(v, j)
    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.plot(v[ok], j[ok], "o", ms=4, label="j")
    _style(ax, xlabel=_U_AXIS, ylabel=_J_AXIS)
    return fig, ax


def plot_smoothed_jv(result, *, window):
    """Measured current density, smoothed against unsmoothed on the same axes.

    ``window`` is ``numerics["window"]``. A centred ``window``-point mean of
    ``J`` (`asclc_backend.centered_mean`) is drawn over the raw markers at the
    same voltages, so the two can be compared directly. The points at each end
    without a full window, and nonpositive readings, do not appear.
    """
    v, j = np.asarray(result.V, float), np.asarray(result.J, float)
    j_s = centered_mean(j, window=window)
    raw, smooth = _finite_positive(v, j), _finite_positive(v, j_s)
    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.plot(v[raw], j[raw], "o", ms=4, alpha=0.35, label="j")
    ax.plot(v[smooth], j_s[smooth], "-", lw=1.6,
            label=f"j, {window}-point mean")
    _style(ax, xlabel=_U_AXIS, ylabel=_J_AXIS)
    return fig, ax
