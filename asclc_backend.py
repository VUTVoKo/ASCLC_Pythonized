"""Reusable A-SCLC functions, independent of the notebook frontend."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from asclc_model import EPS0, E_CHARGE


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


def logarithmic_slope(voltage, density):
    """Return interval midpoint V and m = delta ln(J) / delta ln(V).

    Adjacent pairs stay in acquisition order. Nonpositive/nonfinite endpoints
    and repeated voltages produce NaN; no absolute values or smoothing are used.
    Midpoints are geometric means. Results have len(voltage) - 1 entries.
    """
    v, j = np.asarray(voltage, dtype=float), np.asarray(density, dtype=float)
    if v.ndim != 1 or j.shape != v.shape or v.size < 2:
        raise ValueError("Supply matching one-dimensional arrays with at least two points.")
    valid = np.isfinite(v) & np.isfinite(j) & (v > 0) & (j > 0)
    pair = valid[:-1] & valid[1:] & (v[:-1] != v[1:])
    midpoint = np.full(v.size - 1, np.nan)
    slope = np.full(v.size - 1, np.nan)
    idx = np.flatnonzero(pair)
    lv0, lv1 = np.log(v[idx]), np.log(v[idx + 1])
    midpoint[idx] = np.exp((lv0 + lv1) / 2)
    slope[idx] = (np.log(j[idx + 1]) - np.log(j[idx])) / (lv1 - lv0)
    return midpoint, slope


def effective_mobility(voltage, density, *, thickness, eps_r):
    """A3: return V_mid, J_mid, m, gamma and mu_eff (m^2/(V s)).

    Eq. (4) uses gamma = 1/m from adjacent measured log differences.
    J is interpolated geometrically to the same geometric voltage midpoint.
    Preserve interval order; invalid pairs and slopes outside 0 < gamma < 1
    have NaN mobility. Numerically ohmic slopes are singular, not clipped.
    No smoothing or microscopic mobility inference is performed.
    """
    if not all(np.isfinite(x) and x > 0 for x in (thickness, eps_r)):
        raise ValueError("thickness and eps_r must be positive and finite.")
    v_mid, m = logarithmic_slope(voltage, density)
    j = np.asarray(density, dtype=float)
    j_mid = np.full(m.shape, np.nan)
    valid = np.isfinite(v_mid)
    j_mid[valid] = np.exp((np.log(j[:-1][valid]) + np.log(j[1:][valid])) / 2)
    gamma = np.full(m.shape, np.nan)
    np.divide(1., m, out=gamma, where=np.isfinite(m) & (m != 0))
    usable = (valid & (gamma > 0) & (gamma < 1)
              & ~np.isclose(gamma, 1., rtol=0., atol=1e-12))
    mu = np.full(m.shape, np.nan)
    g = gamma[usable]
    with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
        mu[usable] = (thickness**3 * j_mid[usable]
                      / (EPS0 * eps_r * (1-g) * (2-g)**2 * v_mid[usable]**2))
    mu[~np.isfinite(mu)] = np.nan
    return v_mid, j_mid, m, gamma, mu


def carrier_densities(voltage, density, gamma, *, thickness, eps_r, mobility):
    """A4: return p_f, p_t, p_s (m^-3), and Theta at supplied A3 midpoints.

    Use Eq. (5) for p_f and interpret Eq. (6) as p_t, as stated in SI A4.
    Then p_s = p_f + p_t and Theta = p_f/p_s (not mu_eff/mobility).
    Inputs must be matching 1-D arrays; retain all intervals and use NaN for
    nonpositive/nonfinite data or gamma outside A3's nonsingular (0, 1) range.
    mobility is the explicitly supplied microscopic mobility in m^2/(V s).
    """
    if not all(np.isfinite(x) and x > 0 for x in (thickness, eps_r, mobility)):
        raise ValueError("thickness, eps_r and mobility must be positive and finite.")
    v, j, g = (np.asarray(x, dtype=float) for x in (voltage, density, gamma))
    if v.ndim != 1 or j.shape != v.shape or g.shape != v.shape:
        raise ValueError("Supply matching one-dimensional voltage, density and gamma arrays.")
    valid = (np.isfinite(v) & np.isfinite(j) & np.isfinite(g)
             & (v > 0) & (j > 0) & (g > 0) & (g < 1)
             & ~np.isclose(g, 1., rtol=0., atol=1e-12))
    pf, pt = np.full(v.shape, np.nan), np.full(v.shape, np.nan)
    with np.errstate(over='ignore', divide='ignore', invalid='ignore', under='ignore'):
        pf[valid] = thickness*j[valid]/(E_CHARGE*mobility*(2-g[valid])*v[valid])
        pt[valid] = EPS0*eps_r*(1-g[valid])*(2-g[valid])*v[valid]/(E_CHARGE*thickness**2)
        ps = pf + pt
        theta = pf / ps
    finite = (np.isfinite(pf) & np.isfinite(pt) & np.isfinite(ps)
              & np.isfinite(theta) & (pf > 0) & (pt > 0))
    for values in (pf, pt, ps, theta):
        values[~finite] = np.nan
    return pf, pt, ps, theta
