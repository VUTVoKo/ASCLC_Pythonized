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
    """Trailing mean of values[i:i+window]. Incomplete or nonfinite windows
    return nan. Signs and row order are preserved; window=1 is identity."""
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

    Unlike ``trailing_mean`` , this keeps the
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
    """Local gamma = d ln|U| / d ln|j| from moving least-squares fits.

    Fit ln|U| on ln|j| directly, using raw readings. The default window
    starts at each row; centered=True centres it. Zero/nonfinite readings,
    incomplete windows, and constant log-current windows return nan."""
    v = np.asarray(voltage, dtype=float)
    i = np.asarray(current, dtype=float)
    if v.ndim != 1 or v.shape != i.shape:
        raise ValueError("voltage and current must be 1-D arrays of equal length.")
    window = int(window)
    if window < 2:
        raise ValueError("window must be at least 2 points.")
    with np.errstate(divide="ignore"):    # zeros become -inf, handled per window below
        x = np.log(np.abs(i))     # ln|j|
        y = np.log(np.abs(v))     # ln|U|
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


# Vacuum permittivity, F/m (CODATA 2018). The only physical constant the
# effective-mobility step needs.
EPSILON_0 = 8.8541878128e-12


def effective_mobility(voltage, current_density, gamma, *, thickness, eps_r):
    """Effective mobility from Eq. (4) of the paper (guide step A3).

        mu_eff = mu_0 * Theta = L**3 / (eps_0 eps_r (1 - g) (2 - g)**2) * j / U**2

    with ``L`` the crystal thickness, ``eps_r`` the relative permittivity and
    ``g = gamma = 1 / m`` the local logarithmic slope from ``local_gamma``.
    Returns m^2 V^-1 s^-1 for ``j`` in A/m^2, ``U`` in V and ``L`` in m; the
    grouping is fixed by the Mott-Gurney limit, where ``gamma = 1/2`` makes the
    shape factor ``(1 - g)(2 - g)**2 = 9/8`` and Eq. (4) reduces to
    ``mu = (8/9) L**3 j / (eps_0 eps_r U**2)``.

    Eq. (4) is written for the forward branch (``j``, ``U`` both positive), so
    ``|j|`` is used and the reverse branch of a symmetric device maps onto the
    same positive mobility -- the same absolute-value convention as
    ``local_gamma``. ``U**2`` is already sign-free. The sign that survives is
    the shape factor's: ``gamma > 1`` (a sub-ohmic slope, ``m < 1``) makes
    ``mu_eff`` negative, which is Eq. (4) reporting that the point lies outside
    the transport regime the equation describes, not a defect here.

    ``nan`` is returned wherever the inputs are nonfinite (including the ends
    without a full ``local_gamma`` window), where ``U = 0``, and where the shape
    factor vanishes -- ``gamma = 1`` (ohmic) and ``gamma = 2``. Near ``gamma =
    1`` the formula genuinely diverges; that flank is the rising side of the
    U-shaped curve of guide step A3, not an artefact to be clipped.

    All three arrays must be 1-D and the same length; order and length are
    preserved.
    """
    v = np.asarray(voltage, dtype=float)
    j = np.asarray(current_density, dtype=float)
    g = np.asarray(gamma, dtype=float)
    if v.ndim != 1 or j.shape != v.shape or g.shape != v.shape:
        raise ValueError("voltage, current_density and gamma must be 1-D arrays "
                         "of equal length.")
    if not np.isfinite(thickness) or thickness <= 0:
        raise ValueError("thickness must be finite and positive (m).")
    if not np.isfinite(eps_r) or eps_r <= 0:
        raise ValueError("eps_r must be finite and positive.")
    shape_factor = (1.0 - g) * (2.0 - g) ** 2
    denominator = EPSILON_0 * eps_r * shape_factor * v ** 2
    mu = np.full(v.shape, np.nan)
    usable = np.isfinite(denominator) & (denominator != 0.0) & np.isfinite(j)
    mu[usable] = thickness ** 3 * np.abs(j[usable]) / denominator[usable]
    return mu


# Elementary charge, C (CODATA 2018, exact). Equations (5) and (6) convert a
# space charge per unit volume into a carrier concentration with it.
E_CHARGE = 1.602176634e-19

THETA_MODELS = ("free_over_total", "absolute_over_total", "mobility_ratio")


def model_current_density(charge_density, free_density, *, gamma, thickness,
                          eps_r, mu_0):
    """M4 parametric voltage (V) and drift current density (A/m^2).

    U = e L^2 q / [epsilon (1-gamma)(2-gamma)]
    J = e mu_0 f (2-gamma) U / L

    Densities are in m^-3. For the field profile F(x)=F(L)(x/L)^(1-gamma),
    q is the total space-charge density at x=L, and f the free density
    there. This assumes steady single-carrier drift, uniform permittivity,
    and ideal injection F(0)=0. It is not a spatial transport solver.
    Signed densities are retained as an algebraic continuation; negative
    charge density does not describe the assumed injection branch.
    Nonfinite density pairs give nan for both outputs. No clipping or sorting.
    Gamma is supplied explicitly and must lie strictly between zero and one.
    """
    q, f = (np.asarray(a, dtype=float) for a in (charge_density, free_density))
    if q.ndim != 1 or f.shape != q.shape:
        raise ValueError("Densities must be 1-D arrays of equal length.")
    for name, value in (("thickness", thickness), ("eps_r", eps_r), ("mu_0", mu_0)):
        if not np.isscalar(value) or not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive.")
    if not np.isscalar(gamma) or not np.isfinite(gamma) or not 0 < gamma < 1:
        raise ValueError("gamma must be finite and strictly between zero and one.")
    u, j = (np.full(q.shape, np.nan) for _ in range(2))
    ok = np.isfinite(q) & np.isfinite(f)
    u[ok] = (E_CHARGE * thickness**2 * q[ok]
             / (EPSILON_0 * eps_r * (1-gamma) * (2-gamma)))
    j[ok] = E_CHARGE * mu_0 * f[ok] * (2-gamma) * u[ok] / thickness
    return u, j


def carrier_densities(voltage, current_density, gamma, *, thickness, eps_r, mu_0,
                      theta_model="free_over_total"):
    """Free and trapped carrier concentrations, Eqs. (5) and (6) (guide step A4).

        p_f = L j / (e mu_0 (2 - g) U)                              Eq. (5)
        p_6 = eps_0 eps_r (1 - g) (2 - g) U / (e L**2)              Eq. (6)

    with ``g = gamma`` from ``local_gamma`` and ``mu_0`` the microscopic
    mobility read off the step-A3 plot. Returns ``(p_f, p_t, p_s, theta)`` in
    m^-3 (``theta`` dimensionless), arrays the length of the input.

    Eq. (5) is Ohm's law in the form the SI uses for these devices (Eq. (S14),
    ``j = e mu_0 p_f (2 - g) U / L``) solved for ``p_f``. Eq. (6) follows from
    dividing Eq. (2) by Eq. (S14):

        eps_0 eps_r Theta (1 - g)(2 - g) U / L**2 = e p_f
        =>  p_f / Theta = eps_0 eps_r (1 - g)(2 - g) U / (e L**2) = p_6

    so with ``Theta = p_f / p_s`` the quantity Eq. (6) computes is the **total**
    space charge ``p_s``, and the trapped part is ``p_s - p_f``. The SI labels it
    ``p_t`` and then defines ``Theta = p_f / (p_f + p_t)`` (A4, Eq. (S12)), which
    cannot hold at the same time as Eqs. (2) and (4): reading Eq. (6) as ``p_t``
    makes ``mu_eff / mu_0 = p_f / p_t``, so exact Mott-Gurney data (``g = 1/2``,
    ``mu_eff = mu_0``) would give ``Theta = 1/2`` rather than the ``Theta = 1``
    the article states for that region.

    ``theta_model`` selects the reading; all three return the same ``p_f``:

    - ``"free_over_total"`` (default): Eq. (6) is ``p_s``, ``p_t = p_s - p_f``,
      ``Theta = p_f / p_s``. Eqs. (2), (4), (5), (6) and (S12) are then one
      system with a single ``Theta``, identically equal to ``mu_eff / mu_0``,
      equal to 1 in the Mott-Gurney limit and to 1/2 where ``p_t = p_f`` -- the
      article's own statements about both regions.
    - ``"absolute_over_total"``: the SI read literally -- ``p_t`` is Eq. (6),
      ``p_s = |p_f| + |p_t|``, ``Theta = |p_f| / p_s`` (Eq. (S12)). Absolute
      values keep the noise-floor rows, where the measured current changes sign
      inside a window, finite.
    - ``"mobility_ratio"``: ``p_t`` is Eq. (6),
      ``p_s = p_f + p_t``, ``Theta = mu_eff / mu_0``, which equals ``p_f / p_t``
      identically and so needs no mobility argument.

    Eqs. (5) and (6) are written for the forward branch, so ``|j|`` and ``|U|``
    are used -- the convention of ``effective_mobility``, under which the
    reverse branch of a symmetric device maps onto the same positive
    concentrations. The sign that survives is the shape factor's: ``1 < g < 2``
    makes Eq. (6) negative, as it makes Eq. (4)'s mobility negative, and that is
    the equation reporting a point outside the transport regime it describes.
    ``Theta > 1`` (hence a negative ``p_t`` under the default reading) means the
    row needs more mobility than the ``mu_0`` supplied; it is left in the output
    rather than clipped, because that is how ``mu_0`` gets judged.

    A row is analysed or it is not: all four outputs are ``nan`` wherever
    ``voltage``, ``current_density`` or ``gamma`` is nonfinite -- including the
    ends without a full ``local_gamma`` window -- and where ``U = 0`` or
    ``gamma = 2``, at which Eq. (5) is singular. ``theta`` is additionally
    ``nan`` where its denominator vanishes (``gamma = 1``, which leaves no space
    charge at all). All three arrays must be 1-D and the same length; order
    and length are preserved.
    """
    v = np.abs(np.asarray(voltage, dtype=float))
    j = np.abs(np.asarray(current_density, dtype=float))
    g = np.asarray(gamma, dtype=float)
    if v.ndim != 1 or j.shape != v.shape or g.shape != v.shape:
        raise ValueError("voltage, current_density and gamma must be 1-D arrays "
                         "of equal length.")
    if not np.isfinite(thickness) or thickness <= 0:
        raise ValueError("thickness must be finite and positive (m).")
    if not np.isfinite(eps_r) or eps_r <= 0:
        raise ValueError("eps_r must be finite and positive.")
    if not np.isfinite(mu_0) or mu_0 <= 0:
        raise ValueError("mu_0 must be finite and positive (m^2 V^-1 s^-1).")
    if theta_model not in THETA_MODELS:
        raise ValueError(f"Unknown theta_model: {theta_model}. Use {THETA_MODELS}.")

    usable = (np.isfinite(v) & np.isfinite(j) & np.isfinite(g)
              & (v != 0.0) & (g != 2.0))
    p_f, p_6 = (np.full(v.shape, np.nan) for _ in range(2))
    g_u, v_u = g[usable], v[usable]
    p_f[usable] = thickness * j[usable] / (E_CHARGE * mu_0 * (2.0 - g_u) * v_u)
    p_6[usable] = (EPSILON_0 * eps_r * (1.0 - g_u) * (2.0 - g_u) * v_u
                   / (E_CHARGE * thickness ** 2))

    if theta_model == "free_over_total":
        p_s = p_6
        p_t = p_s - p_f
        numerator, denominator = p_f, p_s
    elif theta_model == "absolute_over_total":
        p_t = p_6
        p_s = np.abs(p_f) + np.abs(p_t)
        numerator, denominator = np.abs(p_f), p_s
    else:                                   # mobility_ratio: theta = mu_eff/mu_0
        p_t = p_6
        p_s = p_f + p_t
        numerator, denominator = p_f, p_t
    theta = np.full(v.shape, np.nan)
    ratio = np.isfinite(numerator) & np.isfinite(denominator) & (denominator != 0.0)
    theta[ratio] = numerator[ratio] / denominator[ratio]
    return p_f, p_t, p_s, theta


def equal_density_voltage(voltage, p_f, p_t):
    """Voltages where the trapped and free concentrations cross, ``p_t = p_f``.

    Guide step A4: "if the ``p_t`` equals ``p_f`` (usually occurs at higher
    voltages), the charge starts to occupy the transport band, and the
    Mott-Gurney's law applies". Under the default reading of equation (6) that
    crossing is exactly where ``Theta = 1/2`` and so, by equation (4), where
    ``mu_eff = mu_0 / 2`` -- the test of ``mu_0`` the article proposes.

    Rows where either concentration is nonpositive or nonfinite (the noise
    floor, the ends without a full window) take no part. Between two
    consecutive surviving rows that straddle the crossing, ``ln(p_f / p_t)`` is
    interpolated linearly in ``ln|U|`` and the voltage at which it vanishes is
    returned; an exact ``p_t = p_f`` row is returned as that row's voltage.
    Returns the crossing voltages in row order, positive and in the units of
    ``voltage``; an empty array means the scan never reaches the transport band
    for the ``mu_0`` in use.
    """
    v = np.abs(np.asarray(voltage, dtype=float))
    f = np.asarray(p_f, dtype=float)
    t = np.asarray(p_t, dtype=float)
    if v.ndim != 1 or f.shape != v.shape or t.shape != v.shape:
        raise ValueError("voltage, p_f and p_t must be 1-D arrays of equal length.")
    usable = np.flatnonzero(np.isfinite(v) & (v > 0.0) & (f > 0.0) & (t > 0.0))
    if usable.size < 1:
        return np.empty(0)
    x = np.log(v[usable])
    d = np.log(f[usable]) - np.log(t[usable])     # > 0 once free carriers lead
    crossings = []
    for a, b in zip(range(usable.size - 1), range(1, usable.size)):
        if d[a] == 0.0:
            crossings.append(x[a])
        elif d[a] * d[b] < 0.0 and x[a] != x[b]:
            crossings.append(x[a] - d[a] * (x[b] - x[a]) / (d[b] - d[a]))
    if d[-1] == 0.0:
        crossings.append(x[-1])
    return np.exp(np.array(crossings, dtype=float))


# Boltzmann constant (J/K), Planck constant (J s) and the electron rest mass
# (kg), CODATA 2018; the first two are exact. Equations (S4) and (7) need them.
K_B = 1.380649e-23
H = 6.62607015e-34
M_E = 9.1093837015e-31

# The mono-energetic level of Eq. (S4) replaces the parabolic band only while
# the Fermi level stays at least this many kT from the band edge (SI,
# Supplementary Note 1).
NONDEGENERATE_KT = 3.0


def effective_dos(m_eff, temperature):
    """Effective band-state concentration, 2 (2 pi m* k_B T / h**2)**1.5.

    The integral of the parabolic DOS weighted by Boltzmann occupation gives
    this expression. m_eff is the effective mass in units of electron mass,
    temperature is in K, and the result is in m^-3, scaling as T**1.5."""
    m = M_E * np.asarray(m_eff, dtype=float)
    t = np.asarray(temperature, dtype=float)
    return 2.0 * (2.0 * np.pi * m * K_B * t / H ** 2) ** 1.5


def fermi_level(p_f, *, E_v, m_eff_h, temperature):
    """Quasi-Fermi level from the free carriers, Eq. (7) (guide step A5).

        p_f = N_v exp(-dE_F / (k_B T)),   dE_F = E_F - E_v

    solved for the Fermi level shift ``dE_F`` and returned as the level itself
    on the energy scale of ``E_v``:

        dE_F = k_B T ln(N_v / p_f)
        E_F  = E_v + dE_F

    ``p_f`` is the free-carrier concentration of step A4 (m^-3), ``E_v`` the
    valence band edge (eV, vacuum scale), ``m_eff_h`` the hole effective mass
    ratio feeding ``effective_dos``, and ``temperature`` either one value for
    the whole scan or one per point.
    Returns ``(E_F, N_v, nondegenerate)``, arrays the length of ``p_f``, with
    ``E_F`` in eV, ``N_v`` in m^-3 and ``nondegenerate`` a boolean mask.

    The shift is a **shift, not a fit**: ``E_F0``, the thermodynamic Fermi
    level at 0 V, stays a hand-selected parameter and does not enter here.
    Equation (7) fixes ``E_F`` against ``E_v`` alone, so the computed curve and
    the hand-selected ``E_F0`` are independent statements, and comparing the
    low-voltage end of one with the other is the test of ``E_F0`` the article
    makes in Fig. 5.

    Signs follow the vacuum scale of the parameter files, on which ``E_c > E_v``
    (-3.36 eV against -5.58 eV for MAPbBr3). ``p_f < N_v`` then puts ``E_F``
    above ``E_v``, inside the gap, and injecting free holes (larger ``p_f``)
    moves it back down toward the valence band -- the direction the article
    reports for MAPbBr3, its Fermi level approaching the transport band as the
    voltage rises. Degenerate cases retain the signed shift and are reported
    through ``nondegenerate``.

    ``nondegenerate`` is ``dE_F >= NONDEGENERATE_KT k_B T`` (equivalently
    ``p_f <= N_v exp(-3) ~ 0.0498 N_v``), the condition of Supplementary Note 1
    under which Eq. (S4) may stand in for the parabolic band and Boltzmann
    statistics for Fermi-Dirac. Where it is ``False``, Eq. (7) is being read
    outside its own approximation and the returned ``E_F`` understates how far
    the level has moved; the point is returned rather than dropped. The
    boundary itself, ``3 kT``, is a convention.

    A row is analysed or it is not: ``E_F`` and ``N_v`` are ``nan`` and
    ``nondegenerate`` is ``False`` wherever ``p_f`` is nonpositive or nonfinite
    -- the noise floor, the sign flips, the ends without a full ``local_gamma``
    window -- or the temperature is nonpositive or nonfinite. ``p_f`` must be
    1-D, and ``temperature`` scalar or 1-D of the same length; order and length
    are preserved.
    """
    p = np.asarray(p_f, dtype=float)
    t = np.asarray(temperature, dtype=float)
    if p.ndim != 1:
        raise ValueError("p_f must be a 1-D array.")
    if t.ndim == 0:
        t = np.full(p.shape, float(t))
    elif t.ndim != 1 or t.shape != p.shape:
        raise ValueError("temperature must be a scalar or a 1-D array as long "
                         "as p_f.")
    if not np.isfinite(E_v):
        raise ValueError("E_v must be finite (eV).")
    if not np.isfinite(m_eff_h) or m_eff_h <= 0:
        raise ValueError("m_eff_h must be finite and positive (m*/m_e).")

    usable = np.isfinite(p) & (p > 0.0) & np.isfinite(t) & (t > 0.0)
    e_f, n_v = (np.full(p.shape, np.nan) for _ in range(2))
    nondegenerate = np.zeros(p.shape, dtype=bool)
    n_v[usable] = effective_dos(m_eff_h, t[usable])
    kt = K_B * t[usable] / E_CHARGE                  # thermal energy in eV
    shift = kt * np.log(n_v[usable] / p[usable])     # dE_F = E_F - E_v
    e_f[usable] = E_v + shift
    nondegenerate[usable] = shift >= NONDEGENERATE_KT * kt
    return e_f, n_v, nondegenerate


def trap_dos(energy, *, N_t, E_t, T_t):
    """Localized trap density of states (M1).

    With a = k_B T_t in eV and u = (E-E_t)/a,
    g_t(E) = N_t/(4 a cosh(u)). Evaluate as
    N_t/(2a) exp(-abs(u))/(1+exp(-2abs(u))) to avoid overflow.

    The profile is symmetric, peaks at N_t/(4a), and integrates to
    (pi/4) N_t. N_t is a concentration scale, not the integrated state count.
    Energy and E_t are in eV, T_t in K, and the result in m^-3 eV^-1.
    N_t must be finite/nonnegative, E_t finite, and T_t finite/positive.
    Output follows the input shape; nan stays nan, infinite energy gives zero."""
    e = np.asarray(energy, dtype=float)
    if not np.isfinite(N_t) or N_t < 0:
        raise ValueError("N_t must be finite and nonnegative (m^-3).")
    if not np.isfinite(E_t):
        raise ValueError("E_t must be finite (eV).")
    if not np.isfinite(T_t) or T_t <= 0:
        raise ValueError("T_t must be finite and positive (K).")
    kt_t = K_B * T_t / E_CHARGE                  # trap width in eV
    decay = np.exp(-np.abs(e - E_t) / kt_t)      # stable decaying exponential
    return N_t / (2.0 * kt_t) * decay / (1.0 + decay * decay)


# 1E-30 keeps its logarithmic plots drawable where the density is really zero.
DOS_FLOOR = 1e-30


def valence_band_dos(energy, *, E_v, m_eff_h):
    """Parabolic valence-band DOS, C_h sqrt(E_v-E) for E < E_v.

    C_h = 4 pi (2 m_h* e)**1.5 / h**3 converts the parabolic DOS to
    m^-3 eV^-1 for energy in eV and m_eff_h in units of electron mass.
    Its Boltzmann integral is effective_dos(m_eff_h, temperature).
    DOS_FLOOR is returned outside the band, including the edge, to support
    logarithmic plots. Nonfinite energy returns nan. Mass must be positive."""
    e = np.asarray(energy, dtype=float)
    if not np.isfinite(E_v):
        raise ValueError("E_v must be finite (eV).")
    if not np.isfinite(m_eff_h) or m_eff_h <= 0:
        raise ValueError("m_eff_h must be finite and positive (m*/m_e).")
    prefactor = 4.0 * np.pi * (2.0 * M_E * m_eff_h * E_CHARGE) ** 1.5 / H ** 3
    inside = e < E_v                              # nonfinite compares False
    depth = np.where(inside, E_v - e, 0.0)        # eV below the band edge
    dos = np.where(inside, prefactor * np.sqrt(depth), DOS_FLOOR)
    return np.where(np.isfinite(e), dos, np.nan)


def conduction_band_dos(energy, *, E_c, m_eff_e):
    """Parabolic conduction-band DOS, C_e sqrt(E-E_c) for E > E_c.

    C_e = 4 pi (2 m_e* e)**1.5 / h**3; energy is in eV, m_eff_e is
    in electron-mass units, and output is in m^-3 eV^-1. DOS_FLOOR is
    returned outside the band, including the edge. Nonfinite energy returns
    nan. The effective mass must be finite and positive."""
    e = np.asarray(energy, dtype=float)
    if not np.isfinite(E_c):
        raise ValueError("E_c must be finite (eV).")
    if not np.isfinite(m_eff_e) or m_eff_e <= 0:
        raise ValueError("m_eff_e must be finite and positive (m*/m_e).")
    prefactor = 4.0 * np.pi * (2.0 * M_E * m_eff_e * E_CHARGE) ** 1.5 / H ** 3
    inside = e > E_c                              # nonfinite compares False
    depth = np.where(inside, e - E_c, 0.0)        # eV above the band edge
    dos = np.where(inside, prefactor * np.sqrt(depth), DOS_FLOOR)
    return np.where(np.isfinite(e), dos, np.nan)


def total_dos(energy, *, E_v, E_c, m_eff_h, m_eff_e, N_t, E_t, T_t,
              valence=True, conduction=True, trap=True):
    """Sum of the valence, conduction, and localized trap DOS (M1).

    Each component can be enabled independently. Disabled terms contribute
    zero. Returns m^-3 eV^-1 with the shape of energy; nonfinite energies
    return nan. model_occupations integrates occupations of these states."""
    e = np.asarray(energy, dtype=float)
    total = np.zeros(e.shape) if e.ndim else np.float64(0.0)
    if valence:
        total = total + valence_band_dos(e, E_v=E_v, m_eff_h=m_eff_h)
    if conduction:
        total = total + conduction_band_dos(e, E_c=E_c, m_eff_e=m_eff_e)
    if trap:
        total = total + trap_dos(e, N_t=N_t, E_t=E_t, T_t=T_t)
    return np.where(np.isfinite(e), total, np.nan)


@dataclass
class ModelOccupations:
    """M2 populations: s is excess total, f absolute free, t = s - f0."""

    E_F: np.ndarray
    E_F0: float
    n_s: np.ndarray
    n_f: np.ndarray
    n_t: np.ndarray
    theta_n: np.ndarray
    p_s: np.ndarray
    p_f: np.ndarray
    p_t: np.ndarray
    theta_p: np.ndarray
    n_s0: float
    n_f0: float
    p_s0: float
    p_f0: float


def model_occupations(energy, E_F, *, g_total, g_conduction, g_valence,
                      E_F0, thermal_energy):
    """M2 rectangular occupation sums on a descending energy grid.

    Supply the grid including its final interval boundary, DOS arrays in
    m^-3 eV^-1, and E_F, E_F0, thermal_energy in eV. The last DOS sample
    is not integrated. E_F is a separate finite 1-D sweep; neither grid
    is resampled.

    Each state contributes g*dE*f electrons or g*dE*(1-f) holes. Sum
    the total DOS for absolute populations A and the band DOS for free
    populations. Define s=A(E_F)-A(E_F0), t=s-free(E_F0), and
    theta=abs(free(E_F)/s). These definitions mix excess and absolute
    populations and do not form an absolute free/trapped decomposition.

    Subtract equilibrium totals after summation. Large backgrounds can
    cause cancellation. Zero denominators yield inf (or nan for 0/0).
    Negative populations and theta > 1 are retained without clipping.
    Stable logistic evaluations avoid cancellation in 1-f."""
    e = np.asarray(energy, dtype=float)
    ef = np.asarray(E_F, dtype=float)
    if (e.ndim != 1 or e.size < 2 or not np.isfinite(e).all()
            or not np.all(np.diff(e) < 0)):
        raise ValueError("energy must be a finite, strictly descending 1-D grid.")
    if ef.ndim != 1 or not np.isfinite(ef).all():
        raise ValueError("E_F must be a finite 1-D sweep.")
    if not np.isfinite(E_F0):
        raise ValueError("E_F0 must be finite (eV).")
    if not np.isfinite(thermal_energy) or thermal_energy <= 0:
        raise ValueError("thermal_energy must be finite and positive (eV).")
    components = [np.asarray(g, dtype=float)
                  for g in (g_total, g_conduction, g_valence)]
    for g in components:
        if g.shape != e.shape or not np.isfinite(g).all() or np.any(g < 0):
            raise ValueError("DOS arrays must match energy and be finite/nonnegative.")
    weights = -np.diff(e)

    def populations(level):
        x = (e[:-1] - level) / thermal_energy
        electrons = np.exp(-np.logaddexp(0.0, x))
        holes = np.exp(-np.logaddexp(0.0, -x))
        total, conduction, valence = components
        return np.array([
            np.sum(weights * total[:-1] * electrons),
            np.sum(weights * conduction[:-1] * electrons),
            np.sum(weights * total[:-1] * holes),
            np.sum(weights * valence[:-1] * holes)])

    ns0, nf0, ps0, pf0 = populations(E_F0)
    values = np.empty((4, ef.size))
    for i, level in enumerate(ef):
        values[:, i] = populations(level)
    ns, nf, ps, pf = values
    ns = ns - ns0
    ps = ps - ps0
    with np.errstate(divide="ignore", invalid="ignore"):
        qn, qp = np.abs(nf / ns), np.abs(pf / ps)
    return ModelOccupations(ef.copy(), float(E_F0), ns, nf, ns - nf0, qn,
                            ps, pf, ps - pf0, qp,
                            float(ns0), float(nf0), float(ps0), float(pf0))
