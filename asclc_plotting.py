"""Matplotlib figures for measured and calculated A-SCLC results.

Formatting is plain and publication-style (default colour cycle, light grid,
power-of-ten tick labels on logarithmic axes, frameless legends). Axis titles
and series names follow the labelling used by the charts in the reference
workbook so the figures can be compared against them panel by panel:

===================  =========================  ==================================
this module          workbook chart             labels (axis / series)
===================  =========================  ==================================
plot_measured_jv     Data-calculations, Graf 5  Voltage, U (V) / Current density,
                                                j (A/m2); "j"
plot_calculated_jv   MODEL, Graf 2              same axes; "j", "jm (pf)",
                                                "-jm (nf)", "m = 1", "m = 2"
plot_effective_...    MODEL, Graf 1              Drift mobility, mu (m2/V/s);
                                                "md", "mdm (p)", "mu0 ="
plot_carrier_dens...  MODEL, Graf 6              charge density, nf, nt (m-3);
                                                "pf", "pt", "ps"
plot_carrier_frac...  MODEL, Graf 7             Theta, Q (-); "Theta"
plot_slopes           (MODEL columns L/O)       m (-); "m", "m = 1", "m = 2"
===================  =========================  ==================================

Greek symbols the workbook draws with the Symbol font are written as Greek here
(mu, Theta); plain multi-letter names from the worksheet headers are kept as is.
Every function takes the ``Calculation`` returned by ``run_calculation`` and
returns ``(figure, axes)``.
"""
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterSciNotation
import numpy as np

_FIGSIZE = (7.2, 4.5)

# Axis titles, in the workbook's wording.
_U_AXIS = "Voltage, $U$ (V)"
_J_AXIS = "Current density, $j$ (A/m$^2$)"
_MU_AXIS = r"Drift mobility, $\mu$ (m$^2$/V/s)"
_N_AXIS = "charge density, $n_f$, $n_t$ (m$^{-3}$)"
_THETA_AXIS = r"Theta, $\Theta$ (-)"


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
    """Measured current density against voltage on logarithmic axes."""
    v, j = np.asarray(result.V, float), np.asarray(result.J, float)
    ok = _finite_positive(v, j)
    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.plot(v[ok], j[ok], "o", ms=4, label="j")
    _style(ax, xlabel=_U_AXIS, ylabel=_J_AXIS)
    return fig, ax


def plot_calculated_jv(result):
    """Measured points, the modelled hole and electron branches, and slope guides.

    The view is kept on the measured data. Model curves and guides are drawn only
    where they fall within that window, and a branch with no points there is left
    out of the legend.
    """
    (vm, jm), (vp, jp), (vn, jn), (vg, g_ohm), (_, g_sclc) = result.series
    vm, jm = np.asarray(vm, float), np.asarray(jm, float)
    ok = _finite_positive(vm, jm)
    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.plot(vm[ok], jm[ok], "o", ms=4, label="j")
    if ok.any():
        xlo, xhi = vm[ok].min() * 0.7, vm[ok].max() * 1.4
        ylo, yhi = jm[ok].min() * 0.4, jm[ok].max() * 2.5
    else:
        xlo = xhi = ylo = yhi = None
    for x, y, label, color in ((vp, jp, "jm (pf)", None),
                               (vn, jn, "-jm (nf)", None),
                               (vg, g_ohm, "m = 1", "0.35"),
                               (vg, g_sclc, "m = 2", "0.6")):
        x, y = np.asarray(x, float), np.asarray(y, float)
        m = _finite_positive(x, y)
        if xlo is not None:
            m &= (x >= xlo) & (x <= xhi) & (y >= ylo) & (y <= yhi)
        if m.any():
            ax.plot(x[m], y[m], "--" if color else "-", lw=1.2,
                    color=color, label=label)
    if xlo is not None:
        ax.set_xlim(xlo, xhi)
        ax.set_ylim(ylo, yhi)
    _style(ax, xlabel=_U_AXIS, ylabel=_J_AXIS)
    return fig, ax


def _model_drift_mobility(result):
    """Modelled drift mobility mu0*|p_f/delta_p| (workbook `mdm (p)`), clipped to
    the measured voltage range."""
    model = result.model
    v = np.asarray(model["V"], float)
    with np.errstate(divide="ignore", invalid="ignore"):
        mu = model["mobility"] * np.abs(np.asarray(model["p_f"], float)
                                        / np.asarray(model["delta_p"], float))
    v_mid = np.asarray(result.V_mid, float)
    measured = v_mid[np.isfinite(v_mid) & (v_mid > 0)]
    if measured.size:
        inside = np.isfinite(v) & (v >= measured.min()) & (v <= measured.max())
    else:
        inside = np.isfinite(v) & (v > 0)
    order = np.argsort(v[inside])
    return v[inside][order], mu[inside][order]


def plot_effective_mobility(result, *, mu0_estimate=None):
    """Measured drift mobility `md` (Eq. 4), the model `mdm (p)` curve, and the
    `mu0 =` reference line (workbook MODEL Graf 1).

    Voltage is linear; mobility is logarithmic. The reference line defaults to
    the model's microscopic mobility. Pass ``mu0_estimate`` in m^2 V^-1 s^-1 to
    override the line without changing the model.
    """
    mu0 = result.model["mobility"]
    reference = mu0 if mu0_estimate is None else float(mu0_estimate)
    if not np.isfinite(reference) or reference <= 0:
        raise ValueError("mu0_estimate must be positive and finite, or None.")
    v_mid = np.asarray(result.V_mid, float)
    mu_eff = np.asarray(result.mu_eff, float)
    ok = _finite_positive(v_mid, mu_eff)

    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.set_yscale("log")
    ax.plot(v_mid[ok], mu_eff[ok], "o", ms=4, label="md")
    model_v, model_mu = _model_drift_mobility(result)
    if model_v.size:
        ax.plot(model_v, model_mu, "-", lw=1.2, label="mdm (p)")
    ax.axhline(reference, color="0.35", ls="--", lw=1,
               label=rf"$\mu_0$ = {reference:.2e}")
    ax.set_xlim(left=0)
    _style(ax, xlabel=_U_AXIS, ylabel=_MU_AXIS)
    return fig, ax


def plot_carrier_densities(result):
    """Measured free, trapped, and total carrier densities `pf`, `pt`, `ps`
    against voltage (workbook MODEL Graf 6)."""
    v = np.asarray(result.V_mid, float)
    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.set_xscale("log")
    ax.set_yscale("log")
    for values, label in ((result.p_f, "pf (m$^{-3}$)"),
                          (result.p_t, "pt (m$^{-3}$)"),
                          (result.p_s, "ps (m$^{-3}$)")):
        y = np.asarray(values, float)
        ok = _finite_positive(v, y)
        ax.plot(v[ok], y[ok], "o", ms=4, label=label)
    _style(ax, xlabel=_U_AXIS, ylabel=_N_AXIS)
    return fig, ax


def plot_carrier_fraction(result):
    """Measured free-carrier fraction Theta = pf / (pf + pt) against voltage
    (workbook MODEL Graf 7)."""
    v = np.asarray(result.V_mid, float)
    theta = np.asarray(result.theta, float)
    ok = np.isfinite(v) & np.isfinite(theta) & (v > 0)
    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.set_xscale("log")
    ax.plot(v[ok], theta[ok], "o", ms=4, label=r"$\Theta$")
    ax.axhline(0.5, color="0.35", ls="--", lw=1, label="pf = pt")
    ax.set_ylim(0, 1)
    _style(ax, xlabel=_U_AXIS, ylabel=_THETA_AXIS)
    return fig, ax


def plot_slopes(result):
    """Measured log-log slope m = d ln j / d ln U against voltage, with the
    m = 1 (ohmic) and m = 2 (SCLC) references."""
    v = np.asarray(result.V_mid, float)
    m = np.asarray(result.m, float)
    ok = np.isfinite(v) & np.isfinite(m) & (v > 0)
    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.set_xscale("log")
    ax.plot(v[ok], m[ok], "o", ms=4, label="m")
    for level in (1.0, 2.0):
        ax.axhline(level, color="0.35", ls="--", lw=1, label=f"m = {level:.0f}")
    _style(ax, xlabel=_U_AXIS, ylabel="$m$ (-)")
    return fig, ax
