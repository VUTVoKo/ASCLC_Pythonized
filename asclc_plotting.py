"""Matplotlib figure for measured A-SCLC results.

Formatting is plain and publication-style (default colour cycle, light grid,
power-of-ten tick labels on logarithmic axes, frameless legends). Axis titles
and series names use the notation of the published equations: Voltage, U (V) /
Current density, j (A/m2); series "j". Each function takes the ``Calculation``
returned by ``run_calculation`` and returns ``(figure, axes)``.
"""
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterSciNotation
from matplotlib.transforms import Bbox
import numpy as np

from asclc_backend import (centered_mean, equal_density_voltage,
                           valence_band_dos, conduction_band_dos,
                           trap_dos, total_dos, DOS_FLOOR)

_FIGSIZE = (7.2, 4.5)

_U_AXIS = "Voltage, $U$ (V)"
_J_AXIS = "Current density, $j$ (A/m$^2$)"
_MU_AXIS = r"Effective mobility, $\mu_\mathrm{eff}$ (m$^2$V$^{-1}$s$^{-1}$)"
_P_AXIS = r"Carrier concentration (m$^{-3}$)"
_THETA_AXIS = r"Fraction of free charge, $\Theta$ (-)"
_E_AXIS = r"Energy, $E$ (eV)"
_G_AXIS = r"Density of states, $g(E)$ (m$^{-3}$eV$^{-1}$)"


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


def plot_effective_mobility(result, mobility, *, mu_0=None):
    """Effective mobility from Eq. (4) against voltage on logarithmic axes.

    ``mobility`` is the ``Mobility`` returned by ``run_effective_mobility``;
    the voltage is its trailing-mean ``U``, the value each row describes.
    Linear voltage against logarithmic mobility, the axes of the published
    Fig. 3. ``mu_0`` (m^2 V^-1 s^-1), if given, is drawn as a dashed reference
    line, as in that figure: equation (4) gives ``mu_eff = mu_0 * Theta`` with
    ``Theta <= 1``, so the line is the limit the data approaches in the ohmic
    and Mott-Gurney regions. Points with no value (the ends without a full window, ``U = 0``, a
    vanishing shape factor) and points where Eq. (4) returns a negative mobility
    (``gamma > 1``) are kept in the data but cannot appear on a logarithmic
    axis.
    """
    v, mu = np.asarray(mobility.U, float), np.asarray(mobility.mu_eff, float)
    ok = _finite_positive(v, mu)
    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.set_yscale("log")
    ax.plot(v[ok], mu[ok], "o", ms=4, label=r"$\mu_\mathrm{eff}$")
    if mu_0 is not None:
        ax.axhline(mu_0, ls="--", lw=1.2, color="0.35", label=r"$\mu_0$")
    _style(ax, xlabel=_U_AXIS, ylabel=_MU_AXIS)
    return fig, ax


def plot_local_gamma(result, mobility):
    """Local slope parameter gamma against the trailing-mean voltage.

    The horizontal lines mark the reference slopes of the model: ``gamma = 1``
    (ohmic, where Eq. (4) diverges) and ``gamma = 1/2`` (Mott-Gurney).
    """
    v, g = np.asarray(mobility.U, float), np.asarray(mobility.gamma, float)
    ok = _finite_positive(v) & np.isfinite(g)
    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.set_ylim(-0.2, 1.6)
    for level, name in ((1.0, "ohmic"), (0.5, "Mott-Gurney")):
        ax.axhline(level, ls="--", lw=1.0, color="0.5")
        ax.annotate(name, (0.01, level), xycoords=("axes fraction", "data"),
                    va="bottom", fontsize="small", color="0.4")
    ax.plot(v[ok], g[ok], "o", ms=4, label=r"$\gamma$")
    _style(ax, xlabel=_U_AXIS, ylabel=r"Slope parameter, $\gamma$ (-)")
    return fig, ax


def plot_carrier_densities(carriers):
    """Free and trapped concentrations against voltage, the axes of Fig. 4a.

    ``carriers`` is the ``Carriers`` returned by ``run_carrier_densities``:
    ``p_f`` from equation (5), ``p_t`` from equation (6) under the reading named
    by ``carriers.theta_model``. Where the two cross, the charge starts to
    occupy the transport band (guide step A4); each crossing is drawn as a
    vertical line. Linear voltage against logarithmic concentration, as in
    Fig. 4a and as in ``plot_effective_mobility``. Nonpositive values -- the
    rows where equation (6) turns negative, and, under the default reading, rows
    needing more mobility than ``mu_0`` -- stay in the data but cannot appear on
    a logarithmic axis.
    """
    v = np.asarray(carriers.U, float)
    p_f, p_t = np.asarray(carriers.p_f, float), np.asarray(carriers.p_t, float)
    free, trapped = _finite_positive(v, p_f), _finite_positive(v, p_t)
    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.set_yscale("log")
    ax.plot(v[trapped], p_t[trapped], "o", ms=4, label=r"$p_\mathrm{t}$")
    ax.plot(v[free], p_f[free], "o", ms=4, label=r"$p_\mathrm{f}$")
    for crossing in equal_density_voltage(v, p_f, p_t):
        ax.axvline(crossing, ls="--", lw=1.0, color="0.5")
        ax.annotate(rf"$p_\mathrm{{t}} = p_\mathrm{{f}}$, {crossing:.2f} V",
                    (crossing, 0.02), xycoords=("data", "axes fraction"),
                    rotation=90, va="bottom", ha="right", fontsize="small",
                    color="0.4")
    _style(ax, xlabel=_U_AXIS, ylabel=_P_AXIS)
    return fig, ax


def plot_carrier_fraction(carriers, *, model=None, model_current=None):
    """Parameter theta against voltage, the axes of Supplementary Fig. S19.

    ``Theta`` is the fraction of free charge in the total, equation (S12); under
    the default reading of equation (6) it is identically ``mu_eff / mu_0``, so
    the dashed lines mark the two values the article names: ``Theta = 1`` in the
    ohmic and Mott-Gurney regions, and ``Theta = 1/2`` where ``p_t = p_f``.
    Points above 1 need more mobility than the ``mu_0`` in use; they are drawn,
    not clipped. Linear voltage against logarithmic ``Theta``; nonpositive
    values -- the rows where equation (6) turns negative -- stay in the data but
    cannot appear on a logarithmic axis.
    """
    if (model is None) != (model_current is None):
        raise ValueError('Supply both model occupations and model current.')
    if model is not None and not np.array_equal(model.E_F, model_current.E_F):
        raise ValueError('Model fractions and voltages must share the E_F sweep.')
    v, theta = np.asarray(carriers.U, float), np.asarray(carriers.theta, float)
    ok = np.isfinite(v) & _finite_positive(theta)
    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.set_yscale("log")
    references = ((1.0, r"$p_\mathrm{t} \ll p_\mathrm{f}$"),
                  (0.5, r"$p_\mathrm{t} = p_\mathrm{f}$")) if model is None else ()
    for level, name in references:
        ax.axhline(level, ls="--", lw=1.0, color="0.5")
        ax.annotate(name, (0.01, level), xycoords=("axes fraction", "data"),
                    va="bottom", fontsize="small", color="0.4")
    ax.plot(np.where(ok, v, np.nan), np.where(ok, theta, np.nan),
            "o", ms=4, label=r"$\Theta$")
    if model is not None:
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        valid = np.isfinite(model_current.U_p) & _finite_positive(model.theta_p)
        ax.plot(np.where(valid, model_current.U_p, np.nan),
                np.where(valid, model.theta_p, np.nan), label=r'$\Theta_m$')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
    _style(ax, xlabel=_U_AXIS, ylabel=r'Theta, $\Theta$ (-)')
    return fig, ax


def plot_model_occupations(model):
    """Free and trapped electrons and holes versus EF on one chart.

    Nonpositive/nonfinite values leave gaps on the logarithmic axis;
    signed populations remain available in model. Returns (figure, axes).
    """
    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.set_yscale('log')
    for carrier, suffix, colour, style in (
            ('p', 't', 'C0', '-'), ('n', 't', 'C1', '-'),
            ('p', 'f', 'C0', '--'), ('n', 'f', 'C1', '--')):
        values = getattr(model, f'{carrier}_{suffix}')
        ax.plot(model.E_F, np.where(
            np.isfinite(values) & (values > 0), values, np.nan),
            color=colour, ls=style, label=rf'${carrier}_{{{suffix}}}$')
    ax.axvline(model.E_F0, color='0.5', ls=':', label=r'$E_{F0}$')
    _style(ax, xlabel=r'Quasi-Fermi level, $E_F$ (eV)', ylabel=_P_AXIS)
    fig.tight_layout()
    return fig, ax


def plot_trapped_charge(carriers, model):
    """Compare trapped and free populations on two independent axis pairs.

    Primary coordinates are (p_f, p_t), (p_f, p_f), and their model
    counterparts, including model (p_f, n_t). The secondary pair is
    (n_f, p_f), with a reversed linear horizontal axis and its own log Y.
    The retained top label is Energy, E_F, although its data are n_f.
    No energy conversion is applied. Nonpositive log values leave gaps.
    Limits follow measured p_t and p_f only; the top axis follows model
    rows within that measured free-charge window.
    """
    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.set_xscale('log')
    ax.set_yscale('log')
    for x, y, label, style in (
            (carriers.p_f, carriers.p_t, r'$p_t$ (m$^{-3}$)', 'o'),
            (carriers.p_f, carriers.p_f, r'$p_f$ (m$^{-3}$)', 'o'),
            (model.p_f, model.p_t, r'$p_{tm}$', '-'),
            (model.p_f, model.p_f, r'$p_{fm}$', '-'),
            (model.p_f, model.n_t, r'$n_{tm}$', '-')):
        valid = _finite_positive(x, y)
        ax.plot(np.where(valid, x, np.nan), np.where(valid, y, np.nan),
                style, ms=4, label=label)
        if len(ax.lines) == 2:
            # Capture measured bounds before adding the unbounded model curves.
            measured_xlim, measured_ylim = ax.get_xlim(), ax.get_ylim()
    ax.set_xlim(measured_xlim)
    ax.set_ylim(measured_ylim)

    secondary = fig.add_subplot(111, label='free_population_pair', frameon=False)
    secondary.set_yscale('log')
    secondary.xaxis.tick_top()
    secondary.xaxis.set_label_position('top')
    secondary.yaxis.set_visible(False)
    secondary.set_xlabel(r'Energy, $E_F$ (eV)')
    valid = np.isfinite(model.n_f) & _finite_positive(model.p_f)
    secondary.plot(np.where(valid, model.n_f, np.nan),
                   np.where(valid, model.p_f, np.nan),
                   color='C5', label=r'$p_{fm}$')
    # Fit the secondary coordinates to the same free-charge window, while
    # retaining the complete model line for clipping at the display bounds.
    window = valid & (model.p_f >= measured_xlim[0]) & (model.p_f <= measured_xlim[1])
    secondary.dataLim = Bbox.null()
    if window.any():
        secondary.update_datalim(np.column_stack((model.n_f[window], model.p_f[window])))
    secondary.autoscale_view()
    secondary.set_ylim(measured_ylim)
    secondary.invert_xaxis()
    _style(ax, xlabel=r'Free charge, $n_f$ (m$^{-3}$)',
           ylabel=r'Trapped charge, $n_t$ (m$^{-3}$)')
    handles, labels = ax.get_legend_handles_labels()
    extra_handles, extra_labels = secondary.get_legend_handles_labels()
    ax.legend(handles + extra_handles, labels + extra_labels, frameon=False)
    return fig, (ax, secondary)


def plot_model_fermi_level(model, fermi):
    """M3 hole energy separation versus the M2 absolute free-hole density.

    Preserve the full signed energy range. Nonpositive/nonfinite densities
    leave gaps on the logarithmic concentration axis.
    """
    if not np.array_equal(model.E_F, fermi.E_F):
        raise ValueError("Occupations and energy differences must share the E_F sweep.")
    pf = np.asarray(model.p_f, dtype=float)
    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.set_xscale('log')
    ax.plot(np.where(_finite_positive(pf), pf, np.nan), fermi.delta_E_v,
            label=r'$E_F-E_v$')
    ax.axhline(0, color='0.5', ls=':', label=r'$E_F=E_v$')
    _style(ax, xlabel=r'Free-hole concentration, $p_f$ (m$^{-3}$)',
           ylabel=r'Fermi level separation, $E_F-E_v$ (eV)')
    return fig, ax


def plot_fermi_energy(measured, model_current, model_levels):
    """Absolute measured/model Fermi energies versus hole-branch voltage.

    Linear axes fit the measured points and equilibrium reference; the full
    model curve remains available but does not widen the display limits.
    """
    if not np.array_equal(model_current.E_F, model_levels.E_F):
        raise ValueError('Model voltage and energy must share the E_F sweep.')
    fig, ax = plt.subplots(figsize=_FIGSIZE)
    valid = np.isfinite(measured.U) & np.isfinite(measured.E_F)
    ax.plot(np.where(valid, measured.U, np.nan),
            np.where(valid, measured.E_F, np.nan), 'o', ms=4,
            label=r'$E_F$')
    ax.axhline(model_levels.E_F0, color='0.5', ls='--', label=r'$E_{F0}$')
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    valid = np.isfinite(model_current.U_p) & np.isfinite(model_levels.E_F)
    ax.plot(np.where(valid, model_current.U_p, np.nan),
            np.where(valid, model_levels.E_F, np.nan), label=r'$E_{Fm}$')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    _style(ax, xlabel=_U_AXIS, ylabel=r'Energy, $E_F$ (eV)')
    return fig, ax


def plot_density_energy(energy, dos, carriers, measured_levels, model, derivatives):
    """DOS, density derivatives and carrier populations versus absolute energy.

    Keep the declared series-to-Y-axis grouping; both Y axes are logarithmic.
    The common linear energy window follows measured energies. Primary Y
    bounds follow measured derivatives, secondary Y measured holes.
    """
    fig, ax = plt.subplots(figsize=(9, 5.5))
    right = ax.twinx()
    ax.set_yscale('log')
    right.set_yscale('log')
    e_meas = measured_levels.E_F
    series = [
        (ax, energy, dos, r'$g(E)$', '-'),
        (ax, e_meas, derivatives.measured, r'$dp/dE_F$', 'o'),
        (right, model.E_F, model.n_t, r'$n_{tm}$', '-'),
        (ax, model.E_F, derivatives.model, r'$dp_m/dE_F$', '-'),
        (ax, model.E_F, model.p_t, r'$p_{tm}$', '-'),
        (ax, model.E_F, model.p_f, r'$p_{fm}$', '-'),
        (ax, model.E_F, model.n_f, r'$n_{fm}$', '-'),
        (right, e_meas, carriers.p_t, r'$p_t$ (m$^{-3}$)', 'o'),
        (right, e_meas, carriers.p_f, r'$p_f$ (m$^{-3}$)', 'o')]
    handles = []
    for i, (target, x, y, label, style) in enumerate(series):
        x, y = np.asarray(x, float), np.asarray(y, float)
        if x.shape != y.shape:
            raise ValueError('Each energy series must match its density array.')
        valid = np.isfinite(x) & _finite_positive(y)
        # Keep the positive DOS floor as a below-axis endpoint, so the band
        # edges descend out of view instead of stopping at the last sample.
        line, = target.plot(np.where(valid, x, np.nan),
                            np.where(valid, y, np.nan), style,
                            color=f'C{i}', ms=3, label=label)
        handles.append(line)
    finite = e_meas[np.isfinite(e_meas)]
    if not finite.size:
        raise ValueError('No finite measured energies to determine display limits.')
    span = np.ptp(finite)
    pad = .05 * (span if span > 0 else max(abs(finite[0]), 1.))
    lo, hi = finite.min() - pad, finite.max() + pad
    # Rebuild bounds only from relevant data, without shortening plotted arrays.
    for target, indices in ((ax, (1,)), (right, (7, 8))):
        target.dataLim = Bbox.null()
        for i in indices:
            x, y = handles[i].get_data()
            valid = np.isfinite(x) & _finite_positive(y) & (x >= lo) & (x <= hi)
            if valid.any():
                target.update_datalim(np.column_stack((x[valid], y[valid])))
        target.autoscale_view(scalex=False, scaley=True)
        target.yaxis.set_major_formatter(LogFormatterSciNotation())
    ax.set_xlim(lo, hi)
    ax.set_xlabel(r'Energy, $E$, $\Delta E_F$ (eV)')
    ax.set_ylabel(r'$g(E),\ dn/dE_F$ (m$^{-3}$eV$^{-1}$)')
    right.set_ylabel(r'$n_f,\ n_t$ (m$^{-3}$)')
    ax.grid(True, which='major', alpha=.3)
    ax.legend(handles=handles, frameon=False, loc='upper left',
              bbox_to_anchor=(1.16, 1.))
    fig.tight_layout()
    return fig, (ax, right)


def plot_model_current(model, *, measured=None):
    """Separate M4 injection curves; invalid log points leave gaps."""
    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.set_xscale("log")
    ax.set_yscale("log")
    if measured is not None:
        ok = _finite_positive(measured.V, measured.J)
        ax.plot(measured.V[ok], measured.J[ok], "o", ms=4, label="Measured")
    for u, j, label in ((model.U_n, model.J_n, "Electrons"),
                        (model.U_p, model.J_p, "Holes")):
        ok = _finite_positive(u, j)
        ax.plot(np.where(ok, u, np.nan), np.where(ok, j, np.nan), label=label)
    _style(ax, xlabel=_U_AXIS, ylabel=_J_AXIS)
    return fig, ax


def plot_model_mobility(mobility):
    """M5 electron and hole effective mobilities on the existing energy sweep."""
    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.set_yscale('log')
    for values, label, colour in (
            (mobility.mu_p, r'$\mu_{\mathrm{eff},p}$', 'C0'),
            (mobility.mu_n, r'$\mu_{\mathrm{eff},n}$', 'C1')):
        ax.plot(mobility.E_F, np.where(_finite_positive(values), values, np.nan),
                label=label, color=colour)
    ax.axhline(mobility.mu_0, color='0.5', ls=':', label=r'$\mu_0$')
    _style(ax, xlabel=r'Quasi-Fermi level, $E_F$ (eV)', ylabel=_MU_AXIS)
    return fig, ax


def plot_dos(energy, *, E_v, E_c, m_eff_h, m_eff_e, N_t, E_t, T_t):
    """Plot the two parabolic bands, localized trap DOS, and their sum.

    Energy is in eV on the vacuum scale. Mark the band edges and trap
    maximum. DOS_FLOOR and underflowed zeros are omitted on logarithmic
    axes. All numerical parameters are passed through to the backend.
    Returns (figure, axes)."""
    e = np.asarray(energy, dtype=float)
    bands = dict(E_v=E_v, E_c=E_c, m_eff_h=m_eff_h, m_eff_e=m_eff_e)
    traps = dict(N_t=N_t, E_t=E_t, T_t=T_t)
    total = total_dos(e, **bands, **traps)
    # Fixed colours: a band absent from this energy range must not repaint
    # the series that remain.
    series = ((valence_band_dos(e, E_v=E_v, m_eff_h=m_eff_h),
               r"$g_\mathrm{h}$, valence band", "C0"),
              (conduction_band_dos(e, E_c=E_c, m_eff_e=m_eff_e),
               r"$g_\mathrm{e}$, conduction band", "C1"),
              (trap_dos(e, **traps), r"$g_\mathrm{t}$, trap", "C2"))
    fig, ax = plt.subplots(figsize=_FIGSIZE)
    ax.set_yscale("log")
    ax.plot(e, np.where(total > DOS_FLOOR * 2, total, np.nan), lw=4.0, alpha=0.22,
            color="0.35", solid_capstyle="round", label=r"$g$, total")
    for values, label, colour in series:
        drawn = np.where(values > DOS_FLOOR, values, np.nan)
        if not np.isfinite(drawn).any():      # wholly outside this energy range
            continue
        ax.plot(e, drawn, lw=1.7, color=colour, label=label)
    low, high = np.nanmin(e), np.nanmax(e)
    for level, name in ((E_v, "$E_v$"), (E_t, "$E_t$"), (E_c, "$E_c$")):
        if not low <= level <= high:      # off the plotted range, so not marked
            continue
        ax.axvline(level, ls="--", lw=1.0, color="0.5")
        ax.annotate(f"{name} = {level:g} eV", (level, 0.02),
                    xycoords=("data", "axes fraction"), rotation=90,
                    va="bottom", ha="right", fontsize="small", color="0.4")
    ax.set_xlim(low, high)                    # the edge markers must not widen it
    top = np.nanmax(total)
    if np.isfinite(top) and top > 0:      # else every row is floor or nan
        ax.set_ylim(top * 1e-21, top * 10)
    _style(ax, xlabel=_E_AXIS, ylabel=_G_AXIS)
    return fig, ax
