"""Figures for the A-SCLC model.

Reproduces the plots the prototype spreadsheet embeds and the article publishes,
overlaying the model branch (lines) on the analysis branch (markers). That
comparison is the method: with the five parameters chosen by hand, agreement
across every projection at once is what says the choice was right.

Requires ``matplotlib``, which is an optional extra so that importing
:mod:`asclc` stays dependency-light::

    pip install -e '.[plot]'

Every series is drawn through the ``valid`` mask. Rejected points are shown
faded rather than dropped: on a real measurement about half the points fail, and
seeing which ones is usually the fastest way to judge whether a sweep is usable.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from asclc import (
    CODATA,
    DEFAULT_TRAP_PROFILE,
    AnalysisResult,
    CarrierDensities,
    Constants,
    Device,
    EnergyGrid,
    Material,
    ModelCurve,
    ModelParams,
    TrapProfile,
    carrier_densities,
    conduction_band_dos,
    trap_dos,
    valence_band_dos,
)

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as exc:  # pragma: no cover - exercised only without the extra
    raise ImportError(
        "asclc_plot needs matplotlib, which is an optional extra. "
        "Install it with: pip install -e '.[plot]'"
    ) from exc

if TYPE_CHECKING:
    from matplotlib.figure import Figure

__all__ = [
    "FIGURES",
    "plot_bandgap_map",
    "plot_concentrations",
    "plot_fermi_level",
    "plot_jv",
    "plot_mobility",
    "plot_pt_vs_pf",
    "plot_theta",
    "write_figures",
]

#: Hole quantities red, electron blue, density of states green: the article's
#: convention in Figs. 4 and 6, kept so the two can be read side by side.
HOLE = "#c0392b"
ELECTRON = "#2471a3"
DOS = "#27ae60"
DATA = "#e67e22"

_MODEL_KW = {"lw": 1.8, "zorder": 3}
_DATA_KW = {"ls": "none", "marker": "o", "ms": 4.0, "zorder": 4}
#: Points failing `valid` are drawn in this style rather than omitted.
_REJECT_KW = {
    "ls": "none", "marker": "o", "ms": 3.0,
    "mfc": "none", "alpha": 0.35, "zorder": 2,
}


def _drawable(ax, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Points a log-or-linear axis pair can actually render.

    Positivity is required only on axes that are log-scaled. Testing it on both
    regardless drops every point of any plot whose abscissa is an energy, since
    those are negative on the vacuum scale.
    """
    keep = np.isfinite(x) & np.isfinite(y)
    if ax.get_xscale() == "log":
        keep &= x > 0.0
    if ax.get_yscale() == "log":
        keep &= y > 0.0
    return keep


def _model(ax, x, y, mask, *, color, label, ls="-"):
    """Draw a model curve, solid where valid and faded where not."""
    drawable = _drawable(ax, x, y)
    ok, bad = mask & drawable, (~mask) & drawable
    if bad.any():
        ax.plot(x[bad], y[bad], color=color, ls=ls, alpha=0.3, lw=1.2, zorder=2)
    if ok.any():
        ax.plot(x[ok], y[ok], color=color, ls=ls, label=label, **_MODEL_KW)


def _data(ax, x, y, valid, *, label, color=DATA):
    """Draw extracted points: filled where valid, hollow where rejected."""
    drawable = _drawable(ax, x, y)
    ok, bad = valid & drawable, (~valid) & drawable
    if bad.any():
        ax.plot(x[bad], y[bad], color=color, **_REJECT_KW)
    if ok.any():
        ax.plot(x[ok], y[ok], color=color, label=label, **_DATA_KW)


def _finish(ax, xlabel, ylabel, title):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10, loc="left")
    ax.grid(True, which="major", alpha=0.25, lw=0.5)
    ax.grid(True, which="minor", alpha=0.12, lw=0.4)
    handles, _ = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=8, framealpha=0.9)
    ax.figure.tight_layout()


def plot_jv(curve: ModelCurve, analysis: AnalysisResult | None = None) -> Figure:
    """Current density against voltage, log-log. The article's Fig. 2."""
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    ax.set_xscale("log")
    ax.set_yscale("log")
    _model(ax, curve.V, curve.J, curve.valid, color=HOLE, label="model $J$")
    if analysis is not None:
        _data(ax, analysis.V, np.abs(analysis.J), analysis.valid, label="measured $J$")
    _finish(ax, "Voltage $U$ (V)", "Current density $j$ (A m$^{-2}$)", "J-V")
    return fig


def plot_mobility(
    curve: ModelCurve, analysis: AnalysisResult | None = None, mu_0: float | None = None
) -> Figure:
    """Effective mobility against voltage. The article's Fig. 3."""
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    ax.set_xscale("log")
    ax.set_yscale("log")
    _model(ax, curve.V, curve.mu_eff, curve.valid, color=HOLE, label=r"model $\mu_{eff}$")
    if analysis is not None:
        _data(ax, analysis.V, analysis.mu_eff, analysis.valid, label=r"measured $\mu_{eff}$")
    if mu_0 is not None:
        ax.axhline(mu_0, color="0.35", ls="--", lw=1.0, label=rf"$\mu_0$ = {mu_0:.3g}")
    _finish(
        ax,
        "Voltage $U$ (V)",
        r"Effective mobility $\mu_{eff}$ (m$^2$ V$^{-1}$ s$^{-1}$)",
        "Effective mobility",
    )
    return fig


def plot_concentrations(
    curve: ModelCurve, analysis: AnalysisResult | None = None
) -> Figure:
    """Free and space-charge hole concentrations against voltage. Fig. 4a,b.

    The model's space charge is ``p_space_charge`` (the article's ``ptm``), not
    ``p_t``: ``p_t`` carries the equilibrium trapped population and misses the
    band contribution at high injection.
    """
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    ax.set_xscale("log")
    ax.set_yscale("log")
    _model(ax, curve.V, curve.p_space_charge, curve.valid, color=HOLE, label="$p_{tm}$")
    _model(
        ax, curve.V, curve.p_f, curve.valid, color=HOLE, ls="--", label="$p_{fm}$"
    )
    if analysis is not None:
        _data(ax, analysis.V, analysis.p_t, analysis.valid, label="$p_t$")
        _data(
            ax, analysis.V, analysis.p_f, analysis.valid, label="$p_f$", color="#8e44ad"
        )
    _finish(
        ax, "Voltage $U$ (V)", "Charge density (m$^{-3}$)", "Carrier concentrations"
    )
    return fig


def plot_theta(curve: ModelCurve, analysis: AnalysisResult | None = None) -> Figure:
    r"""Theta against voltage. Supplementary Fig. S19.

    Theta is capped at 1 by definition of :math:`\mu_{eff} = \mu_0\Theta`; points
    above it are outside the space-charge picture and clear ``valid``.
    """
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    ax.set_xscale("log")
    _model(ax, curve.V, curve.theta, curve.valid, color=HOLE, label=r"model $\Theta$")
    if analysis is not None:
        _data(ax, analysis.V, analysis.theta, analysis.valid, label=r"measured $\Theta$")
    ax.axhline(1.0, color="0.35", ls=":", lw=1.0)
    ax.set_ylim(0.0, 1.15)
    _finish(ax, "Voltage $U$ (V)", r"$\Theta$ (-)", r"$\Theta$")
    return fig


def plot_fermi_level(
    curve: ModelCurve,
    analysis: AnalysisResult | None = None,
    material: Material | None = None,
    params: ModelParams | None = None,
) -> Figure:
    """Fermi level against voltage. The article's Fig. 5."""
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    ax.set_xscale("log")
    _model(ax, curve.V, curve.E_F, curve.valid, color=HOLE, label="model $E_F$")
    if analysis is not None:
        _data(ax, analysis.V, analysis.E_F, analysis.valid, label="measured $E_F$")
    if params is not None:
        ax.axhline(
            params.E_F0, color="0.35", ls="--", lw=1.0, label=f"$E_{{F0}}$ = {params.E_F0} eV"
        )
        ax.axhline(params.E_t, color=DOS, ls=":", lw=1.0, label=f"$E_t$ = {params.E_t} eV")
    if material is not None:
        ax.axhline(material.E_v, color="0.2", lw=1.0, label=f"$E_v$ = {material.E_v} eV")
    _finish(ax, "Voltage $U$ (V)", "Energy $E_F$ (eV)", "Fermi level")
    return fig


def plot_pt_vs_pf(curve: ModelCurve, analysis: AnalysisResult | None = None) -> Figure:
    """Space charge against free charge. Supplementary Figs. S20, S22."""
    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    ax.set_xscale("log")
    ax.set_yscale("log")
    _model(
        ax, curve.p_f, curve.p_space_charge, curve.valid, color=HOLE, label="model"
    )
    if analysis is not None:
        _data(ax, analysis.p_f, analysis.p_t, analysis.valid, label="measured")
    both = np.isfinite(curve.p_f) & (curve.p_f > 0)
    both &= np.isfinite(curve.p_space_charge) & (curve.p_space_charge > 0)
    if both.any():
        lo = float(min(curve.p_f[both].min(), curve.p_space_charge[both].min()))
        hi = float(max(curve.p_f[both].max(), curve.p_space_charge[both].max()))
        ax.plot([lo, hi], [lo, hi], color="0.5", ls=":", lw=1.0, label="$p_t = p_f$")
    _finish(ax, "Free charge $p_f$ (m$^{-3}$)", "Space charge $p_t$ (m$^{-3}$)", "$p_t$ vs $p_f$")
    return fig


def plot_bandgap_map(
    curve: ModelCurve,
    params: ModelParams,
    material: Material,
    device: Device,
    analysis: AnalysisResult | None = None,
    constants: Constants = CODATA,
    *,
    profile: TrapProfile = DEFAULT_TRAP_PROFILE,
    densities: CarrierDensities | None = None,
) -> Figure:
    """Concentrations against energy, with the density of states. Figs. 4c,d and 6.

    This is the plot the method exists for: it maps where in the bandgap the
    charge sits. The hole curves come from the model branch; the electron
    counterparts need the full :func:`carrier_densities` result, which is
    computed here if not supplied.
    """
    if densities is None:
        densities = carrier_densities(
            curve.E_F, material, device, params, constants, profile=profile
        )

    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    ax.set_yscale("log")
    ones = np.ones(curve.E_F.size, dtype=bool)

    # ntm is not drawn: n_s_injected is the exact negative of p_s_injected, so on
    # a hole-injection sweep |ntm| and ptm are the same curve and the electron
    # line would simply occlude the hole one. The paper's Fig. 4c,d separates
    # them because its sweep runs both sides of E_F0.
    _model(
        ax, curve.E_F, curve.p_space_charge, ones, color=HOLE,
        label="$p_{tm}$  ($=|n_{tm}|$)",
    )
    _model(ax, curve.E_F, curve.p_f, ones, color=HOLE, ls="--", label="$p_{fm}$")
    _model(ax, curve.E_F, densities.n_f, ones, color=ELECTRON, ls="--", label="$n_{fm}$")

    if analysis is not None:
        _data(ax, analysis.E_F, analysis.p_t, analysis.valid, label="$p_t$")
        _data(
            ax, analysis.E_F, analysis.p_f, analysis.valid, label="$p_f$", color="#8e44ad"
        )

    # Mark E_t and E_F0 only when the sweep actually reaches them; E_t often sits
    # above E_F0 and so outside a hole-injection sweep, and an annotation pinned
    # to a coordinate off the axes floats loose in the margin.
    lo, hi = float(curve.E_F.min()), float(curve.E_F.max())
    for energy, colour, ls, text in (
        (params.E_t, DOS, ":", "$E_t$"),
        (params.E_F0, "0.35", "--", "$E_{F0}$"),
    ):
        if lo <= energy <= hi:
            ax.axvline(energy, color=colour, ls=ls, lw=1.0)
            ax.annotate(text, (energy, 1.0), xycoords=("data", "axes fraction"),
                        ha="center", va="bottom", fontsize=8, color=colour)

    # Density of states on its own axis: it is m^-3 eV^-1, not a concentration.
    grid = EnergyGrid.build(material, params, constants)
    E = grid.energies(material)
    inside = (E >= curve.E_F.min() - 0.05) & (E <= curve.E_F.max() + 0.05)
    E = E[inside]
    g = (
        valence_band_dos(E, material, constants)
        + conduction_band_dos(E, material, constants)
        + trap_dos(E, params, constants, profile=profile)
    )
    ax.set_xlim(float(curve.E_F.min()), float(curve.E_F.max()))
    ax2 = ax.twinx()
    ax2.set_yscale("log")
    ax2.set_xlim(ax.get_xlim())
    pos = g > 0
    ax2.plot(E[pos], g[pos], color=DOS, ls="--", lw=1.2, alpha=0.8, label="$g(E)$")
    ax2.set_ylabel("$g(E)$ (m$^{-3}$ eV$^{-1}$)", color=DOS)
    ax2.tick_params(axis="y", colors=DOS)

    _finish(ax, "Energy $E$ (eV)", "Charge density (m$^{-3}$)", "Bandgap map")
    return fig


#: Filenames written by :func:`write_figures`, in the order they are produced.
FIGURES = (
    "jv",
    "mobility",
    "concentrations",
    "theta",
    "fermi_level",
    "pt_vs_pf",
    "bandgap_map",
)


def write_figures(
    directory: str,
    curve: ModelCurve,
    params: ModelParams,
    material: Material,
    device: Device,
    analysis: AnalysisResult | None = None,
    constants: Constants = CODATA,
    *,
    profile: TrapProfile = DEFAULT_TRAP_PROFILE,
    dpi: int = 150,
    formats: tuple[str, ...] = ("png",),
) -> list[str]:
    """Write every figure into ``directory``. Returns the paths written."""
    import os

    os.makedirs(directory, exist_ok=True)
    figures = {
        "jv": lambda: plot_jv(curve, analysis),
        "mobility": lambda: plot_mobility(curve, analysis, params.mu_0),
        "concentrations": lambda: plot_concentrations(curve, analysis),
        "theta": lambda: plot_theta(curve, analysis),
        "fermi_level": lambda: plot_fermi_level(curve, analysis, material, params),
        "pt_vs_pf": lambda: plot_pt_vs_pf(curve, analysis),
        "bandgap_map": lambda: plot_bandgap_map(
            curve, params, material, device, analysis, constants, profile=profile
        ),
    }
    written: list[str] = []
    for name in FIGURES:
        fig = figures[name]()
        for ext in formats:
            path = os.path.join(directory, f"{name}.{ext}")
            fig.savefig(path, dpi=dpi)
            written.append(path)
        plt.close(fig)
    return written
