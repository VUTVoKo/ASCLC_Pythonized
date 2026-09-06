"""Advanced space-charge-limited current (A-SCLC) model.

A Python implementation of the A-SCLC model of Gavranovic, Zmeskal, Weiter and
Pospisil, *Communications Physics* **8**, 280 (2025),
doi:10.1038/s42005-025-02202-1, ported from the reference Excel prototype.

The model has two independent branches that meet only when the results are
plotted together:

**Analysis branch** (:func:`analyse_jv`, steps A1-A7 of Supplementary Note 2)
    Takes a measured current-voltage curve. The local logarithmic slope
    ``m = d ln J / d ln V`` gives ``gamma = 1/m``, and the effective mobility and
    the free and trapped carrier concentrations follow algebraically from
    Eqs (4)-(7). Voltage is the independent variable.

**Model branch** (:func:`model_curve`, steps M1-M5)
    Takes a set of five material parameters and sweeps the Fermi level. For each
    ``E_F`` the occupation integrals give the free and trapped concentrations,
    the voltage follows by inverting Eq (6) and the current from Eq (S14).
    **The Fermi level is the independent variable, not the voltage.** No
    self-consistent solve is involved anywhere.

Fitting the five parameters is done by hand at present. Every model entry point
is a pure function of a :class:`ModelParams` instance, so wrapping an optimiser
around one later requires no change to this module.

Units
-----
Energies are in eV throughout. Concentrations are m^-3, mobilities m^2 V^-1 s^-1,
lengths m, current densities A m^-2. Temperatures K.

Unresolved questions and known limitations
------------------------------------------
Points where the prototype, the article and the SI disagree are resolved by an
explicit, switchable default rather than silently, so that changing a decision is
a keyword argument and not a rewrite. See :class:`TrapProfile`, :class:`GammaModel`
and :class:`SpaceCharge`; ``ASCLC_spec.md`` section 6 lists what is still open and
how each could be settled from the data.

Validation against the prototype spreadsheet: with the default settings the
occupation integrals reproduce the prototype's own ``n(E)!L`` and ``n(E)!O``
columns to within 0.1-2 % over five orders of magnitude. The residual is the
prototype's rectangle rule on a 3 meV grid across the square-root band edge,
where this module uses a substitution that removes the singularity; the values
here are the more accurate ones. See ``test_asclc.py``.

Sweeping ``E_F`` to the band edge drives the injected charge, and the voltage
with it, to meaningless values. Pass ``V_max`` to :func:`model_curve` to bound
the sweep at the top of the measured range.
"""

from __future__ import annotations

import argparse
import enum
import sys
import warnings
from dataclasses import dataclass, fields, replace
from typing import Final

import numpy as np
import numpy.typing as npt
from scipy.special import expit

__all__ = [
    "MAPBBR3",
    "MAPBBR3_S2",
    "MAPBI3",
    "AnalysisResult",
    "CarrierDensities",
    "Configuration",
    "Constants",
    "Device",
    "EnergyGrid",
    "GammaModel",
    "Material",
    "Measurement",
    "ModelCurve",
    "ModelParams",
    "SpaceCharge",
    "ThetaModel",
    "TrapProfile",
    "analyse_jv",
    "carrier_densities",
    "conduction_band_dos",
    "effective_dos",
    "gamma_from_trap_temperature",
    "load_config",
    "load_jv",
    "local_loglog_slope",
    "model_curve",
    "trap_dos",
    "valence_band_dos",
]

FloatArray = npt.NDArray[np.float64]


# --------------------------------------------------------------------------- #
# Physical constants
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Constants:
    """Physical constants.

    Defaults are CODATA 2018. The prototype spreadsheet uses slightly rounded
    values; :meth:`prototype` reproduces them, which matters only when comparing
    against the prototype digit for digit. The largest relative difference is
    in ``k_B`` at 4.7e-4, which propagates to roughly the same relative shift in
    every Boltzmann factor.

    Attributes
    ----------
    e : float
        Elementary charge, C.
    k_B : float
        Boltzmann constant, J K^-1.
    h : float
        Planck constant, J s.
    eps_0 : float
        Vacuum permittivity, F m^-1.
    m_0 : float
        Electron rest mass, kg.
    """

    e: float = 1.602176634e-19
    k_B: float = 1.380649e-23
    h: float = 6.62607015e-34
    eps_0: float = 8.8541878128e-12
    m_0: float = 9.1093837015e-31

    @classmethod
    def prototype(cls) -> Constants:
        """Return the rounded constants hard-coded in the prototype spreadsheet."""
        return cls(
            e=1.602e-19,
            k_B=1.38e-23,
            h=6.62607004e-34,
            eps_0=8.854e-12,
            m_0=9.109e-31,
        )

    def kT_eV(self, temperature: float) -> float:
        """Thermal energy in eV at ``temperature`` in K."""
        return self.k_B * temperature / self.e


CODATA: Final[Constants] = Constants()


# --------------------------------------------------------------------------- #
# Modelling choices that the sources do not settle
# --------------------------------------------------------------------------- #


class TrapProfile(enum.Enum):
    """Functional form of the localized-state distribution.

    ``SI_BIEXPONENTIAL`` is Eq (S5) as printed, ``N_t/(k_B T_t) e^u/(1+e^u)^2``
    with ``u = (E - E_t)/(k_B T_t)``. It integrates to exactly ``N_t``.

    ``PROTOTYPE_SECH`` is what the prototype computes. Its cell
    formula reads ``.../(1 + EXP(u)^2)``; since ``^`` binds before ``+`` in
    Excel this is ``1 + e^{2u}`` rather than ``(1 + e^u)^2``, giving a
    ``sech(u)/2`` profile. Combined with the prototype's ``N_t/(2 k_B T_t)``
    prefactor it integrates to ``(pi/4) N_t``, not ``N_t``. Whether that is an
    error or a deliberate departure is unresolved. It is retained so pre-port
    results stay reproducible; ``SI_BIEXPONENTIAL`` is the default, being both
    what the SI documents and the only variant normalized to ``N_t``. See
    ``docs/prototype_differences.md``.

    ``GAUSSIAN`` is Eq (S6), with ``sigma = 2 k_B T_t``. Normalized to ``N_t``.
    Not used by the prototype. Note that Eq (S6) is printed with the exponent
    positive, ``exp(+(E - E_t)^2 / 2 sigma^2)``, which diverges; the sign is
    taken as negative here.
    """

    SI_BIEXPONENTIAL = "si_biexponential"
    PROTOTYPE_SECH = "prototype_sech"
    GAUSSIAN = "gaussian"


class SpaceCharge(enum.Enum):
    """Which charge drives the voltage when Eq (6) is inverted.

    ``INJECTED_TOTAL`` uses ``p_s(E_F) - p_s(E_F0)``, the total injected charge,
    free plus trapped, referenced to equilibrium. This is the default and it is
    what the prototype spreadsheet computes: it reproduces the prototype's own
    ``n(E)!O`` column to within 0.1-1.5 % across five orders of magnitude, the
    residual being the prototype's coarse rectangle rule on the square-root band
    edge. It vanishes at ``E_F = E_F0`` as zero bias requires.

    ``INJECTED_TRAPPED`` uses ``p_t(E_F) - p_t(E_F0)``, trapped charge only. It
    agrees with the above while traps dominate but falls away once the valence
    band contributes, reaching 140 times too small by ``E_F = E_v + 0.38 eV``,
    where the modelled voltage saturates near 1.7 V. For comparison only.

    ``TOTAL`` uses ``p_t(E_F)`` with no equilibrium reference, so the voltage does
    not vanish at zero bias: the trap is already 68 % hole-occupied at ``E_F0``
    for these parameters. For comparison only.

    Notes
    -----
    Eq (S11) writes the total concentration as a single integral whose lower limit
    is ``E_F0``, which acts as an equilibrium reference rather than a literal
    range of integration.
    """

    INJECTED_TOTAL = "injected_total"
    INJECTED_TRAPPED = "injected_trapped"
    TOTAL = "total"


class GammaModel(enum.Enum):
    """Choice of ``gamma`` for the model branch.

    ``TT_OVER_T`` is ``T_t / T``, which is what the prototype uses
    (cell ``j(U)!C6 = MODEL!B30/MODEL!B9``). This is the default, so the port
    reproduces the prototype's curves.

    ``SI_NOTE_M4`` is Supplementary Note M4 as printed:

    .. math::
        \\gamma = \\frac{T}{T_t + T} \\;(T_t \\ge T), \\qquad
        \\gamma = 0.5 \\;(T_t < T)

    Note that ``T/(T_t + T) = 1/(1 + T_t/T)``, which is ``1/m`` for the
    Mark-Helfrich trap-filled-limit exponent ``m = 1 + T_t/T`` of an exponential
    trap distribution. It is therefore the reading that squares Note M4 with the
    article's own definition ``gamma = 1/m``. For a cold trap (``T_t < T``, which
    includes every configuration in the article and the prototype) it selects the
    constant 0.5 branch, so it does *not* reproduce the prototype.

    ``TT_OVER_T_PLUS_TT`` is ``T_t / (T + T_t)``. **This appears in none of the
    sources.** It is the complement ``1 - T/(T_t + T)`` of the Note M4
    expression, and earlier revisions of this port carried it mislabelled as
    "the SI's" reading. It is retained only so that comparisons made against
    those revisions remain reproducible; prefer ``SI_NOTE_M4`` or ``TT_OVER_T``.

    At ``T_t = 30 K`` and ``T = 299 K`` the three give 0.1003, 0.5 and 0.0912.
    The choice does not affect the shape of a modelled J-V curve; see
    ``ASCLC_spec.md`` section 7.2.
    """

    TT_OVER_T = "Tt/T"
    SI_NOTE_M4 = "SI-M4"
    TT_OVER_T_PLUS_TT = "Tt/(T+Tt)"


class ThetaModel(enum.Enum):
    """Which of the three readings of Eq (S12) defines ``Theta``.

    The readings agree at high injection and differ by up to a factor of six
    through the trap-filling region, so this is a modelling choice rather than a
    detail. See ``ASCLC_spec.md`` section 7.3.

    ``ABSOLUTE_OVER_INJECTED`` is ``|p_f / p_s_injected|``: the absolute free
    concentration over the injected total. This is what the prototype computes
    (``n(E)!N = ABS(L/J)``, where ``L`` is the absolute band integral and ``J``
    the equilibrium-referenced total), and it keeps the two branches consistent,
    since the analysis branch's ``Theta = mu_eff/mu_0 = p_f/p_t`` likewise puts
    an absolute ``p_f`` from Eq (5) over the Eq (6) space charge. It is the
    default. **It is not bounded by 1**: where the equilibrium free population
    still rivals the injected charge it exceeds 1, which clears ``valid``.

    ``ABSOLUTE_OVER_TOTAL`` is ``p_f / p_s = p_f / (p_f + p_t)``, Eq (S12) as
    printed in both the SI and the article. Bounded by 1.

    ``INJECTED_OVER_INJECTED`` is ``|dp_f / p_s_injected|``, referencing
    numerator and denominator to equilibrium alike. Bounded by 1. Used by
    neither source, but it is the reading under which ``Theta`` is a fraction of
    one consistently defined population.
    """

    ABSOLUTE_OVER_INJECTED = "absolute_over_injected"
    ABSOLUTE_OVER_TOTAL = "absolute_over_total"
    INJECTED_OVER_INJECTED = "injected_over_injected"


#: Default trap profile: Eq (S5) as printed, the only variant whose integral is
#: ``N_t``, which is what makes ``N_t`` a concentration rather than a scale
#: factor. :attr:`TrapProfile.PROTOTYPE_SECH` reproduces the prototype
#: spreadsheet instead; see ``docs/prototype_differences.md``.
DEFAULT_TRAP_PROFILE: Final[TrapProfile] = TrapProfile.SI_BIEXPONENTIAL

#: Default gamma model: Supplementary Note M4 as printed. For every published
#: configuration ``T_t < T``, so this selects the constant branch, gamma = 0.5.
#: gamma is not identifiable from J-V data and only rescales the axes; see
#: ``ASCLC_spec.md`` section 7.2 before reading anything into the value.
DEFAULT_GAMMA_MODEL: Final[GammaModel] = GammaModel.SI_NOTE_M4

#: Default space-charge reference. See :class:`SpaceCharge`.
DEFAULT_SPACE_CHARGE: Final[SpaceCharge] = SpaceCharge.INJECTED_TOTAL

#: Default reading of Eq (S12) as printed, ``p_f / (p_f + p_t)``. Bounded by 1
#: by construction, so ``valid`` never fires on Theta under it. See
#: :class:`ThetaModel` and ``ASCLC_spec.md`` section 7.3.
DEFAULT_THETA_MODEL: Final[ThetaModel] = ThetaModel.ABSOLUTE_OVER_TOTAL

#: Default window, in points, for the local log-log slope in the analysis
#: branch. The prototype's binning cell is set to 0, which degenerates to a
#: single-point regression and returns gamma = 0 for every row; no such silent
#: default is offered here.
DEFAULT_SLOPE_WINDOW: Final[int] = 7


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Material:
    """Fixed material properties.

    Attributes
    ----------
    name : str
        Label used in plots and reports.
    E_c, E_v : float
        Conduction and valence band edges, eV, on the vacuum scale.
    eps_r : float
        Relative permittivity.
    m_eff_h, m_eff_e : float
        Hole and electron effective masses in units of the electron rest mass.
    """

    name: str
    E_c: float
    E_v: float
    eps_r: float
    m_eff_h: float
    m_eff_e: float

    @property
    def E_g(self) -> float:
        """Band gap, eV."""
        return self.E_c - self.E_v

    def __post_init__(self) -> None:
        if self.E_c <= self.E_v:
            raise ValueError(
                f"{self.name}: E_c ({self.E_c}) must lie above E_v ({self.E_v})"
            )


@dataclass(frozen=True)
class Device:
    """Sample geometry and measurement temperature.

    Attributes
    ----------
    thickness : float
        Inter-electrode spacing ``L``, m.
    area : float
        Active area ``S``, m^2. Used only to convert measured current to current
        density; the model itself depends on thickness alone.
    temperature : float
        Measurement temperature ``T``, K.
    """

    thickness: float
    area: float
    temperature: float = 300.0

    def __post_init__(self) -> None:
        for name in ("thickness", "area", "temperature"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive, got {getattr(self, name)}")


@dataclass(frozen=True)
class ModelParams:
    """The five parameters varied when fitting.

    These are the quantities listed as "Variable" in Supplementary Table S2.

    Attributes
    ----------
    mu_0 : float
        Microscopic mobility, m^2 V^-1 s^-1.
    N_t : float
        Concentration of localized (trap) states, m^-3.
    E_t : float
        Energy of the trap distribution maximum, eV.
    T_t : float
        Characteristic trap temperature, K. Sets the width of the distribution.
    E_F0 : float
        Thermodynamic Fermi level at zero bias, eV.
    """

    mu_0: float
    N_t: float
    E_t: float
    T_t: float
    E_F0: float

    def __post_init__(self) -> None:
        if self.mu_0 <= 0:
            raise ValueError(f"mu_0 must be positive, got {self.mu_0}")
        if self.N_t <= 0:
            raise ValueError(f"N_t must be positive, got {self.N_t}")
        if self.T_t <= 0:
            raise ValueError(f"T_t must be positive, got {self.T_t}")


# --------------------------------------------------------------------------- #
# Energy grid
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class EnergyGrid:
    """Quadrature grids for the occupation integrals.

    The band and trap integrands need different treatment, so they get separate
    grids rather than one shared mesh.

    The square-root band edge has an infinite derivative, and the trapezoid rule
    on it converges only as :math:`O(h^{3/2})`. On a uniform 1 meV mesh the
    resulting error in ``p_f`` is around 2e-3 relative, which is larger than the
    Fermi-Dirac correction to the Boltzmann limit that the integral exists to
    capture. Substituting :math:`E = E_v - t^2` absorbs the singularity into the
    Jacobian, leaving a smooth polynomial integrand:

    .. math::
        \\int_{-\\infty}^{E_v} C\\sqrt{E_v - E}\\,[1-f]\\,dE
        = \\int_0^{t_{max}} 2 C t^2 \\,[1-f]\\,dt

    which converges to machine precision at a couple of thousand points.

    The trap integrand is smooth, so a uniform mesh local to ``E_t`` is fine; it
    only has to be fine relative to ``k_B T_t``, which can be well under 1 meV.

    Attributes
    ----------
    t_band : FloatArray
        Substitution variable for the band integrals, running from 0 to
        ``sqrt(band_span)``.
    E_trap : FloatArray
        Uniform energies covering the localized states, eV.
    """

    t_band: FloatArray
    E_trap: FloatArray

    @classmethod
    def build(
        cls,
        material: Material,
        params: ModelParams,
        constants: Constants = CODATA,
        *,
        band_span: float = 3.0,
        band_points: int = 2001,
        trap_halfwidth: float = 30.0,
        trap_points: int = 1001,
    ) -> EnergyGrid:
        """Construct the band and trap quadrature grids.

        Parameters
        ----------
        material, params
            Define the band edges and the trap position and width.
        constants
            Used for ``k_B T_t``.
        band_span
            How far beyond each band edge to integrate, eV. The hole occupation
            falls off as ``exp(-(E_F - E)/k_B T)``, so 3 eV is many hundreds of
            ``k_B T`` and the truncation is negligible.
        band_points
            Points in the substituted band variable. The default is converged to
            about 1e-5 relative; raising it costs time without buying accuracy,
            and past ~4000 the intermediate matrices stop fitting in cache and
            the whole evaluation slows by nearly an order of magnitude.
        trap_halfwidth
            Half-width of the trap mesh in units of ``k_B T_t``. Both
            biexponential forms fall off as ``e^{-|u|}``, so 30 truncates at
            order 1e-13.
        trap_points
            Points in the trap mesh.

        Notes
        -----
        The two grids are treated as having disjoint support: the band integrals
        see only the band DOS and the trap integral only the localized states.
        The trap's exponential tails do reach the band edges in principle, but
        for the parameters in this system ``|E_t - E_v|`` is some hundreds of
        ``k_B T_t``, making the overlap around ``e^{-300}``.
        """
        if band_span <= 0:
            raise ValueError(f"band_span must be positive, got {band_span}")

        t_band = np.linspace(0.0, np.sqrt(band_span), band_points)
        kT_t = constants.kT_eV(params.T_t)
        E_trap = np.linspace(
            params.E_t - trap_halfwidth * kT_t,
            params.E_t + trap_halfwidth * kT_t,
            trap_points,
        )
        return cls(t_band=t_band, E_trap=E_trap)

    @property
    def band_span(self) -> float:
        """Energy range covered by the band integrals, eV."""
        return float(self.t_band[-1] ** 2)

    def energies(self, material: Material) -> FloatArray:
        """Union of all sampled energies, sorted. For plotting the DOS."""
        return np.unique(
            np.concatenate(
                [
                    material.E_v - self.t_band**2,
                    material.E_c + self.t_band**2,
                    self.E_trap,
                ]
            )
        )


# --------------------------------------------------------------------------- #
# Density of states
# --------------------------------------------------------------------------- #


def _band_prefactor(m_eff: float, constants: Constants) -> float:
    """Square-root band prefactor ``C`` of Eqs (S1), (S2), in m^-3 eV^-3/2."""
    return (
        4.0
        * np.pi
        * (2.0 * constants.m_0 * m_eff * constants.e) ** 1.5
        / constants.h**3
    )


def valence_band_dos(
    E: FloatArray, material: Material, constants: Constants = CODATA
) -> FloatArray:
    """Valence-band density of states, Eq (S2).

    ``g(E) = C_h sqrt(E_v - E)`` for ``E <= E_v``, zero above.
    """
    C = _band_prefactor(material.m_eff_h, constants)
    return np.where(E <= material.E_v, C * np.sqrt(np.maximum(material.E_v - E, 0.0)), 0.0)


def conduction_band_dos(
    E: FloatArray, material: Material, constants: Constants = CODATA
) -> FloatArray:
    """Conduction-band density of states, Eq (S1).

    ``g(E) = C_e sqrt(E - E_c)`` for ``E >= E_c``, zero below.
    """
    C = _band_prefactor(material.m_eff_e, constants)
    return np.where(E >= material.E_c, C * np.sqrt(np.maximum(E - material.E_c, 0.0)), 0.0)


def effective_dos(
    m_eff: float, temperature: float, constants: Constants = CODATA
) -> float:
    """Effective density of states ``N_c`` or ``N_v``, Eqs (S3), (S4).

    ``N = 2 (2 pi m* m_0 k_B T / h^2)^{3/2}``, in m^-3.

    Notes
    -----
    Supplementary Table S4 lists values that do not follow from this expression
    (4.52e16 cm^-3 for MAPbBr3, 5.17e10 cm^-3 for MAPbI3, against roughly
    4.2e18 cm^-3 from the formula). This function follows the formula, which is
    also what the prototype computes and what the standard
    2.5e19 cm^-3 (m*/m_0)^{3/2} room-temperature value reproduces. The MAPbI3
    entry in the table is six orders from anything Eq (S4) can give.
    """
    return 2.0 * (
        2.0
        * np.pi
        * m_eff
        * constants.m_0
        * constants.k_B
        * temperature
        / constants.h**2
    ) ** 1.5


def trap_dos(
    E: FloatArray,
    params: ModelParams,
    constants: Constants = CODATA,
    *,
    profile: TrapProfile = DEFAULT_TRAP_PROFILE,
) -> FloatArray:
    """Density of localized states, Eq (S5) or Eq (S6).

    Parameters
    ----------
    E
        Energies, eV.
    params
        Supplies ``N_t``, ``E_t`` and ``T_t``.
    constants
        Used for ``k_B T_t``.
    profile
        Which functional form to use. See :class:`TrapProfile`; the variants do
        not all carry the same normalization.

    Returns
    -------
    FloatArray
        Density of states, m^-3 eV^-1.

    Notes
    -----
    ``SI_BIEXPONENTIAL`` and ``GAUSSIAN`` integrate to ``N_t``.
    ``PROTOTYPE_SECH`` integrates to ``(pi/4) N_t``; this is reproduced
    deliberately, not accidentally, and :func:`trap_dos_norm` reports it.
    """
    kT_t = constants.kT_eV(params.T_t)
    u = (E - params.E_t) / kT_t

    if profile is TrapProfile.SI_BIEXPONENTIAL:
        # e^u/(1+e^u)^2 == sech^2(u/2)/4, written via expit for overflow safety.
        s = expit(u)
        return (params.N_t / kT_t) * s * (1.0 - s)

    if profile is TrapProfile.PROTOTYPE_SECH:
        # e^u/(1+e^{2u}) == 1/(2 cosh u). Evaluated at -|u| so the exponentials
        # decay rather than overflow; cosh itself overflows for |u| > ~710.
        a = -np.abs(u)
        return (params.N_t / (2.0 * kT_t)) * np.exp(a) / (1.0 + np.exp(2.0 * a))

    if profile is TrapProfile.GAUSSIAN:
        sigma = 2.0 * kT_t
        return (
            params.N_t
            / (sigma * np.sqrt(2.0 * np.pi))
            * np.exp(-0.5 * ((E - params.E_t) / sigma) ** 2)
        )

    raise ValueError(f"unknown trap profile: {profile!r}")


def trap_dos_norm(
    params: ModelParams,
    constants: Constants = CODATA,
    *,
    profile: TrapProfile = DEFAULT_TRAP_PROFILE,
) -> float:
    """Analytic value of ``int g_t(E) dE`` for the given profile, m^-3.

    Useful as a check that a numerical grid resolves the trap, and as an
    explicit statement of which profiles are normalized to ``N_t``.
    """
    if profile is TrapProfile.PROTOTYPE_SECH:
        return float(np.pi / 4.0 * params.N_t)
    return float(params.N_t)


# --------------------------------------------------------------------------- #
# Occupation
# --------------------------------------------------------------------------- #


def fermi_dirac(
    E: FloatArray, E_F: float | FloatArray, temperature: float,
    constants: Constants = CODATA,
) -> FloatArray:
    """Electron occupation ``f(E - E_F)``, Eq (S9). Overflow-safe."""
    return expit(-(E - E_F) / constants.kT_eV(temperature))


@dataclass(frozen=True)
class CarrierDensities:
    """Carrier concentrations as a function of Fermi level.

    Attributes
    ----------
    E_F : FloatArray
        Fermi levels the quantities are evaluated at, eV.
    p_f, p_t, p_s : FloatArray
        Free, trapped and total hole concentrations, m^-3. Absolute, not
        referenced to equilibrium.
    p_s_injected : FloatArray
        ``p_s(E_F) - p_s(E_F0)``: the charge actually injected, which is what
        forms the space charge. Zero at ``E_F = E_F0``.
    p_f_injected : FloatArray
        ``p_f(E_F) - p_f(E_F0)``: the free part of the injected charge alone.
        The numerator of ``ThetaModel.INJECTED_OVER_INJECTED``.
    n_f, n_t, n_s, n_s_injected, n_f_injected : FloatArray
        The electron counterparts.
    theta_p, theta_n : FloatArray
        ``Theta``, Eq (S12), under whichever :class:`ThetaModel` was requested.
        Bounded by 1 for every reading except the default; see ``valid``.
    theta_model : ThetaModel
        Which reading ``theta_p`` and ``theta_n`` carry.
    valid : FloatArray
        Boolean mask of points whose outputs are physically meaningful. **Read
        ``theta`` and anything derived from it through this.** False where
        ``Theta > 1``, since ``mu_eff = mu_0 Theta`` cannot exceed the
        microscopic mobility. That happens at the low-bias end, where the
        equilibrium free population still rivals the injected charge and the
        space-charge picture does not yet apply.
    """

    E_F: FloatArray
    p_f: FloatArray
    p_t: FloatArray
    p_s: FloatArray
    p_s_injected: FloatArray
    p_f_injected: FloatArray
    n_f: FloatArray
    n_t: FloatArray
    n_s: FloatArray
    n_s_injected: FloatArray
    n_f_injected: FloatArray
    theta_p: FloatArray
    theta_n: FloatArray
    valid: FloatArray
    theta_model: ThetaModel = DEFAULT_THETA_MODEL


#: Rough cap on the number of float64 entries in one intermediate matrix, used
#: to size the chunks in :func:`_absolute_densities`. 4e6 entries is about 32 MB,
#: which stays comfortably in cache-friendly territory; the previous unchunked
#: form allocated ~100 MB per temporary and ran an order of magnitude slower.
_MATRIX_BUDGET: Final[int] = 4_000_000


def _trapezoid_weights(x: FloatArray) -> FloatArray:
    """Composite-trapezoid quadrature weights for the sample points ``x``."""
    w = np.empty_like(x)
    w[1:-1] = 0.5 * (x[2:] - x[:-2])
    w[0] = 0.5 * (x[1] - x[0])
    w[-1] = 0.5 * (x[-1] - x[-2])
    return w


def _absolute_densities(
    E_F: FloatArray,
    material: Material,
    device: Device,
    params: ModelParams,
    constants: Constants,
    grid: EnergyGrid,
    profile: TrapProfile,
) -> tuple[
    FloatArray, FloatArray, FloatArray, FloatArray,
    FloatArray, FloatArray, FloatArray,
]:
    """Occupation integrals, m^-3.

    Returns ``(p_f, p_t, n_f, n_t, p_s_injected, p_f_injected, n_f_injected)``.
    The first four are absolute; the last three are referenced to ``E_F0``.

    Band integrals use the square-root substitution described in
    :class:`EnergyGrid`; the trap integral uses a uniform local mesh. Each is a
    single matrix contraction, O(len(E_F) x len(grid)) in time and memory.
    """
    t = grid.t_band
    E_v_grid = material.E_v - t**2
    E_c_grid = material.E_c + t**2

    # dE = 2t dt and sqrt(E_v - E) = t, so g_v(E) dE -> 2 C t^2 dt.
    jac_v = 2.0 * _band_prefactor(material.m_eff_h, constants) * t**2
    jac_c = 2.0 * _band_prefactor(material.m_eff_e, constants) * t**2

    E_t_grid = grid.E_trap
    g_t = trap_dos(E_t_grid, params, constants, profile=profile)

    # Fold the trapezoid weights into the integrand so each integral becomes a
    # matrix-vector product. This is a BLAS call rather than a chain of large
    # temporaries, which matters: at 3000 Fermi levels the naive form allocates
    # ~100 MB per intermediate and spends most of its time in the memory system.
    w_band = _trapezoid_weights(t)
    w_trap = _trapezoid_weights(E_t_grid)
    wv, wc = w_band * jac_v, w_band * jac_c
    wt = w_trap * g_t

    n_EF = E_F.size
    p_f = np.empty(n_EF)
    p_t = np.empty(n_EF)
    n_f = np.empty(n_EF)
    n_t = np.empty(n_EF)
    # Injected charge, accumulated with the occupancy difference taken INSIDE
    # the integral. The equivalent "integrate both, then subtract" form differences
    # two enormous far-band occupancies and loses every significant digit: the
    # conduction-band term alone is ~1e28 with an ulp near 1e12, comparable to
    # the answer. See the note in carrier_densities.
    injected = np.empty(n_EF)
    injected_pf = np.empty(n_EF)
    injected_nf = np.empty(n_EF)
    f0_v = fermi_dirac(E_v_grid, params.E_F0, device.temperature, constants)
    f0_c = fermi_dirac(E_c_grid, params.E_F0, device.temperature, constants)
    f0_t = fermi_dirac(E_t_grid, params.E_F0, device.temperature, constants)

    # Chunk so peak memory stays bounded no matter how long the sweep is.
    chunk = max(1, _MATRIX_BUDGET // max(t.size, E_t_grid.size))
    for lo in range(0, n_EF, chunk):
        sl = slice(lo, min(lo + chunk, n_EF))
        ef = E_F[sl][:, None]
        T = device.temperature

        f_v = fermi_dirac(E_v_grid[None, :], ef, T, constants)
        p_f[sl] = (1.0 - f_v) @ wv
        d_pf = (f0_v[None, :] - f_v) @ wv
        injected_pf[sl] = d_pf
        acc = d_pf
        del f_v

        f_c = fermi_dirac(E_c_grid[None, :], ef, T, constants)
        n_f[sl] = f_c @ wc
        injected_nf[sl] = (f_c - f0_c[None, :]) @ wc
        acc = acc + (f0_c[None, :] - f_c) @ wc
        del f_c

        f_t = fermi_dirac(E_t_grid[None, :], ef, T, constants)
        p_t[sl] = (1.0 - f_t) @ wt
        n_t[sl] = f_t @ wt
        acc = acc + (f0_t[None, :] - f_t) @ wt
        del f_t

        injected[sl] = acc

    return p_f, p_t, n_f, n_t, injected, injected_pf, injected_nf


def carrier_densities(
    E_F: FloatArray,
    material: Material,
    device: Device,
    params: ModelParams,
    constants: Constants = CODATA,
    *,
    grid: EnergyGrid | None = None,
    profile: TrapProfile = DEFAULT_TRAP_PROFILE,
    theta_model: ThetaModel = DEFAULT_THETA_MODEL,
) -> CarrierDensities:
    """Evaluate the occupation integrals of Eqs (S7), (S8), (S10), (S11).

    For holes,

    .. math::
        p_f(E_F) = \\int g_\\mathrm{VB}(E)\\,[1 - f(E - E_F)]\\,dE

    over the valence band, and the trapped concentration is the same integral
    over the localized states. The absolute total is ``p_s = p_f + p_t``, and
    the injected total, which is what forms the space charge, is
    ``p_s(E_F) - p_s(E_F0)``.

    Parameters
    ----------
    E_F
        Fermi levels to evaluate at, eV.
    material, device, params
        Model inputs. ``params.E_F0`` fixes the equilibrium reference.
    grid
        Energy grid. Built by :meth:`EnergyGrid.build` if omitted.
    profile
        Trap functional form, see :class:`TrapProfile`.
    theta_model
        Which reading of Eq (S12) to return as ``theta_p``/``theta_n``. See
        :class:`ThetaModel`; the default is the prototype's.

    Returns
    -------
    CarrierDensities

    Notes
    -----
    Eqs (S8) and (S11) as printed give the total concentration as one integral
    with ``E_F0`` as a limit. The prototype spreadsheet implements this as an
    explicit subtraction inside the column,
    ``SUMPRODUCT(...) - ps0``, so ``E_F0`` acts as an equilibrium reference
    rather than a range of integration. That is what is done here, and it
    reproduces the prototype's ``n(E)!J`` and ``n(E)!L`` to 0.1-2 %.

    The sources support three readings of Eq (S12) which agree at high injection
    and differ by up to a factor of six through the trap-filling region, so the
    choice is not a detail. It is exposed as ``theta_model``; the default divides
    the absolute ``p_f`` by the **injected** total, matching both the prototype
    and the analysis branch, where Eq (5) likewise yields an absolute ``p_f``.
    Only that reading is unbounded, and where it exceeds 1 it clears ``valid``.
    See :class:`ThetaModel` and ``ASCLC_spec.md`` section 7.3.
    """
    E_F = np.atleast_1d(np.asarray(E_F, dtype=float))
    if grid is None:
        grid = EnergyGrid.build(material, params, constants)

    # Evaluate the equilibrium reference in the SAME array as the requested
    # Fermi levels. Computing it in a separate call lets the chunked matrix
    # products round it differently, which leaves p_s_injected at a small
    # nonzero value (sometimes negative) where it must be exactly zero.
    p_f, p_t, n_f, n_t, p_s_injected, p_f_injected, n_f_injected = (
        _absolute_densities(E_F, material, device, params, constants, grid, profile)
    )
    p_s = p_f + p_t
    n_s = n_f + n_t

    # Every injected hole is an electron removed from the same set of states, so
    # the electron injected total is the exact negative. Computing it instead as
    # (n_f + n_t) minus its equilibrium value omits the valence-band term, which
    # for these parameters is 136 times larger than what remains.
    n_s_injected = -p_s_injected

    # Theta is a magnitude ratio throughout, as the prototype's ABS(nf/ns) makes
    # explicit: on a hole sweep the electron populations are depleted, so the
    # electron numerator and denominator are both negative.
    if theta_model is ThetaModel.ABSOLUTE_OVER_INJECTED:
        # The prototype's reading. The mixed reference means Theta can exceed 1
        # where the equilibrium free population is comparable to the injected
        # charge; that is a domain violation, reported through `valid` rather
        # than hidden.
        num_p, den_p, num_n, den_n = p_f, p_s_injected, n_f, n_s_injected
    elif theta_model is ThetaModel.ABSOLUTE_OVER_TOTAL:
        num_p, den_p, num_n, den_n = p_f, p_s, n_f, n_s
    elif theta_model is ThetaModel.INJECTED_OVER_INJECTED:
        num_p, den_p = p_f_injected, p_s_injected
        num_n, den_n = n_f_injected, n_s_injected
    else:
        raise ValueError(f"unknown theta model: {theta_model!r}")

    with np.errstate(divide="ignore", invalid="ignore"):
        theta_p = np.where(den_p != 0, np.abs(num_p / den_p), np.nan)
        theta_n = np.where(den_n != 0, np.abs(num_n / den_n), np.nan)

    with np.errstate(invalid="ignore"):
        valid = ~(theta_p > 1.0)

    return CarrierDensities(
        E_F=E_F, p_f=p_f, p_t=p_t, p_s=p_s, p_s_injected=p_s_injected,
        p_f_injected=p_f_injected,
        n_f=n_f, n_t=n_t, n_s=n_s, n_s_injected=n_s_injected,
        n_f_injected=n_f_injected,
        theta_p=theta_p, theta_n=theta_n, valid=valid, theta_model=theta_model,
    )


# --------------------------------------------------------------------------- #
# Model branch
# --------------------------------------------------------------------------- #


def gamma_from_trap_temperature(
    params: ModelParams,
    device: Device,
    *,
    model: GammaModel = DEFAULT_GAMMA_MODEL,
) -> float:
    """Reverse slope ``gamma`` for the model branch, step M4.

    See :class:`GammaModel` for the three candidate expressions and for why the
    prototype's is the default.
    """
    if model is GammaModel.TT_OVER_T:
        gamma = params.T_t / device.temperature
    elif model is GammaModel.SI_NOTE_M4:
        # Note M4 as printed: T/(T_t + T) for a trap hotter than the sample,
        # otherwise the constant 0.5. Every configuration in the article and the
        # prototype has T_t < T, so they all take the second branch.
        if params.T_t >= device.temperature:
            gamma = device.temperature / (params.T_t + device.temperature)
        else:
            gamma = 0.5
    elif model is GammaModel.TT_OVER_T_PLUS_TT:
        gamma = params.T_t / (device.temperature + params.T_t)
    else:
        raise ValueError(f"unknown gamma model: {model!r}")

    if not 0.0 <= gamma <= 0.5:
        warnings.warn(
            f"gamma = {gamma:.4g} lies outside the [0, 0.5] range that "
            "Supplementary Note M4 states for charge injection",
            RuntimeWarning,
            stacklevel=2,
        )
    return float(gamma)


@dataclass(frozen=True)
class ModelCurve:
    """Output of the model branch, parametrized by Fermi level.

    Attributes
    ----------
    E_F : FloatArray
        Fermi levels swept, eV.
    V : FloatArray
        Voltage from Eq (6) inverted, V.
    J : FloatArray
        Current density from Eq (S14), A m^-2.
    p_f : FloatArray
        Free hole concentration, m^-3. Matches the prototype's ``n(E)!L``.
    p_t : FloatArray
        Absolute trapped hole concentration, m^-3, including the population
        already present at zero bias. **This is not the article's** ``ptm``.
        See ``p_space_charge``.
    p_space_charge : FloatArray
        The charge that drives the voltage through the inverted Eq (6), m^-3.
        By default the total injected charge ``p_s(E_F) - p_s(E_F0)``.

        **This is the article's** ``ptm``, the curve plotted against energy in
        Fig. 4c,d and Fig. 6: it equals the prototype's ``MODEL!AJ``, which reads
        from ``n(E)!O``. Plotting ``p_t`` there instead gives a visibly different
        curve, since ``p_t`` carries the equilibrium trapped population and
        misses the band contribution at high injection.
    theta : FloatArray
        ``Theta``; see :func:`carrier_densities` and :class:`ThetaModel`.
    valid : FloatArray
        Boolean mask of physically meaningful points, false where ``Theta > 1``.
        Read ``theta`` and ``mu_eff`` through it. Carries the same meaning as
        :attr:`AnalysisResult.valid`. Only the default ``theta_model`` can put
        points outside the domain; the other two readings are bounded by 1.
    mu_eff : FloatArray
        Effective mobility ``mu_0 Theta``, m^2 V^-1 s^-1.
    gamma : float
        The constant reverse slope used.
    theta_model : ThetaModel
        Which reading of Eq (S12) ``theta`` and ``mu_eff`` carry.
    """

    E_F: FloatArray
    V: FloatArray
    J: FloatArray
    p_f: FloatArray
    p_t: FloatArray
    p_space_charge: FloatArray
    theta: FloatArray
    valid: FloatArray
    mu_eff: FloatArray
    gamma: float
    theta_model: ThetaModel = DEFAULT_THETA_MODEL


def model_curve(
    params: ModelParams,
    material: Material,
    device: Device,
    constants: Constants = CODATA,
    *,
    n_points: int = 3001,
    E_F_span: float | None = None,
    V_max: float | None = None,
    profile: TrapProfile = DEFAULT_TRAP_PROFILE,
    gamma_model: GammaModel = DEFAULT_GAMMA_MODEL,
    space_charge: SpaceCharge = DEFAULT_SPACE_CHARGE,
    theta_model: ThetaModel = DEFAULT_THETA_MODEL,
    grid: EnergyGrid | None = None,
) -> ModelCurve:
    """Generate a modelled J-V curve, steps M1-M5.

    The Fermi level is swept from ``E_F0`` toward the valence band, which is the
    direction hole injection drives it. For each ``E_F`` the trapped hole
    concentration gives the voltage by inverting Eq (6),

    .. math:: V = e L^2 p_t / [\\varepsilon_0 \\varepsilon_r (1-\\gamma)(2-\\gamma)]

    and the free concentration then gives the current from Eq (S14),

    .. math:: J = e \\mu_0 (2 - \\gamma) p_f V / L

    Parameters
    ----------
    params, material, device
        Model inputs.
    n_points
        Number of Fermi levels in the sweep.
    E_F_span
        How far below ``E_F0`` to sweep, eV. Defaults to reaching ``E_v``.
    V_max
        If given, points above this voltage are dropped. Sweeping E_F all the
        way to the band edge drives the injected charge, and hence the voltage,
        to physically meaningless values, so a bound near the top of the
        measured range is normally what you want.
    profile, gamma_model, space_charge, theta_model
        See :class:`TrapProfile`, :class:`GammaModel`, :class:`SpaceCharge` and
        :class:`ThetaModel`. Every default reproduces the prototype spreadsheet.
    grid
        Energy grid for the occupation integrals.

    Returns
    -------
    ModelCurve

    Notes
    -----
    This is a pure function of its arguments. To fit, wrap a residual around it
    and hand that to ``scipy.optimize.least_squares``; nothing here needs to
    change.
    """
    if E_F_span is None:
        E_F_span = params.E_F0 - material.E_v
    if E_F_span <= 0:
        raise ValueError(
            f"E_F0 ({params.E_F0}) must lie above E_v ({material.E_v}) for hole "
            "injection; pass E_F_span explicitly to override"
        )

    E_F = np.linspace(params.E_F0, params.E_F0 - E_F_span, n_points)
    dens = carrier_densities(
        E_F, material, device, params, constants,
        grid=grid, profile=profile, theta_model=theta_model,
    )
    gamma = gamma_from_trap_temperature(params, device, model=gamma_model)

    if space_charge is SpaceCharge.INJECTED_TOTAL:
        p_t_sc = dens.p_s_injected
    elif space_charge is SpaceCharge.INJECTED_TRAPPED:
        p_t_sc = dens.p_t - dens.p_t[0]
    elif space_charge is SpaceCharge.TOTAL:
        p_t_sc = dens.p_t
    else:
        raise ValueError(f"unknown space-charge reference: {space_charge!r}")

    factor = (1.0 - gamma) * (2.0 - gamma)
    V = (
        constants.e
        * device.thickness**2
        * p_t_sc
        / (constants.eps_0 * material.eps_r * factor)
    )
    J = constants.e * params.mu_0 * (2.0 - gamma) * dens.p_f * V / device.thickness

    if V_max is not None:
        keep = V <= V_max
        if not keep.any():
            raise ValueError(
                f"V_max = {V_max} is below the smallest modelled voltage "
                f"({np.nanmin(V):.4g} V); widen it or move E_F0"
            )
        E_F, V, J, p_t_sc = E_F[keep], V[keep], J[keep], p_t_sc[keep]
        # Filter every per-Fermi-level array, found by type rather than named
        # one at a time: an explicit list silently leaves any field added later
        # at full length, which then misaligns against E_F.
        dens = replace(
            dens,
            **{
                f.name: getattr(dens, f.name)[keep]
                for f in fields(dens)
                if isinstance(getattr(dens, f.name), np.ndarray)
            },
        )

    return ModelCurve(
        E_F=E_F,
        V=V,
        J=J,
        p_f=dens.p_f,
        p_t=dens.p_t,
        p_space_charge=p_t_sc,
        theta=dens.theta_p,
        valid=dens.valid,
        mu_eff=params.mu_0 * dens.theta_p,
        gamma=gamma,
        theta_model=theta_model,
    )


# --------------------------------------------------------------------------- #
# Analysis branch
# --------------------------------------------------------------------------- #


def local_loglog_slope(
    V: FloatArray, J: FloatArray, window: int = DEFAULT_SLOPE_WINDOW
) -> FloatArray:
    """Local logarithmic slope ``m = d ln J / d ln V``, step A2.

    A centred ordinary least-squares fit of ``ln|J|`` against ``ln|V|`` over a
    sliding window of ``window`` points. Windows are truncated near the ends
    rather than padded, so the first and last few slopes rest on fewer points.

    Parameters
    ----------
    V, J
        Voltage and current density. Must be the same length and strictly
        positive in magnitude; non-finite logarithms yield ``nan`` slopes.
    window
        Number of points per fit. Must be odd and at least 3.

    Returns
    -------
    FloatArray
        Slope at each point. ``nan`` where the window holds fewer than two
        usable points or the voltages within it do not vary.

    Notes
    -----
    This regresses ``ln|J|`` on ``ln|V|`` and :func:`analyse_jv` then takes
    ``gamma = 1/m``. The prototype spreadsheet instead regresses the other way
    round, ``LINEST(ln V, ln j)`` in ``Data-calculations!K``, obtaining gamma
    directly. The article defines the two as equal (``gamma = 1/m =
    d ln V/d ln J``) and they are for noise-free data, but ordinary least
    squares is not symmetric: with scatter, ``1/slope(y|x) != slope(x|y)``, and
    the gap widens as the fit degrades. Fitting in the ``J`` direction is kept
    here because it stays conditioned through the trap-filled-limit region,
    where ``ln V`` is nearly constant and the prototype's direction regresses
    against a near-degenerate abscissa.

    Raises
    ------
    ValueError
        If ``window`` is even, smaller than 3, or larger than the data.
    """
    V = np.asarray(V, dtype=float)
    J = np.asarray(J, dtype=float)
    if V.shape != J.shape:
        raise ValueError(f"V and J must have the same shape, got {V.shape} and {J.shape}")
    if window < 3:
        raise ValueError(
            f"window must be at least 3, got {window}; a 1-point window makes the "
            "regression degenerate and silently returns a slope of 0"
        )
    if window % 2 == 0:
        raise ValueError(f"window must be odd so it can be centred, got {window}")
    if window > V.size:
        raise ValueError(f"window ({window}) exceeds the number of points ({V.size})")

    with np.errstate(divide="ignore", invalid="ignore"):
        x = np.log(np.abs(V))
        y = np.log(np.abs(J))

    n = V.size
    half = window // 2
    slope = np.full(n, np.nan)

    for i in range(n):
        lo, hi = max(0, i - half), min(n, i + half + 1)
        xs, ys = x[lo:hi], y[lo:hi]
        ok = np.isfinite(xs) & np.isfinite(ys)
        if ok.sum() < 2:
            continue
        xs, ys = xs[ok], ys[ok]
        if np.ptp(xs) == 0:
            continue
        slope[i] = np.polyfit(xs, ys, 1)[0]

    return slope


@dataclass(frozen=True)
class AnalysisResult:
    """Output of the analysis branch, parametrized by voltage.

    Attributes
    ----------
    V, J : FloatArray
        The input curve, V and A m^-2.
    m : FloatArray
        Local logarithmic slope.
    gamma : FloatArray
        ``1 / m``.
    mu_eff : FloatArray
        Effective mobility from Eq (4), m^2 V^-1 s^-1.
    p_f : FloatArray
        Free hole concentration, Eq (5), m^-3.
    p_t : FloatArray
        The space charge from Eq (6), m^-3. The paper labels this the trapped
        concentration; it is really the total space charge, which is why
        ``theta`` divides by it rather than by ``p_f + p_t``.
    theta : FloatArray
        ``mu_eff / mu_0``, Eq (4), identically ``p_f / p_t``.
    E_F : FloatArray
        Fermi level from Eq (7), eV.
    valid : FloatArray
        Boolean mask of the points whose outputs are physically meaningful.
        **Read every other field only where this is true.** A point is invalid
        when the current changes sign inside its fit window, so the log-log slope
        is not a slope of anything; when ``gamma`` falls outside (0, 1), where
        Eqs (4) to (6) are out of domain and can return negative concentrations;
        or when ``Theta`` exceeds 1, since ``mu_eff = mu_0 Theta`` cannot exceed
        the microscopic mobility. Carries the same meaning as
        :attr:`ModelCurve.valid`.

        Both happen routinely in the ohmic region, where the measured current is
        noise straddling zero. In the reference dataset 113 of 214 points fail
        one or the other. The values at those points are still returned rather
        than blanked, since a plot of what was rejected is often the fastest way
        to see whether a measurement is usable at all.
    """

    V: FloatArray
    J: FloatArray
    m: FloatArray
    gamma: FloatArray
    mu_eff: FloatArray
    p_f: FloatArray
    p_t: FloatArray
    theta: FloatArray
    E_F: FloatArray
    valid: FloatArray

    @property
    def n_valid(self) -> int:
        """Number of points where the extraction is meaningful."""
        return int(np.count_nonzero(self.valid))


def analyse_jv(
    V: FloatArray,
    J: FloatArray,
    material: Material,
    device: Device,
    mu_0: float,
    constants: Constants = CODATA,
    *,
    window: int = DEFAULT_SLOPE_WINDOW,
    N_v: float | FloatArray | None = None,
    temperature: float | FloatArray | None = None,
) -> AnalysisResult:
    """Extract SCLC quantities from a measured J-V curve, steps A2-A5.

    Parameters
    ----------
    V, J
        Measured voltage and current density.
    material, device
        Sample properties.
    mu_0
        Microscopic mobility, needed by Eq (5). Estimated from the ohmic region
        or from the Mott-Gurney limit, per step A3.
    window
        Points per local slope fit, see :func:`local_loglog_slope`.
    N_v
        Effective valence-band density of states. Computed from Eq (S4) at
        ``temperature`` if omitted.
    temperature
        Sample temperature for Eq (7), K. A scalar, or one value per data point.
        Defaults to ``device.temperature``.

        Pass the measured per-point temperature when you have it. The reference
        prototype does, and it is not a refinement: in its dataset the recorded
        temperature ranges over 282-314 K, which moves the extracted ``E_F`` by
        0.075 eV. The Fermi level shifts the article reports are 0.046 eV and
        0.006 eV, so the temperature scatter is larger than the effect being
        measured. Holding T fixed at its nominal value silently attributes that
        scatter to the physics.

    Returns
    -------
    AnalysisResult

    Notes
    -----
    Results carry a ``valid`` mask; read every other field through it. Points in
    the ohmic region, where the measured current is noise straddling zero, produce
    slopes that are not slopes and ``gamma`` outside its domain.

    ``theta`` is ``mu_eff/mu_0`` per Eq (4), which is identically ``p_f/p_t``
    with those taken from Eqs (5) and (6). It is *not* ``p_f/(p_f + p_t)``:
    Eq (6) returns the total space charge despite being labelled ``p_t``, the
    same conflation that appears in the model branch. The two forms agree while
    traps dominate and diverge as ``Theta`` approaches 1, where
    ``p_f/(p_f + p_t)`` saturates at 1/2.

    Every returned quantity except ``E_F`` depends on ``gamma``, so the choice
    of ``window`` propagates into all of them. In particular Eq (6) makes
    ``p_t`` proportional to ``(1-gamma)(2-gamma)``, which is maximal at
    ``gamma = 0``: a degenerate slope fit does not merely add noise, it returns
    the upper envelope of ``p_t`` for every point. This is why
    :func:`local_loglog_slope` refuses a window below 3.
    """
    V = np.asarray(V, dtype=float)
    J = np.asarray(J, dtype=float)

    T = device.temperature if temperature is None else temperature
    T = np.asarray(T, dtype=float)
    if T.ndim and T.shape != V.shape:
        raise ValueError(
            f"temperature must be scalar or match V; got {T.shape} and {V.shape}"
        )
    if np.any(T <= 0):
        raise ValueError("temperature must be positive")

    # N_v and k_B T are taken at the same temperature. The prototype mixes them,
    # evaluating k_B T per row from the measured value while leaving N_v at the
    # nominal temperature; that is worth about 0.002 eV in E_F over its own
    # 282-314 K range, small next to the 0.075 eV the k_B T term contributes,
    # but there is no reason to reproduce the inconsistency.
    if N_v is None:
        N_v = effective_dos(material.m_eff_h, T, constants)
    N_v = np.asarray(N_v, dtype=float)

    m = local_loglog_slope(V, J, window)
    with np.errstate(divide="ignore", invalid="ignore"):
        gamma = 1.0 / m

        mu_eff = (
            device.thickness**3
            * J
            / (
                constants.eps_0
                * material.eps_r
                * (1.0 - gamma)
                * (2.0 - gamma) ** 2
                * V**2
            )
        )
        p_f = device.thickness * J / (constants.e * mu_0 * (2.0 - gamma) * V)
        p_t = (
            constants.eps_0
            * material.eps_r
            * (1.0 - gamma)
            * (2.0 - gamma)
            * V
            / (constants.e * device.thickness**2)
        )
        # Theta = mu_eff/mu_0 by Eq (4). This is identically p_f/p_t with p_f
        # and p_t from Eqs (5) and (6) -- verified in the tests -- and it is
        # NOT p_f/(p_f + p_t). The quantity Eq (6) returns is the total space
        # charge, which the paper labels p_t; in the trap-dominated regime the
        # two agree, but they diverge as Theta approaches 1, where
        # p_f/(p_f + p_t) saturates at 1/2 instead.
        theta = mu_eff / mu_0
        # Eq (7) is p_f = N_v exp(-dE_F/kT) with dE_F = E_F - E_v >= 0, so
        # inverting gives E_F = E_v + kT ln(N_v/p_f). Writing ln(p_f/N_v) here
        # instead flips the sign and puts E_F below the valence band.
        E_F = material.E_v + constants.kT_eV(T) * np.log(N_v / np.abs(p_f))

    # A window spanning a sign change in J has no log-log slope: log|J| traces a
    # V shape and the fit through it describes nothing. Outside 0 < gamma < 1 the
    # SCLC relations are out of domain.
    half = window // 2
    sign_change = np.zeros(V.size, dtype=bool)
    signs = np.sign(J)
    for i in range(V.size):
        w = signs[max(0, i - half) : min(V.size, i + half + 1)]
        nz = w[w != 0]
        sign_change[i] = nz.size > 0 and not np.all(nz == nz[0])

    with np.errstate(invalid="ignore"):
        gamma_in_domain = (gamma > 0.0) & (gamma < 1.0)
    # Theta > 1 means mu_eff exceeds mu_0, which cannot happen. Flagged here on
    # the same footing as in the model branch, so `valid` means the same thing in
    # both: this point's outputs are physically meaningful.
    with np.errstate(invalid="ignore"):
        theta_ok = ~(theta > 1.0)
    # Split so the diagnostics below report on points that fail only the check
    # they describe, rather than on points already excluded for other reasons.
    fittable = (
        gamma_in_domain
        & ~sign_change
        & np.isfinite(mu_eff)
        & np.isfinite(p_f)
        & np.isfinite(E_F)
    )
    valid = fittable & theta_ok

    if not valid.any():
        warnings.warn(
            "no points survived the validity check: every fit window either "
            "spans a sign change in J or gives gamma outside (0, 1). Check the "
            "sign convention, the window size, and whether the sweep reaches "
            "beyond the ohmic region.",
            RuntimeWarning,
            stacklevel=2,
        )
    elif valid.mean() < 0.5:
        warnings.warn(
            f"only {valid.sum()} of {V.size} points are usable; the rest fall in "
            "the ohmic noise, outside 0 < gamma < 1, or give Theta > 1. Read "
            "results through the `valid` mask.",
            RuntimeWarning,
            stacklevel=2,
        )

    finite_t = np.isfinite(theta) & fittable
    if np.any(theta[finite_t] > 1.0):
        worst = float(np.nanmax(theta[finite_t]))
        warnings.warn(
            f"Theta reaches {worst:.3g} > 1, which is unphysical: mu_eff cannot "
            f"exceed the microscopic mobility. The supplied mu_0 = {mu_0:.4g} is "
            f"too small by at least a factor of {worst:.3g}. Step A3 estimates "
            "mu_0 from the ohmic region or from the Mott-Gurney limit. Those "
            "points are excluded from `valid`.",
            RuntimeWarning,
            stacklevel=2,
        )

    finite = np.isfinite(p_f) & fittable
    N_v_at = N_v[finite] if N_v.ndim else N_v
    if np.any(np.abs(p_f[finite]) > N_v_at):
        warnings.warn(
            "p_f exceeds N_v somewhere, so the semiconductor is degenerate there "
            "and Eq (7) does not apply; the reported E_F is unreliable in that "
            "range.",
            RuntimeWarning,
            stacklevel=2,
        )

    return AnalysisResult(
        V=V, J=J, m=m, gamma=gamma, mu_eff=mu_eff,
        p_f=p_f, p_t=p_t, theta=theta, E_F=E_F, valid=valid,
    )



# --------------------------------------------------------------------------- #
# Measurement files
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Measurement:
    """A loaded J-V measurement, ready for :func:`analyse_jv`.

    Attributes
    ----------
    V : FloatArray
        Applied voltage, V.
    J : FloatArray
        Current density, A m^-2.
    temperature : FloatArray | None
        Sample temperature per point, K, if the file carried one. Pass it to
        :func:`analyse_jv`; see the note there on why it matters.
    source : str
        Path the data came from.
    n_raw : int
        Rows read before binning.
    """

    V: FloatArray
    J: FloatArray
    temperature: FloatArray | None
    source: str
    n_raw: int

    def __len__(self) -> int:
        return int(self.V.size)


def _bin_average(x: FloatArray, size: int) -> FloatArray:
    """Average consecutive groups of ``size`` points, dropping any remainder."""
    if size == 1:
        return x
    n = (x.size // size) * size
    return x[:n].reshape(-1, size).mean(axis=1)


def load_jv(
    path: str,
    device: Device,
    *,
    columns: tuple[int, ...] | None = None,
    area: float | None = None,
    current_is_density: bool = False,
    temperature_in_celsius: bool = False,
    voltage_offset: float = 0.0,
    bin_size: int = 1,
    delimiter: str | None = None,
    skip_header: int | None = None,
) -> Measurement:
    """Read a J-V measurement from a delimited text file.

    Expects columns of voltage, current (or current density) and optionally
    temperature, in that order. Blank lines and comment lines beginning with
    ``#`` are ignored, and a single non-numeric header row is skipped
    automatically.

    Parameters
    ----------
    path
        File to read. Delimiter is sniffed unless given.
    device
        Supplies the electrode area for the current-to-density conversion and
        the fallback temperature.
    columns
        Zero-based indices of (voltage, current) or (voltage, current,
        temperature). Defaults to the first two or three columns present.
    area
        Electrode area, m^2, overriding ``device.area``.
    current_is_density
        Set when the current column is already A m^-2 and needs no division by
        area.
    temperature_in_celsius
        Convert the temperature column from degrees Celsius to kelvin.
    voltage_offset
        Added to the voltage column, matching the prototype's ``Vmin`` shift.
    bin_size
        Average this many consecutive points together. The prototype applies the
        same averaging through its binning cell; it reduces noise in the ohmic
        region, where the raw current can change sign.
    delimiter, skip_header
        Override the sniffed delimiter and header handling.

    Returns
    -------
    Measurement

    Raises
    ------
    ValueError
        If the file has too few columns, no usable rows, or a non-positive area.

    Notes
    -----
    Points are returned in file order and are not sorted or filtered. Rows where
    the current is zero or negative are kept: the analysis branch takes
    logarithms of magnitudes, and dropping them silently would bias the ohmic
    region, where noise legitimately straddles zero.
    """
    if bin_size < 1:
        raise ValueError(f"bin_size must be at least 1, got {bin_size}")

    with open(path, "r", encoding="utf-8-sig") as fh:
        lines = [
            ln for ln in (raw.strip() for raw in fh)
            if ln and not ln.startswith("#")
        ]
    if not lines:
        raise ValueError(f"{path}: no data rows")

    if delimiter is None:
        first = lines[0]
        delimiter = max(("\t", ";", ",", None), key=lambda d: len(first.split(d)))

    def parse(line: str) -> list[float] | None:
        parts = line.split(delimiter)
        try:
            return [float(p.replace(",", ".") if delimiter != "," else p) for p in parts]
        except ValueError:
            return None

    if skip_header is None:
        skip_header = 0 if parse(lines[0]) is not None else 1
    rows = [parse(ln) for ln in lines[skip_header:]]
    data = [r for r in rows if r is not None]
    if not data:
        raise ValueError(f"{path}: no numeric rows after the header")

    width = min(len(r) for r in data)
    if width < 2:
        raise ValueError(f"{path}: need at least voltage and current columns")

    if columns is None:
        columns = (0, 1, 2) if width >= 3 else (0, 1)
    if len(columns) not in (2, 3):
        raise ValueError(f"columns must name 2 or 3 fields, got {len(columns)}")
    if max(columns) >= width:
        raise ValueError(
            f"{path}: column index {max(columns)} beyond the {width} columns present"
        )

    arr = np.array([[r[i] for i in columns] for r in data], dtype=float)
    n_raw = arr.shape[0]

    V = _bin_average(arr[:, 0], bin_size) + voltage_offset
    current = _bin_average(arr[:, 1], bin_size)

    if current_is_density:
        J = current
    else:
        a = device.area if area is None else area
        if a <= 0:
            raise ValueError(f"area must be positive, got {a}")
        J = current / a

    temperature = None
    if len(columns) == 3:
        temperature = _bin_average(arr[:, 2], bin_size)
        if temperature_in_celsius:
            temperature = temperature + 273.15
        if np.any(temperature <= 0):
            raise ValueError(
                f"{path}: non-positive temperature; is the column in Celsius? "
                "pass temperature_in_celsius=True"
            )

    return Measurement(
        V=V, J=J, temperature=temperature, source=str(path), n_raw=n_raw
    )


# --------------------------------------------------------------------------- #
# Parameter files
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Configuration:
    """Everything one run needs, as read from a TOML parameter file.

    Groups the sample description, the five hand-selected parameters and the
    modelling and numerical choices into a single value, so a configuration can
    be version-controlled next to the measurement it belongs to instead of
    living as edits to a module constant.

    Attributes
    ----------
    name : str
        Label for reports and plot titles.
    material, device, params
        The model inputs.
    constants : Constants
        CODATA by default; the prototype's rounded values on request.
    trap_profile, gamma_model, theta_model, space_charge
        Modelling choices. Defaults reproduce the prototype spreadsheet.
    n_points, V_max, window, bin_size
        Numerical choices. ``window`` and ``bin_size`` belong to the analysis
        branch, the other two to the model branch.
    """

    name: str
    material: Material
    device: Device
    params: ModelParams
    constants: Constants = CODATA
    trap_profile: TrapProfile = DEFAULT_TRAP_PROFILE
    gamma_model: GammaModel = DEFAULT_GAMMA_MODEL
    theta_model: ThetaModel = DEFAULT_THETA_MODEL
    space_charge: SpaceCharge = DEFAULT_SPACE_CHARGE
    n_points: int = 3001
    V_max: float | None = None
    window: int = DEFAULT_SLOPE_WINDOW
    bin_size: int = 1


_CONFIG_SECTIONS: Final[dict[str, tuple[str, ...]]] = {
    "material": ("name", "E_c", "E_v", "eps_r", "m_eff_h", "m_eff_e"),
    "device": ("thickness", "area", "temperature"),
    "params": ("mu_0", "N_t", "E_t", "T_t", "E_F0"),
    "model": (
        "trap_profile", "gamma_model", "theta_model", "space_charge", "constants",
    ),
    "numerics": ("n_points", "V_max", "window", "bin_size"),
}


def _check_keys(section: str, given: dict, allowed: tuple[str, ...]) -> None:
    """Reject unknown keys.

    Silently ignoring an unrecognised key is the wrong failure mode here: a
    misspelled ``E_t`` would leave the model running on its default while the
    file appears to say otherwise, and nothing downstream would look wrong.
    """
    unknown = set(given) - set(allowed)
    if unknown:
        raise ValueError(
            f"[{section}]: unknown key(s) {sorted(unknown)}; "
            f"allowed keys are {list(allowed)}"
        )


def load_config(path: str) -> Configuration:
    """Read a Configuration from a TOML parameter file.

    The ``[material]``, ``[device]`` and ``[params]`` sections are required and
    every key in them is mandatory. ``[model]`` and ``[numerics]`` are optional
    and every key in them defaults to the value that reproduces the reference
    prototype.

    Raises
    ------
    ValueError
        If a section or key is missing, unknown, or names a modelling choice
        that is not one of the documented options.
    """
    try:
        import tomllib
    except ModuleNotFoundError:  # Python 3.10
        import tomli as tomllib  # type: ignore[no-redef]

    with open(path, "rb") as fh:
        raw = tomllib.load(fh)

    _check_keys("top level", {k: v for k, v in raw.items() if isinstance(v, dict)},
                tuple(_CONFIG_SECTIONS))

    for section in ("material", "device", "params"):
        if section not in raw:
            raise ValueError(f"{path}: missing required section [{section}]")
        _check_keys(section, raw[section], _CONFIG_SECTIONS[section])
        missing = set(_CONFIG_SECTIONS[section]) - set(raw[section])
        if missing:
            raise ValueError(f"{path}: [{section}] is missing {sorted(missing)}")

    model = raw.get("model", {})
    numerics = raw.get("numerics", {})
    _check_keys("model", model, _CONFIG_SECTIONS["model"])
    _check_keys("numerics", numerics, _CONFIG_SECTIONS["numerics"])

    def choice(enum, key, default):
        if key not in model:
            return default
        try:
            return enum(model[key])
        except ValueError:
            raise ValueError(
                f"{path}: [model] {key} = {model[key]!r} is not one of "
                f"{[e.value for e in enum]}"
            ) from None

    constants_name = model.get("constants", "codata").lower()
    if constants_name not in ("codata", "prototype"):
        raise ValueError(
            f"{path}: [model] constants = {model['constants']!r}; "
            "expected 'codata' or 'prototype'"
        )

    return Configuration(
        name=str(raw.get("name", raw["material"]["name"])),
        material=Material(**raw["material"]),
        device=Device(**raw["device"]),
        params=ModelParams(**raw["params"]),
        constants=Constants.prototype() if constants_name == "prototype" else CODATA,
        trap_profile=choice(TrapProfile, "trap_profile", DEFAULT_TRAP_PROFILE),
        gamma_model=choice(GammaModel, "gamma_model", DEFAULT_GAMMA_MODEL),
        theta_model=choice(ThetaModel, "theta_model", DEFAULT_THETA_MODEL),
        space_charge=choice(SpaceCharge, "space_charge", DEFAULT_SPACE_CHARGE),
        n_points=int(numerics.get("n_points", 3001)),
        V_max=numerics.get("V_max"),
        window=int(numerics.get("window", DEFAULT_SLOPE_WINDOW)),
        bin_size=int(numerics.get("bin_size", 1)),
    )


# --------------------------------------------------------------------------- #
# Reference parameter sets
# --------------------------------------------------------------------------- #

#: MAPbBr3 single crystal, Supplementary Table S4.
MAPBBR3: Final[Material] = Material(
    name="MAPbBr3", E_c=-3.36, E_v=-5.58, eps_r=25.5, m_eff_h=0.305, m_eff_e=0.32
)

#: MAPbI3 single crystal, Supplementary Table S4.
MAPBI3: Final[Material] = Material(
    name="MAPbI3", E_c=-3.93, E_v=-5.43, eps_r=32.0, m_eff_h=0.35, m_eff_e=0.35
)

#: The configuration found in the prototype spreadsheet: MAPbBr3 sample "S2",
#: dark, measured 2024-09-27. Note that this is neither of the two samples
#: reported in the article, whose parameters are in Table 1.
MAPBBR3_S2: Final[tuple[Material, Device, ModelParams]] = (
    MAPBBR3,
    Device(thickness=6.0e-4, area=7.7e-6, temperature=299.0),
    ModelParams(mu_0=2.7e-3, N_t=4.7e16, E_t=-4.82, T_t=30.0, E_F0=-4.84),
)


# --------------------------------------------------------------------------- #
# Command line entry point
# --------------------------------------------------------------------------- #


def _summarise(curve: ModelCurve, params: ModelParams, material: Material) -> str:
    ok = np.isfinite(curve.V) & np.isfinite(curve.J) & (curve.V > 0)
    lines = [
        f"gamma            {curve.gamma:.6f}",
        f"E_F sweep        {curve.E_F[0]:.4f} -> {curve.E_F[-1]:.4f} eV",
        f"V range          {curve.V[ok].min():.4g} .. {curve.V[ok].max():.4g} V",
        f"J range          {curve.J[ok].min():.4g} .. {curve.J[ok].max():.4g} A/m2",
        (
            f"p_t max          {curve.p_t.max():.6g} m-3  "
            f"({curve.p_t.max() / params.N_t:.4f} N_t)"
        ),
        f"Theta range      {np.nanmin(curve.theta[curve.valid]):.4g} .. "
        f"{np.nanmax(curve.theta[curve.valid]):.4g}"
        + ("" if curve.valid.all()
           else f"   ({(~curve.valid).sum()} low-bias points flagged invalid)"),
        (
            # Through `valid`, like Theta above: mu_eff = mu_0 * Theta, so an
            # unmasked maximum reports a mobility above mu_0 wherever a
            # low-bias point was flagged, which is what the mask exists to
            # exclude.
            f"mu_eff max       {np.nanmax(curve.mu_eff[curve.valid]):.6g} m2/V/s "
            f"(mu_0 = {params.mu_0:.4g})"
        ),
    ]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    """Run the model branch on the prototype spreadsheet configuration."""
    parser = argparse.ArgumentParser(
        description=(
            "A-SCLC model. With no arguments, runs the model branch on the "
            "prototype spreadsheet's MAPbBr3 S2 configuration and prints a summary."
        )
    )
    parser.add_argument(
        "--params",
        type=str,
        default=None,
        metavar="PATH",
        help="TOML parameter file describing the sample and the five model "
        "parameters (see params/). Without it the bundled prototype S2 "
        "configuration is used.",
    )
    parser.add_argument(
        "--trap-profile",
        choices=[p.value for p in TrapProfile],
        default=DEFAULT_TRAP_PROFILE.value,
        help="localized-state functional form (default: %(default)s)",
    )
    parser.add_argument(
        "--gamma-model",
        choices=[g.value for g in GammaModel],
        default=DEFAULT_GAMMA_MODEL.value,
        help="gamma expression for the model branch (default: %(default)s)",
    )
    parser.add_argument(
        "--theta-model",
        choices=[t.value for t in ThetaModel],
        default=DEFAULT_THETA_MODEL.value,
        help="reading of Eq (S12) for Theta (default: %(default)s)",
    )
    parser.add_argument(
        "--prototype-constants",
        action="store_true",
        help="use the rounded constants from the prototype rather than CODATA",
    )
    parser.add_argument(
        "--v-max",
        type=float,
        default=10.0,
        help="drop modelled points above this voltage (default: %(default)s)",
    )
    parser.add_argument(
        "--data",
        type=str,
        default=None,
        help="measurement file to analyse alongside the model (voltage, current"
        " and optionally temperature columns)",
    )
    parser.add_argument(
        "--celsius",
        action="store_true",
        help="the data file's temperature column is in degrees Celsius",
    )
    parser.add_argument(
        "--bin",
        type=int,
        default=1,
        dest="bin_size",
        help="average this many consecutive data points (default: %(default)s)",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=DEFAULT_SLOPE_WINDOW,
        help="points per local slope fit (default: %(default)s)",
    )
    parser.add_argument(
        "--csv", type=str, default=None, help="write the model curve to this CSV path"
    )
    parser.add_argument(
        "--plot",
        type=str,
        default=None,
        metavar="DIR",
        help="write the model and analysis figures into this directory "
        "(needs the 'plot' extra: pip install -e '.[plot]')",
    )
    parser.add_argument(
        "--plot-format",
        default="png",
        help="comma-separated figure formats, e.g. png,svg (default: %(default)s)",
    )
    args = parser.parse_args(argv)

    if args.params:
        config = load_config(args.params)
    else:
        material, device, params = MAPBBR3_S2
        config = Configuration(
            name="MAPbBr3 S2, dark (bundled)",
            material=material,
            device=device,
            params=params,
        )
    material, device, params = config.material, config.device, config.params

    # An explicitly given flag overrides the file; a flag left at its default
    # does not, so a parameter file can set a non-default choice and still be
    # run without repeating it on the command line.
    given = set(argv if argv is not None else sys.argv[1:])

    def flag(name: str, value, from_file):
        return value if any(a.startswith(name) for a in given) else from_file

    constants = (
        Constants.prototype() if args.prototype_constants else config.constants
    )
    trap_profile = flag(
        "--trap-profile", TrapProfile(args.trap_profile), config.trap_profile
    )
    gamma_model = flag(
        "--gamma-model", GammaModel(args.gamma_model), config.gamma_model
    )
    theta_model = flag(
        "--theta-model", ThetaModel(args.theta_model), config.theta_model
    )
    # The file wins when it sets V_max, otherwise --v-max's own default stands.
    # Falling through to a bare None would sweep E_F to the band edge, where the
    # injected charge and the voltage with it run to 1e8 V.
    v_max = (
        config.V_max
        if config.V_max is not None and not any(a.startswith("--v-max") for a in given)
        else args.v_max
    )
    window = flag("--window", args.window, config.window)
    bin_size = flag("--bin", args.bin_size, config.bin_size)

    curve = model_curve(
        params,
        material,
        device,
        constants,
        n_points=config.n_points,
        profile=trap_profile,
        gamma_model=gamma_model,
        theta_model=theta_model,
        space_charge=config.space_charge,
        V_max=v_max,
    )

    print(f"A-SCLC model branch: {config.name}")
    print(f"  material       {material.name}")
    print(f"  trap profile   {trap_profile.value}")
    print(f"  gamma model    {gamma_model.value}")
    print(f"  theta model    {theta_model.value}")
    print(f"  constants      {'prototype' if constants != CODATA else 'CODATA'}")
    print()
    print(_summarise(curve, params, material))

    analysis = None
    if args.data:
        measurement = load_jv(
            args.data,
            device,
            temperature_in_celsius=args.celsius,
            bin_size=bin_size,
        )
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            analysis = analyse_jv(
                measurement.V,
                measurement.J,
                material,
                device,
                mu_0=params.mu_0,
                temperature=measurement.temperature,
                window=window,
            )
        print()
        print(f"measurement: {measurement.source}")
        print(f"  {measurement.n_raw} rows -> {len(measurement)} points")
        v = analysis.valid
        print(f"  usable           {analysis.n_valid} of {len(measurement)}")
        if analysis.n_valid:
            print(f"  V (usable)       {analysis.V[v].min():.4g} .. "
                  f"{analysis.V[v].max():.4g} V")
            print(f"  slope m          {np.nanmin(analysis.m[v]):.3g} .. "
                  f"{np.nanmax(analysis.m[v]):.3g}")
            print(f"  mu_eff           {np.nanmin(analysis.mu_eff[v]):.4g} .. "
                  f"{np.nanmax(analysis.mu_eff[v]):.4g} m2/V/s")
            print(f"  E_F              {np.nanmin(analysis.E_F[v]):.4f} .. "
                  f"{np.nanmax(analysis.E_F[v]):.4f} eV")
        for w in caught:
            print(f"  note: {w.message}")

    if args.csv:
        header = "E_F_eV,V_V,J_A_per_m2,p_f_m-3,p_t_m-3,theta,mu_eff_m2_per_Vs"
        data = np.column_stack(
            [curve.E_F, curve.V, curve.J, curve.p_f, curve.p_t, curve.theta, curve.mu_eff]
        )
        np.savetxt(args.csv, data, delimiter=",", header=header, comments="")
        print(f"\nwrote {data.shape[0]} rows to {args.csv}")

    if args.plot:
        try:
            from asclc_plot import write_figures
        except ImportError as exc:
            print(f"\ncannot plot: {exc}")
            return 1
        written = write_figures(
            args.plot,
            curve,
            params,
            material,
            device,
            analysis,
            constants,
            profile=trap_profile,
            formats=tuple(f.strip() for f in args.plot_format.split(",") if f.strip()),
        )
        print(f"\nwrote {len(written)} figures to {args.plot}/")
        for path in written:
            print(f"  {path}")

    return 0


if __name__ == "__main__":
    # Re-enter through the module name rather than calling main() directly.
    # Running the file as a script binds it as __main__; a sibling module that
    # imports asclc then gets a second, distinct module object, so the enum
    # members are not the same objects and every `profile is TrapProfile.X`
    # comparison silently fails.
    import asclc

    raise SystemExit(asclc.main())
