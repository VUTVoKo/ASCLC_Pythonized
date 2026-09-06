"""Tests for the A-SCLC model.

Run with ``pytest test_asclc.py``.

The tests fall into three groups: analytic identities that must hold whatever
the parameters (normalizations, the Mott-Gurney limit), consistency between two
independent routes to the same quantity (the DOS integral against the effective
density of states), and regression guards on the specific numbers observed in
the prototype spreadsheet.
"""

from __future__ import annotations

import os
import warnings
from dataclasses import replace

import numpy as np
import pytest

from asclc import (
    CODATA,
    MAPBBR3,
    MAPBBR3_S2,
    MAPBI3,
    Constants,
    Device,
    EnergyGrid,
    GammaModel,
    Material,
    ModelParams,
    ThetaModel,
    TrapProfile,
    analyse_jv,
    carrier_densities,
    effective_dos,
    fermi_dirac,
    gamma_from_trap_temperature,
    load_config,
    load_jv,
    local_loglog_slope,
    model_curve,
    trap_dos,
    trap_dos_norm,
    valence_band_dos,
)

MATERIAL, DEVICE, PARAMS = MAPBBR3_S2


# --------------------------------------------------------------------------- #
# Trap distribution normalization
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    ("profile", "expected_factor"),
    [
        (TrapProfile.SI_BIEXPONENTIAL, 1.0),
        (TrapProfile.GAUSSIAN, 1.0),
        (TrapProfile.PROTOTYPE_SECH, np.pi / 4),
    ],
)
def test_trap_dos_normalization(profile: TrapProfile, expected_factor: float) -> None:
    """Each profile integrates to its documented multiple of N_t.

    The pi/4 for the prototype variant is why it exists as a separate option; if
    this passes with 1.0 the sech form has been swapped for the biexponential.
    """
    kT_t = CODATA.kT_eV(PARAMS.T_t)
    E = np.linspace(PARAMS.E_t - 60 * kT_t, PARAMS.E_t + 60 * kT_t, 200_001)
    total = np.trapezoid(trap_dos(E, PARAMS, profile=profile), E)
    assert total == pytest.approx(expected_factor * PARAMS.N_t, rel=1e-6)
    assert total == pytest.approx(trap_dos_norm(PARAMS, profile=profile), rel=1e-6)


def test_trap_profiles_share_peak_height() -> None:
    """The two biexponential-family profiles peak at the same value.

    This is why the two are indistinguishable on a log-scale DOS plot: they
    differ only in the wings.
    """
    E = np.array([PARAMS.E_t])
    si = trap_dos(E, PARAMS, profile=TrapProfile.SI_BIEXPONENTIAL)[0]
    wb = trap_dos(E, PARAMS, profile=TrapProfile.PROTOTYPE_SECH)[0]
    assert si == pytest.approx(wb, rel=1e-12)


def test_trap_dos_is_symmetric_about_Et() -> None:
    kT_t = CODATA.kT_eV(PARAMS.T_t)
    d = np.array([0.5, 1.0, 3.0, 10.0]) * kT_t
    for profile in TrapProfile:
        lo = trap_dos(PARAMS.E_t - d, PARAMS, profile=profile)
        hi = trap_dos(PARAMS.E_t + d, PARAMS, profile=profile)
        # 1e-10 rather than exact: expit is not bitwise symmetric about 0.
        assert lo == pytest.approx(hi, rel=1e-10), profile


def test_trap_dos_no_overflow_far_from_peak() -> None:
    """The sech and biexponential forms must not overflow in the far tails."""
    E = PARAMS.E_t + np.array([-500.0, -50.0, 0.0, 50.0, 500.0])
    for profile in TrapProfile:
        g = trap_dos(E, PARAMS, profile=profile)
        assert np.all(np.isfinite(g)), profile
        assert np.all(g >= 0.0), profile


# --------------------------------------------------------------------------- #
# Band density of states
# --------------------------------------------------------------------------- #


def test_effective_dos_matches_textbook_value() -> None:
    """N_v should reduce to the standard 2.5e25 m^-3 (m*)^{3/2} at 300 K."""
    m_eff = 1.0
    expected = 2.509e25  # m^-3, the usual quoted value for m* = m_0, T = 300 K
    assert effective_dos(m_eff, 300.0) == pytest.approx(expected, rel=2e-3)


def test_dos_integral_reproduces_effective_dos() -> None:
    """Boltzmann limit: the valence-band integral must equal N_v exp(-dE/kT).

    This checks the square-root prefactor of Eq (S2) and the effective density of
    states of Eq (S4) against each other. They are independent expressions, so
    agreement validates both.
    """
    E_F = MAPBBR3.E_v + 0.20  # well above E_v, so non-degenerate
    kT = CODATA.kT_eV(300.0)
    E = np.linspace(MAPBBR3.E_v - 3.0, MAPBBR3.E_v, 400_001)

    g = valence_band_dos(E, MAPBBR3)
    one_minus_f = 1.0 - fermi_dirac(E, E_F, 300.0)
    numeric = np.trapezoid(g * one_minus_f, E)

    N_v = effective_dos(MAPBBR3.m_eff_h, 300.0)
    analytic = N_v * np.exp(-(E_F - MAPBBR3.E_v) / kT)

    assert numeric == pytest.approx(analytic, rel=1e-3)


# --------------------------------------------------------------------------- #
# Occupation integrals
# --------------------------------------------------------------------------- #


def test_carrier_densities_theta_identity() -> None:
    """By default Theta is Eq (S12) as printed: p_f over the absolute total.

    The prototype spreadsheet divides by the injected total instead, which is a
    different quantity by the whole equilibrium population; the two are more
    than 3x apart over this sweep.
    """
    E_F = np.linspace(PARAMS.E_F0 - 0.02, MATERIAL.E_v + 0.05, 51)
    d = carrier_densities(E_F, MATERIAL, DEVICE, PARAMS)
    assert d.p_s == pytest.approx(d.p_f + d.p_t, rel=1e-12)
    assert d.theta_p == pytest.approx(d.p_f / d.p_s, rel=1e-12)

    proto = carrier_densities(
        E_F, MATERIAL, DEVICE, PARAMS,
        theta_model=ThetaModel.ABSOLUTE_OVER_INJECTED,
    )
    assert proto.theta_p == pytest.approx(d.p_f / d.p_s_injected, rel=1e-12)
    assert np.max(proto.theta_p / d.theta_p) > 3.0


def test_injected_total_is_exactly_zero_at_equilibrium() -> None:
    """Exactly zero, not merely small.

    The equilibrium reference is evaluated inside the same array as the
    requested Fermi levels. Computing it in a separate call lets the chunked
    matrix products round it differently, leaving a residual that is sometimes
    negative and breaks any V > 0 filter or log-scale plot downstream.
    """
    E_F = np.linspace(PARAMS.E_F0, MATERIAL.E_v + 0.1, 257)
    d = carrier_densities(E_F, MATERIAL, DEVICE, PARAMS)
    assert d.p_s_injected[0] == 0.0
    assert d.n_s_injected[0] == 0.0


def test_carrier_densities_monotonic_in_EF() -> None:
    """Moving E_F toward the valence band must increase the hole population."""
    E_F = np.linspace(PARAMS.E_F0, MATERIAL.E_v + 0.05, 51)
    d = carrier_densities(E_F, MATERIAL, DEVICE, PARAMS)
    assert np.all(np.diff(d.p_f) > 0)
    assert np.all(np.diff(d.p_t) > 0)
    # electrons go the other way
    assert np.all(np.diff(d.n_f) < 0)


def test_carrier_densities_free_holes_match_boltzmann() -> None:
    """p_f from the DOS integral must match Eq (7) in the non-degenerate range.

    Eq (7) is the Boltzmann approximation to the Fermi-Dirac integral, so the
    two agree only well away from the band edge. At 10 k_B T the correction is
    of order exp(-10) ~ 5e-5; at 6 k_B T it is already 0.3 %, which is why the
    range starts at 0.25 eV rather than at the 3 k_B T the SI mentions.
    """
    kT = CODATA.kT_eV(DEVICE.temperature)
    E_F = np.linspace(MATERIAL.E_v + 0.25, MATERIAL.E_v + 0.50, 25)
    assert (E_F[0] - MATERIAL.E_v) / kT > 9.0, "range must be non-degenerate"

    d = carrier_densities(E_F, MATERIAL, DEVICE, PARAMS)
    N_v = effective_dos(MATERIAL.m_eff_h, DEVICE.temperature)
    boltzmann = N_v * np.exp(-(E_F - MATERIAL.E_v) / kT)
    assert d.p_f == pytest.approx(boltzmann, rel=1e-3)


def test_boltzmann_overestimates_near_band_edge() -> None:
    """Closer to the edge, Eq (7) must overestimate p_f, and by ~exp(-dE/kT).

    Agreement everywhere would mean the Fermi-Dirac factor had been replaced by
    a Boltzmann one somewhere in the occupation integral.
    """
    kT = CODATA.kT_eV(DEVICE.temperature)
    E_F = np.array([MATERIAL.E_v + 6.0 * kT])
    d = carrier_densities(E_F, MATERIAL, DEVICE, PARAMS)
    N_v = effective_dos(MATERIAL.m_eff_h, DEVICE.temperature)
    boltzmann = N_v * np.exp(-(E_F - MATERIAL.E_v) / kT)

    rel_error = float((boltzmann[0] - d.p_f[0]) / boltzmann[0])
    assert rel_error > 0
    assert rel_error == pytest.approx(np.exp(-6.0) / 2**1.5, rel=0.05)


# --------------------------------------------------------------------------- #
# Space-charge factor and the Mott-Gurney limit
# --------------------------------------------------------------------------- #


def test_mott_gurney_prefactor() -> None:
    """At gamma = 1/2 the Eq (2) factor must be exactly 9/8.

    This pins the grouping (1-g)(2-g)^2 rather than (1-g)(2-g)/2, which is how
    the factor reads if the superscript is lost. Eq (3) is the authority.
    """
    gamma = 0.5
    assert (1.0 - gamma) * (2.0 - gamma) ** 2 == pytest.approx(9.0 / 8.0)


def test_pt_prefactor_is_maximal_at_zero_gamma() -> None:
    """(1-g)(2-g) decreases monotonically on [0, 1], peaking at gamma = 0.

    Consequence: a degenerate slope fit returning gamma = 0 does not add noise
    to p_t, it returns p_t's upper envelope. This is why local_loglog_slope
    requires a window of at least 3.
    """
    g = np.linspace(0.0, 1.0, 1001)
    factor = (1.0 - g) * (2.0 - g)
    assert factor[0] == pytest.approx(2.0)
    assert factor[-1] == pytest.approx(0.0)
    assert np.all(np.diff(factor) < 0)


# --------------------------------------------------------------------------- #
# Analysis branch
# --------------------------------------------------------------------------- #


def test_local_slope_on_exact_power_law() -> None:
    """A synthetic J = k V^2 curve must give m = 2 everywhere."""
    V = np.logspace(-2, 1, 200)
    J = 3.7 * V**2
    m = local_loglog_slope(V, J, window=7)
    assert np.allclose(m, 2.0, rtol=1e-10)


@pytest.mark.parametrize("exponent", [1.0, 1.5, 2.0, 4.0])
def test_local_slope_recovers_arbitrary_exponent(exponent: float) -> None:
    V = np.logspace(-2, 1, 300)
    J = 0.5 * V**exponent
    m = local_loglog_slope(V, J, window=5)
    assert np.allclose(m, exponent, rtol=1e-10)


def test_local_slope_rejects_degenerate_window() -> None:
    """A window of 1 must raise rather than silently return zero slopes.

    The prototype's binning cell is set to 0, making every LINEST call a
    one-point regression that returns 0, which propagates an upper-envelope p_t
    into every downstream quantity.
    """
    V = np.logspace(-2, 1, 50)
    J = V**2
    with pytest.raises(ValueError, match="at least 3"):
        local_loglog_slope(V, J, window=1)


def test_local_slope_rejects_even_window() -> None:
    V = np.logspace(-2, 1, 50)
    with pytest.raises(ValueError, match="odd"):
        local_loglog_slope(V, V**2, window=6)


def test_analysis_recovers_mott_gurney_mobility() -> None:
    """Round trip: build a J-V curve from Eq (3), recover mu_0 from Eq (4).

    In the Mott-Gurney regime Theta = 1, so mu_eff must come back equal to the
    mu_0 that generated the curve.
    """
    mu_0 = 5.0e-3
    V = np.logspace(-1, 1, 400)
    J = (
        9.0 / 8.0
        * CODATA.eps_0
        * MATERIAL.eps_r
        * mu_0
        * V**2
        / DEVICE.thickness**3
    )
    res = analyse_jv(V, J, MATERIAL, DEVICE, mu_0=mu_0, window=7)
    interior = slice(10, -10)
    assert np.allclose(res.m[interior], 2.0, rtol=1e-9)
    assert np.allclose(res.mu_eff[interior], mu_0, rtol=1e-9)


def test_analysis_theta_is_bounded() -> None:
    mu_0 = 5.0e-3
    V = np.logspace(-1, 1, 200)
    J = 1e-6 * V**2
    res = analyse_jv(V, J, MATERIAL, DEVICE, mu_0=mu_0)
    finite = np.isfinite(res.theta)
    assert np.all(res.theta[finite] >= 0.0)
    assert np.all(res.theta[finite] <= 1.0)


# --------------------------------------------------------------------------- #
# Model branch
# --------------------------------------------------------------------------- #


def test_model_curve_is_monotonic_in_voltage() -> None:
    """Sweeping E_F toward E_v must raise both p_t and hence V monotonically."""
    curve = model_curve(PARAMS, MATERIAL, DEVICE, n_points=401)
    assert np.all(np.diff(curve.p_t) > 0)
    assert np.all(np.diff(curve.V) > 0)


def test_model_pt_saturates_at_trap_total() -> None:
    """Deep in the sweep p_t must approach the trap normalization.

    With the biexponential profile that is N_t; with the prototype profile it is
    (pi/4) N_t. This is the property Supplementary Note A6 relies on when it
    reads N_t off the saturation of p_t.
    """
    for profile in (TrapProfile.SI_BIEXPONENTIAL, TrapProfile.PROTOTYPE_SECH):
        curve = model_curve(
            PARAMS, MATERIAL, DEVICE, n_points=801, profile=profile
        )
        expected = trap_dos_norm(PARAMS, profile=profile)
        assert curve.p_t.max() == pytest.approx(expected, rel=5e-3), profile


def test_model_theta_approaches_unity_at_high_bias() -> None:
    """Once free carriers dominate the traps, Theta -> 1 and mu_eff -> mu_0."""
    curve = model_curve(PARAMS, MATERIAL, DEVICE, n_points=801)
    assert np.nanmax(curve.theta) > 0.99
    assert np.nanmax(curve.mu_eff) == pytest.approx(PARAMS.mu_0, rel=1e-2)


def test_model_curve_is_pure() -> None:
    """Calling twice with the same inputs gives bitwise identical output."""
    a = model_curve(PARAMS, MATERIAL, DEVICE, n_points=201)
    b = model_curve(PARAMS, MATERIAL, DEVICE, n_points=201)
    assert np.array_equal(a.V, b.V)
    assert np.array_equal(a.J, b.J)


def test_gamma_models_differ_as_documented() -> None:
    """Tt/T = 0.1003, Note M4 gives 0.5, and Tt/(T+Tt) = 0.0912 here."""
    a = gamma_from_trap_temperature(PARAMS, DEVICE, model=GammaModel.TT_OVER_T)
    b = gamma_from_trap_temperature(
        PARAMS, DEVICE, model=GammaModel.TT_OVER_T_PLUS_TT
    )
    c = gamma_from_trap_temperature(PARAMS, DEVICE, model=GammaModel.SI_NOTE_M4)
    assert a == pytest.approx(30.0 / 299.0, rel=1e-12)
    assert b == pytest.approx(30.0 / 329.0, rel=1e-12)
    assert c == 0.5


def test_si_note_m4_gamma_is_one_over_the_mark_helfrich_exponent() -> None:
    """Note M4 as printed: gamma = T/(T_t + T) for T_t >= T, else 0.5.

    T/(T_t + T) is 1/(1 + T_t/T), i.e. 1/m for the trap-filled-limit exponent
    m = 1 + T_t/T of an exponential trap distribution, which is what reconciles
    Note M4 with the article's own definition gamma = 1/m. Spec section 7.2.

    Note that this is NOT Tt/(T + Tt): that expression appears in no source and
    is the complement of this one. Earlier revisions of the port carried it
    mislabelled as the SI's reading.
    """
    T = DEVICE.temperature
    hot = replace(PARAMS, T_t=2.0 * T)  # T_t >= T selects the formula branch
    g = gamma_from_trap_temperature(hot, DEVICE, model=GammaModel.SI_NOTE_M4)

    assert g == pytest.approx(T / (hot.T_t + T), rel=1e-12)
    m = 1.0 + hot.T_t / T
    assert g == pytest.approx(1.0 / m, rel=1e-12)
    # and it is the complement of the expression it was once confused with
    other = gamma_from_trap_temperature(
        hot, DEVICE, model=GammaModel.TT_OVER_T_PLUS_TT
    )
    assert g + other == pytest.approx(1.0, rel=1e-12)


def test_si_note_m4_takes_the_constant_branch_for_every_published_set() -> None:
    """Every configuration in the article and the prototype has T_t < T.

    So Note M4 as printed selects its gamma = 0.5 branch throughout, and cannot
    reproduce the prototype's T_t/T. This is the substance of spec section 7.2.
    """
    configs = [(DEVICE, PARAMS)] + [
        (device, params) for _, device, params in PUBLISHED.values()
    ]
    for device, params in configs:
        assert params.T_t < device.temperature
        g = gamma_from_trap_temperature(
            params, device, model=GammaModel.SI_NOTE_M4
        )
        assert g == 0.5


def test_gamma_warns_when_out_of_range() -> None:
    """Note M4 states gamma <= 0.5; a hot trap breaks that and must warn."""
    hot = ModelParams(mu_0=2.7e-3, N_t=4.7e16, E_t=-4.82, T_t=400.0, E_F0=-4.84)
    with pytest.warns(RuntimeWarning, match="outside"):
        gamma_from_trap_temperature(hot, DEVICE, model=GammaModel.TT_OVER_T)


# --------------------------------------------------------------------------- #
# Input validation
# --------------------------------------------------------------------------- #


def test_material_rejects_inverted_bands() -> None:
    with pytest.raises(ValueError, match="must lie above"):
        Material(name="bad", E_c=-5.0, E_v=-3.0, eps_r=1.0, m_eff_h=1.0, m_eff_e=1.0)


@pytest.mark.parametrize("bad", [{"mu_0": 0.0}, {"N_t": -1.0}, {"T_t": 0.0}])
def test_model_params_reject_nonpositive(bad: dict[str, float]) -> None:
    kwargs = {"mu_0": 1e-3, "N_t": 1e16, "E_t": -4.8, "T_t": 30.0, "E_F0": -4.84}
    kwargs.update(bad)
    with pytest.raises(ValueError, match="must be positive"):
        ModelParams(**kwargs)


def test_device_rejects_nonpositive() -> None:
    with pytest.raises(ValueError, match="must be positive"):
        Device(thickness=-1.0, area=1.0)


# --------------------------------------------------------------------------- #
# Constants
# --------------------------------------------------------------------------- #


def test_prototype_constants_are_close_to_codata() -> None:
    """The prototype's rounded constants agree with CODATA to under 0.1 %."""
    wb = Constants.prototype()
    for name in ("e", "k_B", "h", "eps_0", "m_0"):
        assert getattr(wb, name) == pytest.approx(getattr(CODATA, name), rel=1e-3), name


def test_energy_grid_resolves_the_trap() -> None:
    """The default trap mesh must integrate the trap to its analytic norm."""
    grid = EnergyGrid.build(MATERIAL, PARAMS)
    for profile in (TrapProfile.SI_BIEXPONENTIAL, TrapProfile.PROTOTYPE_SECH):
        g = trap_dos(grid.E_trap, PARAMS, profile=profile)
        total = np.trapezoid(g, grid.E_trap)
        assert total == pytest.approx(
            trap_dos_norm(PARAMS, profile=profile), rel=1e-4
        ), profile


def test_band_quadrature_is_grid_independent() -> None:
    """The sqrt substitution must make p_f insensitive to the point count.

    A uniform mesh converges as O(h^{3/2}) on the band edge and needs ~1e-5 eV
    steps to reach 1e-3; the substituted grid is converged by 1001 points.
    Doubling the resolution must move p_f by less than 1e-9 relative.
    """
    E_F = np.array([MATERIAL.E_v + 0.30])
    values = []
    for n in (1001, 2001, 4001, 8001):
        grid = EnergyGrid.build(MATERIAL, PARAMS, band_points=n)
        d = carrier_densities(E_F, MATERIAL, DEVICE, PARAMS, grid=grid)
        values.append(float(d.p_f[0]))
    spread = (max(values) - min(values)) / np.mean(values)
    assert spread < 1e-9, f"band quadrature not converged, spread = {spread:.2e}"


# --------------------------------------------------------------------------- #
# Space-charge reference
# --------------------------------------------------------------------------- #


def test_voltage_vanishes_at_zero_bias() -> None:
    """V(E_F0) must be exactly zero.

    At E_F = E_F0 no charge has been injected, so the space charge and hence the
    voltage must vanish. The prototype's own model columns cross zero at E_F0,
    which is the property this reproduces.
    """
    from asclc import SpaceCharge

    for sc in (SpaceCharge.INJECTED_TOTAL, SpaceCharge.INJECTED_TRAPPED):
        curve = model_curve(
            PARAMS, MATERIAL, DEVICE, n_points=401, space_charge=sc
        )
        assert curve.V[0] == 0.0, sc
        assert curve.J[0] == 0.0, sc


def test_total_space_charge_does_not_vanish_at_zero_bias() -> None:
    """The rejected alternative must fail the same check, for the record.

    For the prototype's parameters E_t sits above E_F0, so the trap is already
    68 % hole-occupied at equilibrium and the total-charge reading starts at
    several volts, so the two references are not interchangeable.
    """
    from asclc import SpaceCharge

    curve = model_curve(
        PARAMS, MATERIAL, DEVICE, n_points=401, space_charge=SpaceCharge.TOTAL
    )
    assert curve.V[0] > 1.0


def test_space_charge_is_monotonic_and_starts_at_zero() -> None:
    """The driving charge starts at zero and rises through the sweep.

    It is not bounded by p_t: with the default INJECTED_TOTAL reference it
    overtakes the trapped charge once the valence band contributes, which is
    exactly what lets the model reach the measured voltage range.
    """
    curve = model_curve(PARAMS, MATERIAL, DEVICE, n_points=401)
    assert curve.p_space_charge[0] == 0.0
    assert np.all(np.diff(curve.p_space_charge) > 0)
    assert curve.p_space_charge[-1] > curve.p_t[-1]


# --------------------------------------------------------------------------- #
# Regression against the prototype spreadsheet
# --------------------------------------------------------------------------- #
#
# These are the project's only independent numerical oracle: every other check
# is an analytic identity or an internal consistency relation, and those cannot
# catch an assumption the code and its docstrings share. They are kept for that
# reason, not because the prototype is authoritative -- it is not, and the
# defaults no longer follow it. Each therefore selects the prototype's options
# explicitly. docs/prototype_differences.md lists what those options change.

#: Options that reproduce the prototype spreadsheet's occupation integrals.
PROTOTYPE_DENSITIES = {
    "profile": TrapProfile.PROTOTYPE_SECH,
    "theta_model": ThetaModel.ABSOLUTE_OVER_INJECTED,
}

#: Values read from the prototype's own cached cells: E_F (n(E)!A), p_t (n(E)!O)
#: and p_f (n(E)!L), for the S2 configuration. The port must agree to within the
#: prototype's own quadrature error, which is about 2 % on the band integrals
#: because it uses the rectangle rule on a 3 meV grid across a sqrt band edge.
PROTOTYPE_REFERENCE = [
    # E_F (eV),   p_t (m^-3),   p_f (m^-3)
    (-5.199, 1.5794e18, 1.5677e18),
    (-5.100, 4.5235e16, 3.3572e16),
    (-5.001, 1.2342e16, 7.1893e14),
    (-4.950, 1.1520e16, 9.9255e13),
    (-4.899, 1.0013e16, 1.3703e13),
    (-4.860, 5.1740e15, 3.0145e12),
]


def test_injected_charge_matches_reference_values() -> None:
    """p_s - p_s(E_F0) must reproduce the prototype's n(E)!O column.

    This is the check that settled which quantity Eq (S11) means. Trapped-only
    charge fails it by a factor of 140 at the first row.
    """
    E_F = np.array([row[0] for row in PROTOTYPE_REFERENCE])
    expected = np.array([row[1] for row in PROTOTYPE_REFERENCE])

    d = carrier_densities(E_F, MATERIAL, DEVICE, PARAMS, **PROTOTYPE_DENSITIES)
    d0 = carrier_densities(
        np.array([PARAMS.E_F0]), MATERIAL, DEVICE, PARAMS, **PROTOTYPE_DENSITIES
    )
    injected_total = d.p_s - d0.p_s[0]

    assert injected_total == pytest.approx(expected, rel=0.02)


def test_trapped_only_charge_fails_the_reference_comparison() -> None:
    """Trapped-only charge is 140x low here, so it is not what Eq (S11) means."""
    E_F = np.array([PROTOTYPE_REFERENCE[0][0]])
    expected = PROTOTYPE_REFERENCE[0][1]

    d = carrier_densities(E_F, MATERIAL, DEVICE, PARAMS, **PROTOTYPE_DENSITIES)
    d0 = carrier_densities(
        np.array([PARAMS.E_F0]), MATERIAL, DEVICE, PARAMS, **PROTOTYPE_DENSITIES
    )
    injected_trapped = float((d.p_t - d0.p_t[0])[0])

    assert injected_trapped < expected / 100


def test_free_holes_match_reference_values() -> None:
    """p_f must match the prototype's n(E)!L to within its quadrature error."""
    E_F = np.array([row[0] for row in PROTOTYPE_REFERENCE])
    expected = np.array([row[2] for row in PROTOTYPE_REFERENCE])
    d = carrier_densities(E_F, MATERIAL, DEVICE, PARAMS, **PROTOTYPE_DENSITIES)
    assert d.p_f == pytest.approx(expected, rel=0.025)


def test_model_voltage_covers_the_measured_range() -> None:
    """The model must reach the 0-3 V window the J-V curves are measured over."""
    curve = model_curve(PARAMS, MATERIAL, DEVICE, n_points=4001)
    assert curve.V[0] == 0.0
    assert curve.V.max() > 3.0


# --------------------------------------------------------------------------- #
# Eq (7) inversion
# --------------------------------------------------------------------------- #


def test_fermi_level_lands_inside_the_gap() -> None:
    """E_F from Eq (7) must sit in the band gap, not below E_v.

    Eq (7) is p_f = N_v exp(-(E_F - E_v)/kT), so E_F = E_v + kT ln(N_v/p_f).
    Writing ln(p_f/N_v) flips the sign and puts E_F below the valence band,
    which is what the ABS() in the prototype's MODEL!S9 exists to prevent.
    """
    mu_0 = 2.7e-3
    V = np.logspace(-2, 0.5, 200)
    J = 1e-4 * V**2
    res = analyse_jv(V, J, MATERIAL, DEVICE, mu_0=mu_0)
    finite = np.isfinite(res.E_F)
    assert np.all(res.E_F[finite] > MATERIAL.E_v), "E_F fell below the valence band"
    assert np.all(res.E_F[finite] < MATERIAL.E_c)


def test_fermi_level_round_trips_through_eq7() -> None:
    """Feeding p_f back through Eq (7) must return the same p_f."""
    kT = CODATA.kT_eV(DEVICE.temperature)
    N_v = effective_dos(MATERIAL.m_eff_h, DEVICE.temperature)
    p_f = np.array([1e18, 1e20, 1e22])
    E_F = MATERIAL.E_v + kT * np.log(N_v / p_f)
    assert N_v * np.exp(-(E_F - MATERIAL.E_v) / kT) == pytest.approx(p_f, rel=1e-12)


def test_more_holes_moves_fermi_level_toward_the_valence_band() -> None:
    """Physical direction check on Eq (7)."""
    kT = CODATA.kT_eV(DEVICE.temperature)
    N_v = effective_dos(MATERIAL.m_eff_h, DEVICE.temperature)
    E_F = MATERIAL.E_v + kT * np.log(N_v / np.array([1e18, 1e21]))
    assert E_F[1] < E_F[0], "raising p_f must bring E_F closer to E_v"


def test_v_max_truncates_consistently() -> None:
    """V_max must trim every array in the curve to the same length."""
    curve = model_curve(PARAMS, MATERIAL, DEVICE, n_points=2001, V_max=5.0)
    n = curve.V.size
    for name in ("E_F", "J", "p_f", "p_t", "p_space_charge", "theta", "mu_eff"):
        assert getattr(curve, name).size == n, name
    assert curve.V.max() <= 5.0
    assert curve.V[0] == 0.0


def test_v_max_rejects_an_impossible_bound() -> None:
    with pytest.raises(ValueError, match="below the smallest modelled voltage"):
        model_curve(PARAMS, MATERIAL, DEVICE, n_points=201, V_max=-1.0)


def test_theta_and_mobility_match_reference_values() -> None:
    """Theta and mu_eff must track the prototype's n(E)!N and j(U)!K.

    The denominator is the injected total. Using the absolute p_s would include
    the 2.5e16 m^-3 of trapped charge present at zero bias, suppressing Theta by
    up to a factor of six and mu_eff = mu_0 * Theta with it.
    """
    E_F = np.array([-5.199, -5.001, -4.860])
    expected_theta = np.array([0.99262, 0.058246, 0.00058247])

    d = carrier_densities(E_F, MATERIAL, DEVICE, PARAMS, **PROTOTYPE_DENSITIES)
    assert d.theta_p == pytest.approx(expected_theta, rel=0.03)
    assert PARAMS.mu_0 * d.theta_p == pytest.approx(
        PARAMS.mu_0 * expected_theta, rel=0.03
    )


def test_voltage_and_current_match_reference_values() -> None:
    """End-to-end: V and J against the prototype's j(U)!C and j(U)!E."""
    E_F = np.array([-5.199, -5.001, -4.899, -4.860])
    expected_V = np.array([236.05, 1.8446, 1.4965, 0.77331])
    expected_J = np.array([506.79, 0.0018161, 2.8084e-05, 3.1924e-06])

    d = carrier_densities(E_F, MATERIAL, DEVICE, PARAMS, **PROTOTYPE_DENSITIES)
    from asclc import gamma_from_trap_temperature

    gamma = gamma_from_trap_temperature(
        PARAMS, DEVICE, model=GammaModel.TT_OVER_T
    )
    V = (
        CODATA.e
        * DEVICE.thickness**2
        * d.p_s_injected
        / (CODATA.eps_0 * MATERIAL.eps_r * (1 - gamma) * (2 - gamma))
    )
    J = CODATA.e * PARAMS.mu_0 * (2 - gamma) * d.p_f * V / DEVICE.thickness

    assert V == pytest.approx(expected_V, rel=0.02)
    assert J == pytest.approx(expected_J, rel=0.04)


# --------------------------------------------------------------------------- #
# Electron side, degeneracy, and cost
# --------------------------------------------------------------------------- #


def test_theta_n_is_finite_on_a_hole_sweep() -> None:
    """Injecting holes depletes electrons, so n_s_injected is negative.

    Theta is a magnitude ratio, as the prototype's ABS(nf/ns) makes explicit. A
    bare positivity guard made theta_n NaN for every point of a hole sweep,
    silently removing the electron model curves the paper plots in Fig. 4c,d.
    """
    E_F = np.linspace(PARAMS.E_F0 - 0.01, MATERIAL.E_v + 0.1, 25)
    d = carrier_densities(E_F, MATERIAL, DEVICE, PARAMS)
    assert np.all(d.n_s_injected < 0), "electrons should be depleted"
    assert np.all(np.isfinite(d.theta_n))
    assert np.all(d.theta_n >= 0.0)


def test_theta_is_bounded_by_one() -> None:
    E_F = np.linspace(PARAMS.E_F0 - 0.005, MATERIAL.E_v + 0.05, 101)
    d = carrier_densities(E_F, MATERIAL, DEVICE, PARAMS)
    finite = np.isfinite(d.theta_p)
    assert np.all(d.theta_p[finite] <= 1.0 + 1e-9)


def test_analysis_warns_when_degenerate() -> None:
    """Eq (7) is the Boltzmann limit; p_f > N_v means it does not apply.

    The prototype's ABS() form folds such points back into the gap without
    complaint, which hides the failure rather than reporting it.
    """
    V = np.logspace(-2, 1, 50)
    J = 1e6 * V**2
    with pytest.warns(RuntimeWarning, match="degenerate"):
        analyse_jv(V, J, MATERIAL, DEVICE, mu_0=2.7e-3)


def test_model_curve_is_fast_enough_to_tune_against() -> None:
    """Hand fitting means re-evaluating repeatedly, so cost is a feature.

    Generous bound: it exists to catch a grid or memory regression that pushes
    the evaluation out of cache, which costs nearly an order of magnitude.
    """
    import time

    t0 = time.perf_counter()
    model_curve(PARAMS, MATERIAL, DEVICE, n_points=3001, V_max=10.0)
    assert time.perf_counter() - t0 < 1.5


def test_chunking_does_not_change_results() -> None:
    """Results must not depend on the internal chunk size."""
    import asclc

    E_F = np.linspace(PARAMS.E_F0, MATERIAL.E_v + 0.1, 300)
    ref = carrier_densities(E_F, MATERIAL, DEVICE, PARAMS)
    original = asclc._MATRIX_BUDGET
    try:
        asclc._MATRIX_BUDGET = 5000  # forces many small chunks
        small = carrier_densities(E_F, MATERIAL, DEVICE, PARAMS)
    finally:
        asclc._MATRIX_BUDGET = original
    assert small.p_f == pytest.approx(ref.p_f, rel=1e-12)
    assert small.p_s_injected == pytest.approx(ref.p_s_injected, rel=1e-12)


# --------------------------------------------------------------------------- #
# Theta in the analysis branch, and the gamma degeneracy
# --------------------------------------------------------------------------- #


def test_analysis_theta_equals_pf_over_pt_not_pf_over_ps() -> None:
    """Eq (4) fixes Theta = mu_eff/mu_0, which is identically p_f/p_t.

    Eqs (4), (5) and (6) make mu_eff/mu_0 and p_f/p_t algebraically equal for
    every gamma. They are NOT equal to p_f/(p_f + p_t): Eq (6) returns the total
    space charge despite its label, so dividing by p_f + p_t double-counts the
    free carriers. The error grows with Theta, reaching 2.3x as Theta -> 1.
    """
    mu_0 = 2.7e-3
    V = np.logspace(-1, 0.8, 120)
    J = 3e-3 * V**2.4
    r = analyse_jv(V, J, MATERIAL, DEVICE, mu_0=mu_0, window=7)

    good = np.isfinite(r.theta) & np.isfinite(r.p_t) & (r.p_t != 0)
    assert r.theta[good] == pytest.approx(r.mu_eff[good] / mu_0, rel=1e-12)
    assert r.theta[good] == pytest.approx(r.p_f[good] / r.p_t[good], rel=1e-9)

    naive = r.p_f[good] / (r.p_f[good] + r.p_t[good])
    assert np.max(np.abs(r.theta[good] / naive - 1.0)) > 0.1


def test_analysis_warns_when_theta_exceeds_one() -> None:
    """Theta > 1 means mu_eff > mu_0, which cannot happen: mu_0 is too small."""
    V = np.logspace(-1, 1, 80)
    J = 1e3 * V**2
    with pytest.warns(RuntimeWarning, match="unphysical"):
        analyse_jv(V, J, MATERIAL, DEVICE, mu_0=1e-12)


def test_gamma_only_translates_the_model_curve() -> None:
    """gamma rescales V and J by constants; it does not change the shape.

    V carries 1/[(1-g)(2-g)] and J carries (2-g) times V, both independent of
    E_F, so on a log-log plot changing gamma is a rigid translation and the
    local slope is untouched. This is why the choice between the two candidate
    gamma expressions cannot be settled by fitting a single J-V curve: the 1.5 %
    it moves V is degenerate with mu_0 and with the other scale parameters.
    """
    from asclc import GammaModel

    a = model_curve(
        PARAMS, MATERIAL, DEVICE, n_points=2001, V_max=50.0,
        gamma_model=GammaModel.TT_OVER_T,
    )
    b = model_curve(
        PARAMS, MATERIAL, DEVICE, n_points=2001, V_max=50.0,
        gamma_model=GammaModel.TT_OVER_T_PLUS_TT,
    )
    n = min(a.V.size, b.V.size)
    ok = (a.V[:n] > 1e-6) & (b.V[:n] > 1e-6)
    assert a.gamma != b.gamma

    v_ratio = b.V[:n][ok] / a.V[:n][ok]
    j_ratio = b.J[:n][ok] / a.J[:n][ok]
    assert np.ptp(v_ratio) < 1e-12 * np.mean(v_ratio), "V ratio must be constant"
    assert np.ptp(j_ratio) < 1e-12 * np.mean(j_ratio), "J ratio must be constant"

    ma = local_loglog_slope(a.V[:n][ok], a.J[:n][ok], window=11)
    mb = local_loglog_slope(b.V[:n][ok], b.J[:n][ok], window=11)
    assert np.nanmax(np.abs(ma - mb)) < 1e-9


def test_model_gamma_disagrees_with_its_own_curve_slope() -> None:
    """Recorded, not asserted as correct: the model is not self-consistent in gamma.

    Branch A defines gamma = 1/m from the local slope. Branch M assumes a
    constant gamma = T_t/T. The curve branch M produces has a local slope
    implying a quite different, and varying, gamma. Since gamma only translates
    the curve this does not affect its shape, but it does mean the two branches
    use the symbol differently. Flagged in ASCLC_spec.md section 6.2.
    """
    curve = model_curve(PARAMS, MATERIAL, DEVICE, n_points=4001, V_max=5.0)
    ok = curve.V > 1e-3
    m = local_loglog_slope(curve.V[ok], curve.J[ok], window=11)
    sel = slice(len(m) // 4, 3 * len(m) // 4)
    implied = 1.0 / np.nanmedian(m[sel])
    assert not np.isclose(implied, curve.gamma, rtol=0.3)


# --------------------------------------------------------------------------- #
# Per-point temperature in Eq (7)
# --------------------------------------------------------------------------- #


def test_analysis_accepts_per_point_temperature() -> None:
    """Eq (7) must be able to use the measured temperature of each point.

    The prototype spreadsheet does this, taking T from its own recorded column.
    Over its 282-314 K range the extracted E_F moves by 0.075 eV, which is
    larger than either Fermi level shift the article reports (0.046 and
    0.006 eV). Holding T at its nominal value attributes that scatter to the
    physics instead.
    """
    V = np.logspace(-1, 0.4, 40)
    J = 3e-4 * V**1.8
    rng = np.random.default_rng(0)
    T = 299.0 + rng.uniform(-17.0, 15.0, size=V.size)

    fixed = analyse_jv(V, J, MATERIAL, DEVICE, mu_0=2.7e-3)
    varying = analyse_jv(V, J, MATERIAL, DEVICE, mu_0=2.7e-3, temperature=T)

    # everything except E_F is temperature-independent
    assert varying.p_f == pytest.approx(fixed.p_f, rel=1e-12)
    assert varying.p_t == pytest.approx(fixed.p_t, rel=1e-12)
    assert varying.mu_eff == pytest.approx(fixed.mu_eff, rel=1e-12)
    # E_F moves, and by an amount comparable to the reported shifts
    spread = float(np.nanmax(np.abs(varying.E_F - fixed.E_F)))
    assert spread > 0.01, f"temperature should move E_F appreciably, got {spread}"


def test_scalar_temperature_matches_device_default() -> None:
    V = np.logspace(-1, 0.4, 30)
    J = 3e-4 * V**1.8
    a = analyse_jv(V, J, MATERIAL, DEVICE, mu_0=2.7e-3)
    b = analyse_jv(
        V, J, MATERIAL, DEVICE, mu_0=2.7e-3, temperature=DEVICE.temperature
    )
    assert a.E_F == pytest.approx(b.E_F, rel=1e-12)


def test_temperature_shape_is_validated() -> None:
    V = np.logspace(-1, 0.4, 30)
    J = 3e-4 * V**1.8
    with pytest.raises(ValueError, match="scalar or match V"):
        analyse_jv(V, J, MATERIAL, DEVICE, mu_0=2.7e-3, temperature=np.ones(7))
    with pytest.raises(ValueError, match="must be positive"):
        analyse_jv(V, J, MATERIAL, DEVICE, mu_0=2.7e-3, temperature=-5.0)


def test_model_ptm_is_the_space_charge_not_pt() -> None:
    """The article's ptm is p_space_charge; p_t is a different quantity.

    Fig. 4c,d and Fig. 6 plot ptm, which is the prototype's MODEL!AJ, reading
    from n(E)!O, which is the injected total. Reproducing those figures from
    curve.p_t gives a visibly different curve: p_t carries the equilibrium
    trapped population and misses the band contribution at high injection.
    """
    curve = model_curve(PARAMS, MATERIAL, DEVICE, n_points=1501, V_max=5.0)
    assert curve.p_t[0] > 0.5 * PARAMS.N_t, "p_t starts at its equilibrium value"
    assert curve.p_space_charge[0] == 0.0, "ptm starts at zero"
    ratio = curve.p_space_charge[-1] / curve.p_t[-1]
    assert not np.isclose(ratio, 1.0, rtol=0.05), "the two must not be conflated"


# --------------------------------------------------------------------------- #
# Electron side and charge conservation
# --------------------------------------------------------------------------- #


def test_injected_charge_conserves_exactly() -> None:
    """Every injected hole is an electron removed from the same states.

    n_s_injected must be the exact negative of p_s_injected, bit for bit. It is
    not computed as (n_f + n_t) minus its equilibrium value: that omits the
    valence-band term, which for these parameters is 136 times larger than what
    remains.
    """
    E_F = np.linspace(PARAMS.E_F0, MATERIAL.E_v + 0.05, 64)
    d = carrier_densities(E_F, MATERIAL, DEVICE, PARAMS)
    assert np.array_equal(d.n_s_injected, -d.p_s_injected)


def test_electron_densities_match_reference_values() -> None:
    """n_f against n(E)!E, and the injected electron total against n(E)!H."""
    E_F = np.array([-5.199, -5.001, -4.899, -4.860])
    wb_nf = np.array([4.3899e-07, 0.00095727, 0.050224, 0.2283])
    wb_nt = np.array([-1.5794e18, -1.2348e16, -1.0016e16, -5.1765e15])

    d = carrier_densities(E_F, MATERIAL, DEVICE, PARAMS, **PROTOTYPE_DENSITIES)
    assert d.n_f == pytest.approx(wb_nf, rel=0.04)
    assert d.n_s_injected == pytest.approx(wb_nt, rel=0.02)


def test_injected_charge_is_well_conditioned() -> None:
    """The occupancy difference is taken inside the integral, not after it.

    Integrating each band's absolute occupancy and subtracting afterwards
    differences two far-band numbers of order 1e28 whose ulp is near 1e12,
    comparable to the answer itself, so the result then depends on the chunk
    size. This checks the well-conditioned formulation is still in place.
    """
    import asclc

    E_F = np.linspace(PARAMS.E_F0, MATERIAL.E_v + 0.1, 300)
    ref = carrier_densities(E_F, MATERIAL, DEVICE, PARAMS)
    original = asclc._MATRIX_BUDGET
    try:
        asclc._MATRIX_BUDGET = 5000
        small = carrier_densities(E_F, MATERIAL, DEVICE, PARAMS)
    finally:
        asclc._MATRIX_BUDGET = original

    rel = np.abs(small.p_s_injected[1:] / ref.p_s_injected[1:] - 1.0)
    assert np.nanmax(rel) < 1e-12


# --------------------------------------------------------------------------- #
# Measurement loading
# --------------------------------------------------------------------------- #


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return str(p)


def test_load_jv_converts_current_to_density(tmp_path) -> None:
    path = _write(
        tmp_path,
        "iv.csv",
        "U,I,T\n0.10,1.0e-9,25.0\n0.20,2.0e-9,25.5\n0.30,3.0e-9,26.0\n",
    )
    m = load_jv(path, DEVICE, temperature_in_celsius=True)
    assert len(m) == 3
    assert m.V == pytest.approx([0.10, 0.20, 0.30])
    assert m.J == pytest.approx(np.array([1e-9, 2e-9, 3e-9]) / DEVICE.area)
    assert m.temperature == pytest.approx([298.15, 298.65, 299.15])


def test_load_jv_feeds_analyse_jv(tmp_path) -> None:
    """The loader's output must go straight into the analysis branch."""
    V = np.logspace(-1, 0.4, 40)
    J = 5e-4 * V**1.8
    area = DEVICE.area
    rows = "\n".join(f"{v:.8e}\t{j * area:.8e}\t{26.0:.2f}" for v, j in zip(V, J))
    path = _write(tmp_path, "iv.txt", "U(V)\tI(A)\tT(C)\n" + rows + "\n")

    m = load_jv(path, DEVICE, temperature_in_celsius=True)
    assert m.J == pytest.approx(J, rel=1e-6)

    res = analyse_jv(
        m.V, m.J, MATERIAL, DEVICE, mu_0=2.7e-3, temperature=m.temperature
    )
    interior = slice(5, -5)
    assert np.allclose(res.m[interior], 1.8, rtol=1e-6)


def test_load_jv_bins_and_reports_raw_count(tmp_path) -> None:
    rows = "\n".join(f"{i * 0.1:.4f},{i * 1e-9:.4e}" for i in range(1, 13))
    path = _write(tmp_path, "iv.csv", rows + "\n")
    m = load_jv(path, DEVICE, bin_size=3)
    assert m.n_raw == 12
    assert len(m) == 4
    assert m.V[0] == pytest.approx(np.mean([0.1, 0.2, 0.3]))


def test_load_jv_handles_headers_comments_and_delimiters(tmp_path) -> None:
    path = _write(
        tmp_path,
        "iv.dat",
        "# exported 2024-09-27\n\nU\tI\n0.1\t1e-9\n0.2\t2e-9\n",
    )
    m = load_jv(path, DEVICE)
    assert len(m) == 2
    assert m.temperature is None


def test_load_jv_accepts_density_directly(tmp_path) -> None:
    path = _write(tmp_path, "jv.csv", "0.1,3.0\n0.2,6.0\n")
    m = load_jv(path, DEVICE, current_is_density=True)
    assert m.J == pytest.approx([3.0, 6.0])


def test_load_jv_applies_voltage_offset(tmp_path) -> None:
    path = _write(tmp_path, "iv.csv", "0.1,1e-9\n0.2,2e-9\n")
    m = load_jv(path, DEVICE, voltage_offset=0.05)
    assert m.V == pytest.approx([0.15, 0.25])


def test_load_jv_keeps_negative_currents(tmp_path) -> None:
    """Ohmic-region noise straddles zero; dropping those rows would bias it."""
    path = _write(tmp_path, "iv.csv", "0.1,-1e-12\n0.2,2e-9\n0.3,-3e-12\n")
    m = load_jv(path, DEVICE)
    assert len(m) == 3
    assert np.sum(m.J < 0) == 2


@pytest.mark.parametrize(
    ("text", "kwargs", "match"),
    [
        ("0.1\n0.2\n", {}, "at least voltage and current"),
        ("", {}, "no data rows"),
        ("0.1,1e-9\n", {"bin_size": 0}, "bin_size must be at least 1"),
        ("0.1,1e-9,-300\n", {}, "non-positive temperature"),
        ("0.1,1e-9\n", {"columns": (0, 5)}, "beyond the"),
    ],
)
def test_load_jv_rejects_bad_input(tmp_path, text, kwargs, match) -> None:
    path = _write(tmp_path, "bad.csv", text)
    with pytest.raises(ValueError, match=match):
        load_jv(path, DEVICE, **kwargs)


# --------------------------------------------------------------------------- #
# Validity of extracted points
# --------------------------------------------------------------------------- #


def test_valid_mask_rejects_sign_changes_and_out_of_domain_gamma() -> None:
    """A window spanning a sign change in J has no meaningful slope."""
    V = np.logspace(-1, 0.6, 60)
    J = 2e-4 * V**2.0
    J[10:13] *= -1.0  # a noisy patch straddling zero

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = analyse_jv(V, J, MATERIAL, DEVICE, mu_0=2.7e-3, window=7)

    assert not r.valid[9:14].any(), "points near a sign change must be rejected"
    assert r.valid.any(), "clean points must survive"
    assert np.all((r.gamma[r.valid] > 0) & (r.gamma[r.valid] < 1))


def test_valid_mask_keeps_physical_quantities_physical() -> None:
    """Through the mask, Theta and the concentrations stay in range.

    Ungated, the ohmic region produces negative Theta and negative
    concentrations, which cannot be distinguished from signal on a plot.
    """
    V = np.logspace(-2, 0.5, 120)
    rng = np.random.default_rng(1)
    J = 3e-4 * V**2.2 + rng.normal(0.0, 2e-6, size=V.size)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = analyse_jv(V, J, MATERIAL, DEVICE, mu_0=2.7e-3, window=7)

    v = r.valid
    assert v.any()
    assert np.all(r.p_t[v] > 0)
    assert np.all(r.theta[v] > 0)
    assert r.n_valid == int(np.count_nonzero(v))


def test_analysis_warns_when_most_points_are_unusable() -> None:
    V = np.logspace(-2, 0.0, 80)
    rng = np.random.default_rng(2)
    J = rng.normal(0.0, 1e-9, size=V.size)  # pure noise
    with pytest.warns(RuntimeWarning, match="usable|no points survived"):
        analyse_jv(V, J, MATERIAL, DEVICE, mu_0=2.7e-3, window=7)


def test_degeneracy_check_handles_per_point_temperature() -> None:
    """N_v becomes an array when temperature does; the check must align."""
    V = np.logspace(-1, 0.5, 40)
    J = 3e-4 * V**2.0
    T = 299.0 + np.linspace(-15.0, 15.0, V.size)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = analyse_jv(V, J, MATERIAL, DEVICE, mu_0=2.7e-3, temperature=T)
    assert r.valid.any()
    assert np.isfinite(r.E_F[r.valid]).all()


# --------------------------------------------------------------------------- #
# Published parameter sets
# --------------------------------------------------------------------------- #

CM2 = 1e-4   # cm^2 V^-1 s^-1 -> m^2 V^-1 s^-1
CM3 = 1e6    # cm^-3 -> m^-3

#: The four configurations of Supplementary Table S4.
PUBLISHED = {
    "MAPbBr3 dark": (
        MAPBBR3, Device(0.80e-3, 7.70e-6, 300.0),
        ModelParams(mu_0=50 * CM2, N_t=3.03e10 * CM3, E_t=-4.768, T_t=8.0,
                    E_F0=-4.828)),
    "MAPbBr3 light": (
        MAPBBR3, Device(0.80e-3, 7.70e-6, 300.0),
        ModelParams(mu_0=160 * CM2, N_t=3.54e10 * CM3, E_t=-4.768, T_t=40.0,
                    E_F0=-4.874)),
    "MAPbI3 dark": (
        MAPBI3, Device(1.05e-3, 20.10e-6, 300.0),
        ModelParams(mu_0=90 * CM2, N_t=4.60e10 * CM3, E_t=-4.715, T_t=40.0,
                    E_F0=-4.740)),
    "MAPbI3 light": (
        MAPBI3, Device(1.05e-3, 20.10e-6, 300.0),
        ModelParams(mu_0=230 * CM2, N_t=3.45e10 * CM3, E_t=-4.715, T_t=95.0,
                    E_F0=-4.746)),
}


@pytest.mark.parametrize("name", list(PUBLISHED))
def test_published_parameters_span_the_measured_voltage_range(name) -> None:
    """Each Table S4 set must produce a curve covering the measured 0-3 V."""
    material, device, params = PUBLISHED[name]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c = model_curve(params, material, device, n_points=2001, V_max=3.0)
    usable = c.V > 1e-6
    assert c.V[usable].max() > 2.5, name
    assert c.V[0] == 0.0, name
    assert np.all(np.isfinite(c.J[usable])), name


@pytest.mark.parametrize("name", list(PUBLISHED))
def test_theta_is_bounded_within_the_domain(name) -> None:
    """Theta may exceed 1 only where valid is false.

    mu_eff = mu_0 * Theta cannot exceed the microscopic mobility, so any point
    reporting otherwise is outside the space-charge picture and must be flagged.
    """
    material, device, params = PUBLISHED[name]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c = model_curve(params, material, device, n_points=2001, V_max=3.0)
    good = c.valid & np.isfinite(c.theta)
    assert good.any(), name
    assert np.all(c.theta[good] <= 1.0), name
    assert np.all(c.mu_eff[good] <= params.mu_0 * (1.0 + 1e-9)), name


@pytest.mark.parametrize("n_points", [1001, 2001, 3001, 4001])
def test_out_of_domain_points_sit_at_low_bias(n_points) -> None:
    """Where the mixed reference fails, it fails only near zero bias.

    Whether any grid point lands in that region depends on the sampling, so this
    asserts the location of such points rather than their existence.
    """
    material, device, params = PUBLISHED["MAPbBr3 light"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c = model_curve(params, material, device, n_points=n_points, V_max=3.0)
    bad = ~c.valid
    if bad.any():
        assert c.V[bad].max() < 0.05 * c.V.max()


def test_theta_readings_are_all_reachable_and_agree_at_high_injection() -> None:
    """The three readings of Eq (S12) are selectable and differ where it matters.

    They converge once the injected charge swamps the equilibrium population,
    and separate by a large factor through the trap-filling region. Spec 7.3.
    """
    E_F = np.linspace(PARAMS.E_F0 - 0.02, MATERIAL.E_v + 0.05, 201)
    got = {
        tm: carrier_densities(
            E_F, MATERIAL, DEVICE, PARAMS, theta_model=tm
        ).theta_p
        for tm in ThetaModel
    }
    for tm, theta in got.items():
        assert np.all(np.isfinite(theta)), tm

    # deep in injection every reading agrees
    for tm, theta in got.items():
        assert theta[-1] == pytest.approx(
            got[ThetaModel.ABSOLUTE_OVER_INJECTED][-1], rel=1e-3
        ), tm

    # through the trap-filling region they do not
    spread = max(np.max(a / b) for a in got.values() for b in got.values())
    assert spread > 3.0


def test_theta_model_default_is_the_printed_equation() -> None:
    """The default is Eq (S12) as printed in both the SI and the article."""
    E_F = np.array([-5.199, -5.001, -4.860])
    d = carrier_densities(E_F, MATERIAL, DEVICE, PARAMS)
    explicit = carrier_densities(
        E_F, MATERIAL, DEVICE, PARAMS, theta_model=ThetaModel.ABSOLUTE_OVER_TOTAL
    )
    assert d.theta_model is ThetaModel.ABSOLUTE_OVER_TOTAL
    assert np.array_equal(d.theta_p, explicit.theta_p)
    assert d.theta_p == pytest.approx(d.p_f / d.p_s, rel=1e-12)


@pytest.mark.parametrize(
    "theta_model",
    [ThetaModel.ABSOLUTE_OVER_TOTAL, ThetaModel.INJECTED_OVER_INJECTED],
)
def test_alternative_theta_readings_are_bounded_by_one(theta_model) -> None:
    """Only the prototype's mixed-reference reading can exceed 1.

    The other two put a population over a superset of itself under the same
    reference, so they are fractions by construction and `valid` never fires.
    """
    material, device, params = PUBLISHED["MAPbBr3 light"]
    E_F = params.E_F0 - np.logspace(-6, -0.5, 200)
    d = carrier_densities(
        E_F, material, device, params, theta_model=theta_model
    )
    finite = np.isfinite(d.theta_p)
    assert np.all(d.theta_p[finite] <= 1.0 + 1e-12)
    assert d.valid.all()

    # the prototype reading on the same sweep does go out of domain
    ref = carrier_densities(
        E_F, material, device, params,
        theta_model=ThetaModel.ABSOLUTE_OVER_INJECTED,
    )
    assert np.nanmax(ref.theta_p) > 1.0


def test_injected_free_carriers_are_returned() -> None:
    """p_f_injected is p_f referenced to equilibrium, and it is what the

    INJECTED_OVER_INJECTED reading divides by p_s_injected. It used to be
    computed and discarded.
    """
    E_F = np.linspace(PARAMS.E_F0, MATERIAL.E_v + 0.05, 64)
    d = carrier_densities(
        E_F, MATERIAL, DEVICE, PARAMS,
        theta_model=ThetaModel.INJECTED_OVER_INJECTED,
    )
    assert d.p_f_injected[0] == 0.0
    assert d.p_f_injected == pytest.approx(d.p_f - d.p_f[0], rel=1e-9, abs=1e-6)
    nonzero = d.p_s_injected != 0
    assert d.theta_p[nonzero] == pytest.approx(
        d.p_f_injected[nonzero] / d.p_s_injected[nonzero], rel=1e-12
    )


def test_model_curve_carries_the_theta_model_through_v_max() -> None:
    """The selected reading must survive the V_max trim, arrays and all."""
    curve = model_curve(
        PARAMS, MATERIAL, DEVICE, n_points=2001, V_max=5.0,
        theta_model=ThetaModel.ABSOLUTE_OVER_TOTAL,
    )
    assert curve.theta_model is ThetaModel.ABSOLUTE_OVER_TOTAL
    n = curve.V.size
    for name in ("E_F", "J", "p_f", "p_t", "p_space_charge", "theta", "mu_eff"):
        assert getattr(curve, name).size == n, name
    assert np.all(curve.theta[np.isfinite(curve.theta)] <= 1.0)


def test_gamma_estimator_is_not_the_prototypes_under_noise() -> None:
    """1/slope(lnJ|lnV) is not slope(lnV|lnJ) once the data has scatter.

    The article defines them as equal and they are for a clean power law, but
    ordinary least squares is not symmetric. The prototype regresses the other
    way (Data-calculations!K = LINEST(ln V, ln j)); this port fits in the J
    direction. Documented in local_loglog_slope and spec section 5.
    """
    V = np.logspace(-1, 0.8, 200)
    clean = 3e-3 * V**2.4
    rng = np.random.default_rng(3)
    noisy = clean * np.exp(rng.normal(0.0, 0.25, size=V.size))

    def prototype_gamma(v, j, window):
        x, y = np.log(np.abs(j)), np.log(np.abs(v))  # LINEST(ln V, ln j)
        half, out = window // 2, np.full(v.size, np.nan)
        for i in range(v.size):
            lo, hi = max(0, i - half), min(v.size, i + half + 1)
            out[i] = np.polyfit(x[lo:hi], y[lo:hi], 1)[0]
        return out

    interior = slice(10, -10)
    # on the clean power law the two agree to rounding
    ours = 1.0 / local_loglog_slope(V, clean, window=11)
    theirs = prototype_gamma(V, clean, 11)
    assert ours[interior] == pytest.approx(theirs[interior], rel=1e-9)

    # with scatter they do not, and the gap is not negligible
    ours_n = 1.0 / local_loglog_slope(V, noisy, window=11)
    theirs_n = prototype_gamma(V, noisy, 11)
    assert np.max(np.abs(ours_n[interior] / theirs_n[interior] - 1.0)) > 0.05


def test_theta_exceeds_one_immediately_above_equilibrium() -> None:
    """Sampling close enough to E_F0 must reach the out-of-domain region.

    There the injected charge is still smaller than the equilibrium free
    population, so absolute p_f over injected total exceeds 1.
    """
    material, device, params = PUBLISHED["MAPbBr3 light"]
    E_F = params.E_F0 - np.logspace(-6, -3, 40)
    d = carrier_densities(
        E_F, material, device, params,
        theta_model=ThetaModel.ABSOLUTE_OVER_INJECTED,
    )
    assert np.nanmax(d.theta_p) > 1.0
    assert not d.valid.all()
    assert np.all(d.theta_p[d.valid] <= 1.0)


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #

plt_mod = pytest.importorskip("matplotlib", reason="the 'plot' extra is not installed")


def _curve_and_analysis():
    curve = model_curve(PARAMS, MATERIAL, DEVICE, n_points=401, V_max=5.0)
    V = np.logspace(-1, 0.4, 60)
    J = 3e-4 * V**1.8
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        analysis = analyse_jv(V, J, MATERIAL, DEVICE, mu_0=PARAMS.mu_0)
    return curve, analysis


def test_write_figures_produces_every_figure(tmp_path) -> None:
    import asclc_plot

    curve, analysis = _curve_and_analysis()
    written = asclc_plot.write_figures(
        str(tmp_path), curve, PARAMS, MATERIAL, DEVICE, analysis
    )
    assert len(written) == len(asclc_plot.FIGURES)
    for path in written:
        assert os.path.getsize(path) > 0, path


def test_figures_work_without_a_measurement(tmp_path) -> None:
    """The model branch alone must plot: there is not always a measurement."""
    import asclc_plot

    curve = model_curve(PARAMS, MATERIAL, DEVICE, n_points=401, V_max=5.0)
    written = asclc_plot.write_figures(
        str(tmp_path), curve, PARAMS, MATERIAL, DEVICE, None
    )
    assert len(written) == len(asclc_plot.FIGURES)


def test_drawable_is_axis_aware() -> None:
    """Positivity is required only on log axes.

    The bandgap map's abscissa is an energy, negative on the vacuum scale.
    Requiring every array to be positive drops all of its points and the figure
    comes out empty but for the density of states.
    """
    from matplotlib.figure import Figure

    import asclc_plot

    ax = Figure().add_subplot()
    ax.set_yscale("log")
    E = np.array([-5.1, -5.0, -4.9])
    p = np.array([1e16, 1e17, 1e18])
    assert asclc_plot._drawable(ax, E, p).all(), "negative energies must survive"

    ax.set_xscale("log")
    assert not asclc_plot._drawable(ax, E, p).any(), "log x must reject negatives"


def test_figures_separate_valid_from_rejected() -> None:
    """Rejected points are drawn faded, not dropped: both sets must appear."""
    import asclc_plot

    curve, _ = _curve_and_analysis()
    V = np.logspace(-1, 0.6, 60)
    J = 2e-4 * V**2.0
    J[10:13] *= -1.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        noisy = analyse_jv(V, J, MATERIAL, DEVICE, mu_0=PARAMS.mu_0, window=7)
    assert not noisy.valid.all(), "fixture must contain rejected points"

    fig = asclc_plot.plot_jv(curve, noisy)
    ax = fig.axes[0]
    # model line, rejected markers, accepted markers
    assert len(ax.lines) >= 3


# --------------------------------------------------------------------------- #
# Parameter files
# --------------------------------------------------------------------------- #

PARAM_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "params")


def test_bundled_param_file_round_trips() -> None:
    """params/MAPbBr3_S2.toml must reproduce MAPBBR3_S2 exactly.

    The bundled constant and the shipped file describe the same configuration;
    if they drift apart the regression suite and the CLI stop testing the same
    thing.
    """
    cfg = load_config(os.path.join(PARAM_DIR, "MAPbBr3_S2.toml"))
    assert cfg.material == MATERIAL
    assert cfg.device == DEVICE
    assert cfg.params == PARAMS


@pytest.mark.parametrize(
    "name", ["MAPbBr3_dark", "MAPbBr3_light", "MAPbI3_dark", "MAPbI3_light"]
)
def test_published_param_files_match_table_s4(name: str) -> None:
    """The shipped files must carry the article's own published parameters."""
    key = name.replace("_", " ").replace("dark", "dark").replace("light", "light")
    material, device, params = PUBLISHED[key]
    cfg = load_config(os.path.join(PARAM_DIR, f"{name}.toml"))
    assert cfg.material == material
    assert cfg.device == device
    assert cfg.params == params


def test_param_file_defaults_reproduce_the_prototype() -> None:
    """Omitting [model] must select every prototype-reproducing default."""
    from asclc import (
        DEFAULT_GAMMA_MODEL,
        DEFAULT_SPACE_CHARGE,
        DEFAULT_THETA_MODEL,
        DEFAULT_TRAP_PROFILE,
    )

    cfg = load_config(os.path.join(PARAM_DIR, "MAPbBr3_S2.toml"))
    assert cfg.trap_profile is DEFAULT_TRAP_PROFILE
    assert cfg.gamma_model is DEFAULT_GAMMA_MODEL
    assert cfg.theta_model is DEFAULT_THETA_MODEL
    assert cfg.space_charge is DEFAULT_SPACE_CHARGE
    assert cfg.constants == CODATA


_MINIMAL = """
[material]
name = "X"
E_c = -3.0
E_v = -5.0
eps_r = 25.0
m_eff_h = 0.3
m_eff_e = 0.3

[device]
thickness = 6.0e-4
area = 7.7e-6
temperature = 299.0

[params]
mu_0 = 2.7e-3
N_t = 4.7e16
E_t = -4.4
T_t = 30.0
E_F0 = -4.5
"""


def test_minimal_param_file_loads(tmp_path) -> None:
    path = _write(tmp_path, "p.toml", _MINIMAL)
    cfg = load_config(path)
    assert cfg.params.E_t == -4.4
    assert cfg.name == "X"


def test_param_file_rejects_an_unknown_key(tmp_path) -> None:
    """A misspelled key must fail loudly.

    Ignoring it would leave the model on its default while the file appears to
    say otherwise, and nothing downstream would look wrong.
    """
    path = _write(tmp_path, "p.toml", _MINIMAL.replace("E_t = -4.4", "E_tt = -4.4"))
    with pytest.raises(ValueError, match="unknown key"):
        load_config(path)


def test_param_file_rejects_a_missing_key(tmp_path) -> None:
    path = _write(tmp_path, "p.toml", _MINIMAL.replace("E_t = -4.4\n", ""))
    with pytest.raises(ValueError, match="missing"):
        load_config(path)


def test_param_file_rejects_a_missing_section(tmp_path) -> None:
    path = _write(tmp_path, "p.toml", _MINIMAL.split("[params]")[0])
    with pytest.raises(ValueError, match=r"missing required section \[params\]"):
        load_config(path)


def test_param_file_rejects_an_unknown_modelling_choice(tmp_path) -> None:
    path = _write(
        tmp_path, "p.toml", _MINIMAL + '\n[model]\ntrap_profile = "lorentzian"\n'
    )
    with pytest.raises(ValueError, match="not one of"):
        load_config(path)


def test_param_file_selects_non_default_choices(tmp_path) -> None:
    path = _write(
        tmp_path,
        "p.toml",
        _MINIMAL + '\n[model]\ntrap_profile = "prototype_sech"\n'
        'gamma_model = "Tt/T"\nconstants = "prototype"\n',
    )
    cfg = load_config(path)
    assert cfg.trap_profile is TrapProfile.PROTOTYPE_SECH
    assert cfg.gamma_model is GammaModel.TT_OVER_T
    assert cfg.constants == Constants.prototype()
