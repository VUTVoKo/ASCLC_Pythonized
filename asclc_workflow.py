"""Notebook-facing loading and export functions; no example file dependencies."""
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from asclc_backend import (Measurements, load_measurements, current_density,
                           trailing_mean, local_gamma, effective_mobility,
                           carrier_densities, model_occupations, model_current_density,
                           valence_band_dos, conduction_band_dos, total_dos,
                           K_B, E_CHARGE, fermi_level_absolute_shift,
                           density_energy_derivative)


@dataclass
class Calculation:
    measurements: Measurements
    V: np.ndarray
    J: np.ndarray


def run_calculation(path, *, device):
    """Load the standard CSV and return measured voltage and current density.

    ``V = U + device['voltage_offset']`` and ``J = I / device['area']``.
    Signs, row order and nonfinite values are preserved; ``device`` is not mutated.
    """
    measurements = load_measurements(path)
    if not np.isfinite(device['voltage_offset']):
        raise ValueError("Voltage offset must be finite.")
    v = measurements.U + device['voltage_offset']
    j = current_density(measurements.I, device['area'])
    return Calculation(measurements, v, j)


@dataclass
class Mobility:
    U: np.ndarray
    j: np.ndarray
    gamma: np.ndarray
    mu_eff: np.ndarray
    window: int


def run_effective_mobility(result, *, material, window):
    """Steps A2 and A3 on a common trailing window of raw measurements.

    U and j are trailing means; gamma fits ln|U| on ln|j| using raw
    readings over the same window. Equation (4) gives mu_eff from these
    quantities, thickness and permittivity. Incomplete windows return nan;
    material is not mutated and row order and length are preserved."""
    u = trailing_mean(result.V, window=window)
    j = trailing_mean(result.J, window=window)
    gamma = local_gamma(result.V, result.J, window=window, centered=False)
    mu_eff = effective_mobility(u, j, gamma, thickness=material['L'],
                                eps_r=material['eps_r'])
    return Mobility(u, j, gamma, mu_eff, int(window))


@dataclass
class Carriers:
    U: np.ndarray
    j: np.ndarray
    gamma: np.ndarray
    p_f: np.ndarray
    p_t: np.ndarray
    p_s: np.ndarray
    theta: np.ndarray
    mu_0: float
    theta_model: str
    window: int


def run_carrier_densities(mobility, *, material, mu_0,
                          theta_model="free_over_total"):
    """Guide step A4: equations (5) and (6) on the rows of step A3.

    ``mobility`` is the ``Mobility`` returned by ``run_effective_mobility``, so
    each row keeps the same trailing window, voltage and slope the mobility was
    built from: ``p_f[i]``, ``p_t[i]`` and ``theta[i]`` describe ``U[i]``,
    ``j[i]`` and ``gamma[i]``, and ``theta`` equals ``mu_eff[i] / mu_0`` under
    the default reading of equation (6).

    ``mu_0`` is the microscopic mobility read off the step-A3 plot (m^2 V^-1
    s^-1); it scales ``p_f`` and leaves equation (6) untouched. ``theta_model``
    selects which reading of equation (6) is used -- see
    ``asclc_backend.carrier_densities``. ``material`` is not mutated; arrays
    keep the length and order of ``mobility.U``.
    """
    p_f, p_t, p_s, theta = carrier_densities(
        mobility.U, mobility.j, mobility.gamma, thickness=material['L'],
        eps_r=material['eps_r'], mu_0=mu_0, theta_model=theta_model)
    return Carriers(mobility.U, mobility.j, mobility.gamma, p_f, p_t, p_s, theta,
                    float(mu_0), theta_model, mobility.window)


def run_model_occupations(energy, E_F, *, material, params, temperature):
    """M2 populations from M1 DOS and the supplied equilibrium reference.

    Grids and parameters are supplied explicitly and are not mutated.
    Uses the existing backend physical constants."""
    bands = dict(E_v=material['E_v'], E_c=material['E_c'],
                 m_eff_h=material['m_eff_p'], m_eff_e=material['m_eff_e'])
    traps = {key: params[key] for key in ('N_t', 'E_t', 'T_t')}
    return model_occupations(
        energy, E_F, g_total=total_dos(energy, **bands, **traps),
        g_conduction=conduction_band_dos(
            energy, E_c=bands['E_c'], m_eff_e=bands['m_eff_e']),
        g_valence=valence_band_dos(
            energy, E_v=bands['E_v'], m_eff_h=bands['m_eff_h']),
        E_F0=params['E_F0'], thermal_energy=K_B * temperature / E_CHARGE)


@dataclass
class ModelFermiLevel:
    """M3 absolute quasi-Fermi level and signed energy differences, all in eV."""

    E_F: np.ndarray
    E_F0: float
    delta_E_v: np.ndarray
    delta_E_c: np.ndarray
    delta_E_F0: np.ndarray


@dataclass
class MeasuredFermiLevel:
    U: np.ndarray
    E_F: np.ndarray
    N_v: np.ndarray
    nondegenerate: np.ndarray


def run_measured_fermi_level(carriers, *, material, temperature):
    """Absolute-shift inversion on the same intervals as measured holes."""
    ef, nv, valid = fermi_level_absolute_shift(
        carriers.p_f, E_v=material['E_v'], m_eff_h=material['m_eff_p'],
        temperature=temperature)
    return MeasuredFermiLevel(carriers.U.copy(), ef, nv, valid)


@dataclass
class DensityDerivatives:
    measured: np.ndarray
    model: np.ndarray


def run_density_derivatives(carriers, measured_levels, model):
    """Derivatives of measured pf+pt and model excess total holes."""
    if not np.array_equal(carriers.U, measured_levels.U, equal_nan=True):
        raise ValueError('Measured energies and carriers must share intervals.')
    return DensityDerivatives(
        density_energy_derivative(measured_levels.E_F, carriers.p_f + carriers.p_t),
        density_energy_derivative(model.E_F, model.p_s, model=True))


def run_model_fermi_level(model, *, material):
    """M3 energy coordinates of the existing M2 occupation sweep.

    Occupations already associate each E_F with carrier densities at the
    supplied M2 temperature. Preserve that coordinate and row order; no
    Boltzmann inversion or resampling is needed. The hole separation is
    E_F-E_v (equation 7); E_F-E_F0 is a distinct equilibrium-referenced
    change. Signed differences remain valid inside and outside the bands.
    All energies are in eV on the same reference scale.
    """
    ef = np.asarray(model.E_F, dtype=float)
    ev, ec = material['E_v'], material['E_c']
    if ef.ndim != 1 or not np.isfinite(ef).all():
        raise ValueError("E_F must be a finite 1-D sweep.")
    if not np.isfinite([ev, ec, model.E_F0]).all() or ec <= ev:
        raise ValueError("Band edges and E_F0 must be finite, with E_c > E_v.")
    return ModelFermiLevel(ef.copy(), float(model.E_F0),
                           ef - ev, ef - ec, ef - model.E_F0)


@dataclass
class ModelCurrent:
    """M4 separate electron/hole U (V), J (A/m^2) on the M2 energy sweep."""

    E_F: np.ndarray
    U_n: np.ndarray
    J_n: np.ndarray
    U_p: np.ndarray
    J_p: np.ndarray
    gamma: float


def run_model_current(model, *, material, mu_0, gamma):
    """Map M2 labelled n_t/p_t to charge density, n_f/p_f to drift density.

    Retains the established q = s - f0 convention. These labelled trapped
    populations are not generally the total charge required by Poisson's
    equation. The resulting curves are an algebraic continuation of M4;
    they need not reproduce M5's mu_0 * abs(f/s). Both carriers are separate
    single-carrier constructions with the same supplied gamma and mobility.
    """
    ef = np.asarray(model.E_F, dtype=float)
    if ef.ndim != 1 or not np.isfinite(ef).all():
        raise ValueError("E_F must be a finite 1-D sweep.")
    if any(np.shape(a) != ef.shape for a in
           (model.n_t, model.n_f, model.p_t, model.p_f)):
        raise ValueError("Carrier densities must match the E_F sweep.")
    kwargs = dict(gamma=gamma, thickness=material['L'],
                  eps_r=material['eps_r'], mu_0=mu_0)
    un, jn = model_current_density(model.n_t, model.n_f, **kwargs)
    up, jp = model_current_density(model.p_t, model.p_f, **kwargs)
    return ModelCurrent(ef.copy(), un, jn, up, jp, float(gamma))


@dataclass
class ModelMobility:
    """M5 carrier-resolved mobilities in m^2 V^-1 s^-1 on the M2 sweep."""

    E_F: np.ndarray
    mu_n: np.ndarray
    mu_p: np.ndarray
    mu_0: float


def run_model_mobility(model, *, mu_0):
    """M5: mu_eff = mu_0 * theta for each M2 carrier population.

    The supplied microscopic mobility applies separately to both carriers;
    this is not a combined bipolar mobility. Preserve M2's absolute-free /
    excess-total ratios, including values above one and singularities.
    No slope estimation, fitting, or resampling is performed.
    """
    if not np.isscalar(mu_0) or not np.isfinite(mu_0) or mu_0 <= 0:
        raise ValueError("mu_0 must be finite and positive (m^2 V^-1 s^-1).")
    ef = np.asarray(model.E_F, dtype=float)
    qn, qp = (np.asarray(q, dtype=float) for q in (model.theta_n, model.theta_p))
    if ef.ndim != 1 or not np.isfinite(ef).all():
        raise ValueError("E_F must be a finite 1-D sweep.")
    if qn.shape != ef.shape or qp.shape != ef.shape:
        raise ValueError("Carrier fractions must match the E_F sweep.")
    return ModelMobility(ef.copy(), mu_0 * qn, mu_0 * qp, float(mu_0))


def export_results(result, figures, directory):
    """Save each named figure as PNG/SVG plus the measured CSV.

    ``figures`` maps a base filename (no extension) to a Matplotlib figure.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    for name, figure in figures.items():
        for extension in ('png', 'svg'):
            figure.savefig(directory / f'{name}.{extension}', dpi=240,
                           bbox_inches='tight')
    np.savetxt(directory / 'measured.csv',
               np.column_stack((result.measurements.U, result.measurements.I,
                                result.V, result.J)),
               delimiter=',', comments='',
               header='U_V,I_A,V_V,J_Am2')
    return directory
