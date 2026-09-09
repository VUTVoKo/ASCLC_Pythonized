"""Notebook-facing calculation and export functions; no example file dependencies."""
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from asclc_backend import (Measurements, load_measurements, current_density,
                           effective_mobility, carrier_densities)
from asclc_model import calculate_model, slope_guides


@dataclass
class Calculation:
    measurements: Measurements
    V: np.ndarray
    J: np.ndarray
    T_K: np.ndarray
    V_mid: np.ndarray
    m: np.ndarray
    model: dict
    series: list
    J_mid: np.ndarray
    gamma: np.ndarray
    mu_eff: np.ndarray
    p_f: np.ndarray
    p_t: np.ndarray
    p_s: np.ndarray
    theta: np.ndarray
    analysis_mobility: float


def _grid(start, stop, step, *, endpoint):
    if not np.isfinite([start, stop, step]).all() or step <= 0 or start == stop:
        raise ValueError("Grid endpoints must differ and step must be positive and finite.")
    count = abs(stop-start)/step
    if not np.isclose(count, round(count), rtol=0, atol=1e-8):
        raise ValueError("Grid span must be an integer multiple of its step.")
    return start + np.sign(stop-start)*step*np.arange(round(count) + int(endpoint))


def run_calculation(path, *, material, device, model, numerics, analysis_mobility=None):
    """Load standard CSV and calculate measured slopes, model and plot guides.

    Model temperature is the configured device temperature in K. Measured
    temperatures are returned separately in T_K and do not change that policy.
    Configuration dictionaries are never mutated. Energy stop is exclusive;
    Fermi stop is inclusive. Both grids use the configured energy step.
    A4 treats Eq. (6) as trapped density and uses Theta = p_f/(p_f+p_t).
    analysis_mobility defaults to model['mobility']; an explicit value overrides
    it for measured carrier densities without changing the model.
    """
    measurements = load_measurements(path)
    if not np.isfinite(device['voltage_offset']):
        raise ValueError("Voltage offset must be finite.")
    v = measurements.U + device['voltage_offset']
    j = current_density(measurements.I, device['area'])
    v_mid, j_mid, m, gamma, mu_eff = effective_mobility(
        v, j, thickness=device['thickness'], eps_r=material['eps_r'])
    if analysis_mobility is None:
        analysis_mobility = model['mobility']
    pf, pt, ps, theta = carrier_densities(
        v_mid, j_mid, gamma, thickness=device['thickness'],
        eps_r=material['eps_r'], mobility=analysis_mobility)
    step = numerics['energy_step']
    energy = _grid(*numerics['energy_range'], step, endpoint=False)
    fermi = _grid(*numerics['fermi_range'], step, endpoint=True)
    modeled = calculate_model(
        energy, fermi, energy_step=step, **material,
        thickness=device['thickness'], temperature=device['temperature'], **model,
    )
    guide_v = np.geomspace(1e-3, 100, 201)
    guides = slope_guides(
        guide_v, m_eff_h=material['m_eff_h'], temperature=device['temperature'],
        E_F0=model['E_F0'], E_v=material['E_v'], mobility=model['mobility'],
        thickness=device['thickness'], eps_r=material['eps_r'],
        ohmic_scale=numerics['ohmic_scale'], quadratic_scale=numerics['quadratic_scale'],
    )
    series = [(v, j), (modeled['V'], modeled['J_p']),
              (modeled['V'], -modeled['J_n']), (guide_v, guides[0]), (guide_v, guides[1])]
    return Calculation(measurements, v, j, measurements.T+273.15,
                       v_mid, m, modeled, series, j_mid, gamma, mu_eff,
                       pf, pt, ps, theta, analysis_mobility)


def export_results(result, figures, directory):
    """Save each named figure as PNG/PDF/SVG plus measured and calculated CSVs.

    ``figures`` maps a base filename (no extension) to a Matplotlib figure.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    for name, figure in figures.items():
        for extension in ('png', 'pdf', 'svg'):
            figure.savefig(directory / f'{name}.{extension}', dpi=240,
                           bbox_inches='tight')
    np.savetxt(directory / 'measured_carriers.csv',
               np.column_stack((result.V_mid, result.J_mid, result.gamma,
                                result.p_f, result.p_t, result.p_s, result.theta,
                                np.full(result.V_mid.shape, result.analysis_mobility))),
               delimiter=',', comments='',
               header='V_mid_V,J_mid_Am2,gamma,p_f_m3,p_t_eq6_m3,p_s_pf_plus_pt_m3,'
                      'theta_pf_over_ps,mu0_m2_per_Vs')
    np.savetxt(directory / 'measured_mobility.csv',
               np.column_stack((result.V_mid, result.J_mid, result.m,
                                result.gamma, result.mu_eff)), delimiter=',',
               header='V_mid_V,J_mid_Am2,m_dlnJ_dlnV,gamma,mu_eff_m2_per_Vs', comments='')
    np.savetxt(directory / 'measured_slopes.csv',
               np.column_stack((result.V_mid, result.m)), delimiter=',',
               header='V_mid_V,m_dlnJ_dlnV', comments='')
    keys = ('E_F', 'V', 'J_p', 'J_n', 'p_f', 'n_f', 'delta_p')
    np.savetxt(directory / 'model_calculated.csv',
               np.column_stack([result.model[key] for key in keys]), delimiter=',',
               header='E_F_eV,V_V,J_p_Am2,J_n_Am2,p_f_m3,n_f_m3,delta_p_m3', comments='')
    return directory
