"""Numerical DOS/occupation model. Inputs are arrays and physical parameters only."""
import numpy as np
from scipy.constants import (
    elementary_charge as E_CHARGE,
    Boltzmann as K_B,
    epsilon_0 as EPS0,
    Planck as H,
    electron_mass as M_E,
)


def _occupation(x):
    return np.exp(-np.logaddexp(0.0, -x))


def calculate_model(energy, fermi_levels, *, energy_step, E_c, E_v,
                    m_eff_h, m_eff_e, eps_r, thickness, temperature,
                    mobility, E_F0, E_t, N_t, T_t, gamma,
                    trap_profile, voltage_density):
    """Integrate DOS with rectangular weights and calculate both drift currents.

    trap_profile: 'logistic' (SI Eq. S5, called biexponential there) or
    'gaussian' (SI Eq. S6, sigma = 2*k_B*T_t in energy units).
    Both profiles integrate to N_t over the full energy axis.
    voltage_density: 'reference' subtracts equilibrium free carriers from
    injected total density; 'injected' uses the injected total directly.
    These conventions are explicit because they produce different curves.
    gamma is a prescribed transport parameter, not the measured inverse slope.
    """
    e = np.asarray(energy, dtype=float)
    ef = np.asarray(fermi_levels, dtype=float)
    positive = [energy_step, m_eff_h, m_eff_e, eps_r, thickness,
                temperature, mobility, T_t]
    if not all(np.isfinite(x) and x > 0 for x in positive):
        raise ValueError("Grid step and physical scale parameters must be positive and finite.")
    if not (0 < gamma < 1) or not np.isfinite(N_t) or N_t < 0:
        raise ValueError("Require 0 < gamma < 1 and finite N_t >= 0.")
    if not all(np.isfinite(x) for x in [E_c, E_v, E_F0, E_t]) or E_v >= E_c:
        raise ValueError("Require finite energy parameters and E_v < E_c.")
    if (e.ndim != 1 or e.size < 2 or ef.ndim != 1 or ef.size < 1
            or not np.isfinite(e).all() or not np.isfinite(ef).all()
            or not (np.all(np.diff(e) > 0) or np.all(np.diff(e) < 0))
            or not np.allclose(np.abs(np.diff(e)), energy_step)):
        raise ValueError("Supply a uniform energy grid and a finite 1D Fermi sweep.")
    kt, ktt = K_B * temperature / E_CHARGE, K_B * T_t / E_CHARGE
    cv = 4*np.pi*(2*M_E*m_eff_h*E_CHARGE)**1.5/H**3
    cc = 4*np.pi*(2*M_E*m_eff_e*E_CHARGE)**1.5/H**3
    gv = cv*np.sqrt(np.maximum(E_v-e, 0))
    gc = cc*np.sqrt(np.maximum(e-E_c, 0))
    if trap_profile == "logistic":
        z = np.exp(-np.abs((e-E_t)/ktt))
        gt = N_t/ktt*z/(1+z)**2
    elif trap_profile == "gaussian":
        sigma = 2*ktt
        gt = N_t/(sigma*np.sqrt(2*np.pi))*np.exp(-0.5*((e-E_t)/sigma)**2)
    else:
        raise ValueError("trap_profile must be 'logistic' (SI S5) or 'gaussian' (SI S6).")
    if voltage_density not in ("reference", "injected"):
        raise ValueError("Unknown voltage_density.")
    total = gv + gc + gt
    h0, n0 = _occupation((e-E_F0)/kt), _occupation((E_F0-e)/kt)
    pf0, nf0 = np.sum(gv*h0)*energy_step, np.sum(gc*n0)*energy_step
    pf, nf, delta = (np.empty(ef.size) for _ in range(3))
    for i, f in enumerate(ef):
        holes, electrons = _occupation((e-f)/kt), _occupation((f-e)/kt)
        # Use small complementary occupations when holes are nearly saturated.
        difference = np.where(e > max(f, E_F0), n0-electrons, holes-h0)
        delta[i] = np.sum(total*difference)*energy_step
        pf[i], nf[i] = np.sum(gv*holes)*energy_step, np.sum(gc*electrons)*energy_step
    qp = delta - pf0 if voltage_density == "reference" else delta
    qn = -delta - nf0 if voltage_density == "reference" else -delta
    factor = E_CHARGE*thickness**2/(EPS0*eps_r*(1-gamma)*(2-gamma))
    vp, vn = factor*qp, factor*qn
    drift = E_CHARGE*mobility*(2-gamma)/thickness
    return dict(E_F=ef, V=vp, V_n=vn, J_p=drift*vp*pf, J_n=drift*vn*nf,
                p_f=pf, n_f=nf, delta_p=delta, p_f0=pf0, n_f0=nf0,
                g_v=gv, g_c=gc, g_t=gt, mobility=mobility)


def slope_guides(voltage, *, m_eff_h, temperature, E_F0, E_v,
                 mobility, thickness, eps_r, ohmic_scale=1., quadratic_scale=1.):
    """Ohmic and Mott-Gurney laws with explicit display scale factors."""
    v = np.asarray(voltage, dtype=float)
    nv = 2*(2*np.pi*M_E*m_eff_h*K_B*temperature/H**2)**1.5
    pf0 = nv*np.exp(-(E_F0-E_v)*E_CHARGE/(K_B*temperature))
    return (ohmic_scale*E_CHARGE*mobility*pf0*v/thickness,
            quadratic_scale*9/8*EPS0*eps_r*mobility*v**2/thickness**3)
