"""MODEL!K and MODEL!N reference amplitudes on the log/log current chart."""
import unittest
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import numpy as np

from asclc_backend import E_CHARGE, EPSILON_0
from asclc_plotting import (plot_model_current, _OHMIC_FACTOR,
                            _MOTT_GURNEY_FACTOR)
from asclc_workflow import (ModelCurrent, run_calculation,
                            run_model_occupations, run_model_current)

ROOT = Path(__file__).resolve().parent
MATERIAL = dict(L=6.00e-4, eps_r=25.50)
MU_0 = 2.7e-3
P_F0 = 1399565233278.4204   # MODEL!J2, the workbook's own nf0


def reference_lines(axes):
    """The m = 1 and m = 2 lines, by label, as (voltage, current) arrays."""
    lines = {line.get_label(): line for line in axes.get_lines()}
    return ({name: lines[name].get_xydata().T for name in ('$m=1$', '$m=2$')})


class WorkbookReferenceAmplitudes(unittest.TestCase):
    def setUp(self):
        ef = np.linspace(-4.84, -5.30, 64)
        one = np.ones_like(ef)
        self.model = ModelCurrent(ef, one, one, np.geomspace(1e-2, 4.0, 64),
                                  np.geomspace(1e-7, 1e-2, 64), 0.5)

    def figure(self):
        return plot_model_current(self.model, material=MATERIAL, mu_0=MU_0,
                                  equilibrium_holes=P_F0)

    def test_reference_endpoints_match_the_workbook_cells(self):
        # MODEL!K and MODEL!N evaluated at MODEL!I2 = 1e-3 V and I3 = 100 V.
        # Both already carry L3 = 2 and O3 = 0.4; the constants differ from
        # the workbook's rounded e and eps_0 by under 2e-4 relative.
        figure, axes = self.figure()
        try:
            lines = reference_lines(axes)
        finally:
            matplotlib.pyplot.close(figure)
        for name, cells in (('$m=1$', (2.0178931533408267e-9,
                                       2.0178931533408267e-4)),
                            ('$m=2$', (7.9374726562500008e-9,
                                       79.374726562500015))):
            u, j = lines[name]
            slope = np.diff(np.log(j)) / np.diff(np.log(u))
            exponent = 1.0 if name == '$m=1$' else 2.0
            np.testing.assert_allclose(slope, exponent, rtol=1e-9)
            scale = j[0] / u[0] ** exponent
            np.testing.assert_allclose(
                [scale * 1e-3 ** exponent, scale * 100.0 ** exponent],
                cells, rtol=2e-4)

    def test_factors_are_the_only_departure_from_the_bare_equations(self):
        figure, axes = self.figure()
        try:
            lines = reference_lines(axes)
        finally:
            matplotlib.pyplot.close(figure)
        u = lines['$m=1$'][0]
        bare_ohmic = E_CHARGE * MU_0 * P_F0 * u / MATERIAL['L']
        bare_mg = (9 / 8 * EPSILON_0 * MATERIAL['eps_r'] * MU_0 * u ** 2
                   / MATERIAL['L'] ** 3)
        np.testing.assert_allclose(lines['$m=1$'][1],
                                   bare_ohmic * _OHMIC_FACTOR, rtol=1e-12)
        np.testing.assert_allclose(lines['$m=2$'][1],
                                   bare_mg / _MOTT_GURNEY_FACTOR, rtol=1e-12)

    def test_the_m4_curve_stays_between_both_references(self):
        # J/J_MG = (8/9)(1-g)(2-g)**2 * theta and J/J_ohmic = (2-g) * pf/pf0.
        # Neither ratio is bounded by the bare equations at gamma = T_t/T; the
        # workbook factors are what place the hole branch inside the band over
        # the plotted range. Uses the notebook's parameters and sweep.
        material = dict(L=6.00e-4, S=7.70e-6, eps_r=25.50, E_c=-3.36,
                        E_v=-5.58, m_eff_p=0.305, m_eff_e=0.320)
        params = dict(N_t=4.7e16, E_t=-4.82, T_t=30.0, E_F0=-4.84)
        temperature = 299.0
        measured = run_calculation(
            ROOT / 'data/MAPbBr3_S2_dark.csv',
            device=dict(area=material['S'], voltage_offset=0.0))
        occupations = run_model_occupations(
            -np.arange(3001) * 0.003, -9.0 + np.arange(3001) * 0.003,
            material=material, params=params, temperature=temperature)
        current = run_model_current(occupations, material=material, mu_0=MU_0,
                                    gamma=params['T_t'] / temperature)
        figure, axes = plot_model_current(
            current, measured=measured, material=material, mu_0=MU_0,
            equilibrium_holes=occupations.p_f0)
        try:
            lines = reference_lines(axes)
            low, high = axes.get_xlim()
        finally:
            matplotlib.pyplot.close(figure)
        u, j = current.U_p, current.J_p
        shown = (np.isfinite(u) & np.isfinite(j) & (j > 0)
                 & (u >= low) & (u <= high))
        self.assertGreater(shown.sum(), 50)
        u, j = u[shown], j[shown]
        for name, exponent, compare in (('$m=1$', 1.0, np.greater_equal),
                                        ('$m=2$', 2.0, np.less_equal)):
            ur, jr = lines[name]
            scale = jr[0] / ur[0] ** exponent
            self.assertTrue(compare(j, scale * u ** exponent).all(), name)
        # Without the factors the curve leaves the band on the m = 2 side.
        bare = (9 / 8 * EPSILON_0 * material['eps_r'] * MU_0 * u ** 2
                / material['L'] ** 3)
        self.assertTrue((j > bare).any())


if __name__ == '__main__':
    unittest.main()
