import unittest
import numpy as np
from asclc_backend import logarithmic_slope
from asclc_model import calculate_model, K_B, E_CHARGE


class CalculatorChecks(unittest.TestCase):
    def model(self, step=0.003, **overrides):
        params = dict(E_c=-3.36, E_v=-5.58, m_eff_h=0.305, m_eff_e=0.32,
                      eps_r=25.5, thickness=6e-4, temperature=299,
                      mobility=2.7e-3, E_F0=-4.84, E_t=-4.82, N_t=4.7e16,
                      T_t=30, gamma=30/299, trap_profile="logistic",
                      voltage_density="injected")
        params.update(overrides)
        return calculate_model(-np.arange(round(9/step))*step,
                               [-4.84, -4.899, -5.001], energy_step=step, **params)

    def test_equilibrium_and_injection(self):
        result = self.model()
        self.assertEqual(result['delta_p'][0], 0.)
        self.assertEqual(result['V'][0], 0.)
        self.assertEqual(result['J_p'][0], 0.)
        self.assertTrue(np.all(np.diff(result['p_f']) > 0))
        self.assertTrue(np.all(np.diff(result['V']) > 0))

    def test_mobility_changes_current_only(self):
        a, b = self.model(), self.model(mobility=5.4e-3)
        np.testing.assert_allclose(b['J_p'], 2*a['J_p'])
        np.testing.assert_array_equal(b['V'], a['V'])

    def test_quadrature_converges(self):
        coarse, medium, fine = (self.model(step=s) for s in [0.003, 0.0015, 0.00075])
        for key in ['p_f', 'delta_p', 'J_p']:
            error1 = np.abs(coarse[key][1:]-fine[key][1:])
            error2 = np.abs(medium[key][1:]-fine[key][1:])
            self.assertTrue(np.all(error2 < error1), key)
            np.testing.assert_allclose(medium[key][1:], fine[key][1:], rtol=0.01)

    def test_trap_normalization(self):
        for profile in ['logistic', 'gaussian']:
            model = self.model(step=0.00075, trap_profile=profile)
            self.assertAlmostEqual(np.sum(model['g_t'])*0.00075/4.7e16, 1., places=8)

    def test_published_trap_widths(self):
        step = 0.00015
        energy = -np.arange(round(9/step))*step
        ktt = K_B*30/E_CHARGE
        for profile, variance in [('logistic', np.pi**2*ktt**2/3),
                                  ('gaussian', (2*ktt)**2)]:
            with self.subTest(profile=profile):
                gt = self.model(step=step, trap_profile=profile)['g_t']
                mean = np.sum(energy*gt)*step/4.7e16
                self.assertAlmostEqual(mean, -4.82, places=10)
                actual = np.sum((energy+4.82)**2*gt)*step/4.7e16
                self.assertAlmostEqual(actual/variance, 1., places=8)

    def test_unpublished_trap_profile_rejected(self):
        with self.assertRaisesRegex(ValueError, "SI S5.*SI S6"):
            self.model(trap_profile="sech_reference")

    def test_slopes_and_invalid_pairs(self):
        v = np.geomspace(0.01, 100, 25)
        for exponent in [1, 2, 3.7]:
            np.testing.assert_allclose(logarithmic_slope(v, 7*v**exponent)[1], exponent)
        self.assertTrue(np.isnan(logarithmic_slope([1, 1, 2, 3], [1, 2, -1, 4])[1]).all())


if __name__ == '__main__':
    unittest.main()
