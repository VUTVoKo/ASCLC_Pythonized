"""Independent analytic cases for density-versus-energy differentiation."""
import unittest

import numpy as np

from asclc_backend import density_energy_derivative


class DensityDerivativeChecks(unittest.TestCase):
    def test_constant_and_affine_populations(self):
        energy = np.array([-5.1, -5., -4.8, -4.5])
        for model in (False, True):
            count = 2 if model else 1
            flat = density_energy_derivative(energy, np.full(4, 8e16), model=model)
            np.testing.assert_allclose(flat[:-count], 0)
            slope = density_energy_derivative(energy, -3e16 * energy + 2e16,
                                              model=model)
            np.testing.assert_allclose(slope[:-count], 3e16, rtol=1e-13)
            self.assertTrue(np.isnan(slope[-count:]).all())

    def test_quadratic_center_slope_is_assigned_to_first_row(self):
        energy = np.array([-5.2, -5.1, -5., -4.9])
        actual = density_energy_derivative(energy, 1e16 * energy**2, model=True)
        np.testing.assert_allclose(actual[:2], -2e16 * energy[1:3], rtol=1e-12)

    def test_sign_and_invalid_windows(self):
        energy = np.array([-5., -4.9, -4.8, -4.7])
        density = np.array([1., 2., 3., 4.]) * 1e16
        self.assertTrue((density_energy_derivative(energy, density)[:-1] > 0).all())
        self.assertTrue((density_energy_derivative(energy, density, model=True)[:-2] < 0).all())
        density[1] = np.nan
        self.assertTrue(np.isnan(density_energy_derivative(energy, density, model=True)).all())
        actual = density_energy_derivative(energy, density)
        self.assertTrue(np.isnan(actual[:2]).all())
        self.assertTrue(np.isfinite(actual[2]))
        self.assertTrue(np.isnan(density_energy_derivative([1., 1.], [2., 3.])).all())


if __name__ == '__main__':
    unittest.main()
