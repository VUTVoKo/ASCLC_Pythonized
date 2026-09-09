import unittest

import numpy as np

from asclc_backend import carrier_densities, effective_mobility
from asclc_model import EPS0, E_CHARGE


class CarrierChecks(unittest.TestCase):
    def test_mott_gurney_and_density_convention(self):
        mu, length, eps_r = 2.7e-3, 6e-4, 25.5
        for v in (np.geomspace(.01, 100, 40), np.geomspace(100, .01, 40)):
            j = 9/8 * EPS0 * eps_r * mu * v**2 / length**3
            vm, jm, _, g, mu_eff = effective_mobility(
                v, j, thickness=length, eps_r=eps_r)
            pf, pt, ps, theta = carrier_densities(
                vm, jm, g, thickness=length, eps_r=eps_r, mobility=mu)
            expected = .75 * EPS0 * eps_r * vm / (E_CHARGE * length**2)
            np.testing.assert_allclose(pf, expected, rtol=1e-12)
            np.testing.assert_allclose(pt, expected, rtol=1e-12)
            np.testing.assert_allclose(ps, 2*expected, rtol=1e-12)
            np.testing.assert_allclose(theta, .5, rtol=1e-12)
            np.testing.assert_allclose(pf/pt, mu_eff/mu, rtol=1e-12)

    def test_mobility_changes_free_density_only(self):
        kwargs = dict(thickness=.001, eps_r=10)
        a = carrier_densities([1, 2], [3, 12], [.5, .25], mobility=.002, **kwargs)
        b = carrier_densities([1, 2], [3, 12], [.5, .25], mobility=.004, **kwargs)
        np.testing.assert_allclose(b[0], a[0]/2)
        np.testing.assert_allclose(b[1], a[1])
        np.testing.assert_allclose(b[3], b[0]/(b[0]+b[1]))

    def test_invalid_intervals_retained(self):
        result = carrier_densities(
            [1, 2, 3, 4, 0, 6, 7, 8, 9],
            [1, 4, 9, 16, 25, -1, np.inf, 64, 81],
            [.5, 1, 0, 1-5e-13, .5, .5, .5, np.nan, .25],
            thickness=1, eps_r=1, mobility=1)
        for values in result:
            np.testing.assert_array_equal(np.isfinite(values),
                                          [True, False, False, False, False,
                                           False, False, False, True])

    def test_invalid_parameters_and_shapes(self):
        for key in ('thickness', 'eps_r', 'mobility'):
            for value in (0, -1, np.inf, np.nan):
                kwargs = dict(thickness=1, eps_r=1, mobility=1)
                kwargs[key] = value
                with self.assertRaises(ValueError):
                    carrier_densities([1], [1], [.5], **kwargs)
        for arrays in (([1], [1, 2], [.5]), (1, 1, .5)):
            with self.assertRaises(ValueError):
                carrier_densities(*arrays, thickness=1, eps_r=1, mobility=1)


if __name__ == '__main__':
    unittest.main()
