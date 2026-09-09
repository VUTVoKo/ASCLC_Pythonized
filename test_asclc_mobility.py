import unittest

import numpy as np

from asclc_backend import effective_mobility
from asclc_model import EPS0


class MobilityChecks(unittest.TestCase):
    def test_mott_gurney_recovers_mobility_in_both_sweep_directions(self):
        mu, length, eps_r = 2.7e-3, 6e-4, 25.5
        for v in (np.geomspace(.01, 100, 40), np.geomspace(100, .01, 40)):
            j = 9/8 * EPS0 * eps_r * mu * v**2 / length**3
            vm, jm, m, gamma, actual = effective_mobility(
                v, j, thickness=length, eps_r=eps_r)
            np.testing.assert_allclose(actual, mu, rtol=1e-12)
            np.testing.assert_allclose(gamma, .5)
            np.testing.assert_allclose(jm, 9/8 * EPS0 * eps_r * mu * vm**2 / length**3)

    def test_general_power_law_and_geometry_scaling(self):
        v = np.geomspace(.1, 10, 15)
        vm, jm, m, g, mu = effective_mobility(v, 3*v**4, thickness=.001, eps_r=10)
        np.testing.assert_allclose(mu, .001**3*3*vm**2/(EPS0*10*.75*1.75**2))
        other = effective_mobility(v, 3*v**4, thickness=.002, eps_r=20)[-1]
        np.testing.assert_allclose(other, 4*mu)

    def test_ohmic_sublinear_and_invalid_intervals(self):
        for exponent in [1, .5, 0, -1]:
            v = np.geomspace(.1, 10, 20)
            self.assertTrue(np.isnan(effective_mobility(
                v, v**exponent, thickness=1, eps_r=1)[-1]).all())
        result = effective_mobility([1, 2, 2, 3, 4, 5, 6],
                                    [1, 4, 4, -1, np.nan, 25, 36],
                                    thickness=1, eps_r=1)
        np.testing.assert_array_equal(np.isfinite(result[-1]), [True, False, False, False, False, True])

    def test_invalid_geometry(self):
        for value in (0, -1, np.nan, np.inf):
            for key in ('thickness', 'eps_r'):
                kwargs = dict(thickness=1, eps_r=1)
                kwargs[key] = value
                with self.assertRaises(ValueError):
                    effective_mobility([1, 2], [1, 4], **kwargs)


if __name__ == '__main__':
    unittest.main()
