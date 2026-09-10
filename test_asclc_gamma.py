import unittest

import numpy as np

from asclc_backend import local_gamma


class LocalGammaChecks(unittest.TestCase):
    def test_power_law_gamma_is_one_over_exponent(self):
        # j = c * V**k  =>  ln V = (ln j - ln c) / k  =>  gamma = d lnV/d lnj = 1/k.
        v = np.geomspace(0.01, 10, 40)
        for k, expected in ((1.0, 1.0), (2.0, 0.5), (0.5, 2.0)):
            j = 3.7 * v ** k
            gamma = local_gamma(v, j, window=5)
            np.testing.assert_allclose(gamma[:-4], expected, rtol=1e-9)
            self.assertTrue(np.isnan(gamma[-4:]).all())

    def test_scale_and_sign_of_current_do_not_change_gamma(self):
        v = np.geomspace(0.1, 5, 25)
        j = 1e-6 * v ** 1.5
        base = local_gamma(v, j, window=7)
        np.testing.assert_allclose(local_gamma(v, -2500 * j, window=7), base,
                                   rtol=1e-9, equal_nan=True)

    def test_centered_places_nan_at_both_ends_and_keeps_the_power_law(self):
        v = np.geomspace(0.01, 10, 40)
        j = 3.7 * v ** 2.0
        gamma = local_gamma(v, j, window=5, centered=True)
        np.testing.assert_allclose(gamma[2:-2], 0.5, rtol=1e-9)
        self.assertTrue(np.isnan(gamma[:2]).all() and np.isnan(gamma[-2:]).all())

    def test_nonfinite_ln_makes_the_spanning_windows_nan(self):
        v = np.geomspace(0.1, 5, 12)
        j = v ** 2.0
        j[4] = 0.0                       # ln|j| -> -inf at index 4
        gamma = local_gamma(v, j, window=3)
        self.assertTrue(np.isnan(gamma[2:5]).all())   # windows that include index 4
        np.testing.assert_allclose(gamma[:2], 0.5, rtol=1e-9)

    def test_shape_and_window_validation(self):
        v = np.geomspace(0.1, 5, 10)
        with self.assertRaises(ValueError):
            local_gamma(v, v[:-1], window=3)
        with self.assertRaises(ValueError):
            local_gamma(v.reshape(2, 5), v.reshape(2, 5), window=3)
        with self.assertRaises(ValueError):
            local_gamma(v, v, window=1)


if __name__ == '__main__':
    unittest.main()
