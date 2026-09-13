import unittest
from pathlib import Path

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

    def test_scatter_pins_the_regression_direction(self):
        # gamma is least-squares fit(ln|U| over the window, ln|j| over the window): ln U
        # regressed on ln j. On scattered data that is a different number from
        # 1/(slope of ln j on ln U), so this pins which way the fit runs.
        rng = np.random.default_rng(0)
        v = np.logspace(-1, 0.5, 60)
        j = 3.0 * v ** 2.0 * np.exp(rng.normal(0.0, 0.35, v.size))
        gamma = local_gamma(v, j, window=11, centered=True)
        x, y = np.log(v), np.log(j)
        for i in range(5, v.size - 5):
            xs, ys = x[i-5:i+6], y[i-5:i+6]
            xc, yc = xs - xs.mean(), ys - ys.mean()
            regression = np.dot(yc, xc) / np.dot(yc, yc)          # ln U on ln j
            inverted = np.dot(xc, xc) / np.dot(xc, yc)        # 1 / (ln j on ln U)
            self.assertAlmostEqual(gamma[i], regression, places=12)
            self.assertNotAlmostEqual(gamma[i], inverted, places=2)

    def test_flat_current_has_no_gamma(self):
        # m = 0 makes gamma = 1/m unbounded; report nan rather than an infinity.
        v = np.logspace(-1, 1, 20)
        gamma = local_gamma(v, np.full(v.size, 4e-9), window=5, centered=True)
        self.assertTrue(np.isnan(gamma).all())

    def test_ohmic_and_mott_gurney_reference_slopes(self):
        v = np.logspace(-2, 0.5, 40)
        for exponent, expected in ((1.0, 1.0), (2.0, 0.5)):
            gamma = local_gamma(v, 1e-6 * v ** exponent, window=7, centered=True)
            finite = gamma[np.isfinite(gamma)]
            np.testing.assert_allclose(finite, expected, rtol=1e-9)

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
