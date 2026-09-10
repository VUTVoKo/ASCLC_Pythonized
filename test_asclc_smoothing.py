import unittest

import numpy as np

from asclc_backend import centered_mean, trailing_mean


class TrailingMeanChecks(unittest.TestCase):
    def test_matches_explicit_trailing_windows(self):
        x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        np.testing.assert_allclose(trailing_mean(x, window=3),
                                   [2.0, 3.0, 4.0, np.nan, np.nan])

    def test_window_one_is_identity(self):
        x = np.array([3.0, -1.0, 4.0, 1.5])
        np.testing.assert_array_equal(trailing_mean(x, window=1), x)

    def test_window_covering_a_nonfinite_value_is_nan(self):
        x = np.array([1.0, np.nan, 3.0, 4.0, 5.0, 6.0])
        out = trailing_mean(x, window=2)
        self.assertTrue(np.isnan(out[:2]).all())        # windows touching index 1
        np.testing.assert_allclose(out[2:4], [3.5, 4.5])

    def test_shape_and_window_validation(self):
        with self.assertRaises(ValueError):
            trailing_mean(np.zeros((2, 3)), window=2)
        with self.assertRaises(ValueError):
            trailing_mean(np.zeros(5), window=0)


class CenteredMeanChecks(unittest.TestCase):
    def test_odd_window_is_symmetric_and_nan_at_both_ends(self):
        x = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
        np.testing.assert_allclose(centered_mean(x, window=3),
                                   [np.nan, 2.0, 3.0, 4.0, 5.0, 6.0, np.nan])

    def test_even_window_takes_the_extra_point_from_the_leading_side(self):
        x = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
        # window=4: two points lead, one trails, so index 0-1 and 5 have no window
        np.testing.assert_allclose(centered_mean(x, window=4),
                                   [np.nan, np.nan, 1.5, 2.5, 3.5, np.nan])

    def test_stays_aligned_with_a_linear_series(self):
        x = np.linspace(0.0, 10.0, 21)
        out = centered_mean(x, window=5)
        np.testing.assert_allclose(out[2:-2], x[2:-2])   # mean of a line == centre

    def test_window_one_is_identity_and_validation(self):
        x = np.array([3.0, -1.0, 4.0])
        np.testing.assert_array_equal(centered_mean(x, window=1), x)
        with self.assertRaises(ValueError):
            centered_mean(np.zeros((2, 3)), window=3)


if __name__ == '__main__':
    unittest.main()
