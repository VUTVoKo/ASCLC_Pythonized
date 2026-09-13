import unittest

import numpy as np

from asclc_backend import EPSILON_0, effective_mobility


def mott_gurney_current(voltage, *, mu, thickness, eps_r):
    """J = (9/8) eps_0 eps_r mu U^2 / L^3, the gamma = 1/2 limit of Eq. (4)."""
    return 9 / 8 * EPSILON_0 * eps_r * mu * np.asarray(voltage, float) ** 2 / thickness ** 3


class EffectiveMobilityChecks(unittest.TestCase):
    material = dict(thickness=6.0e-4, eps_r=25.5)

    def test_mott_gurney_current_returns_its_own_mobility(self):
        # Independent analytic case: at gamma = 1/2 Eq. (4) must invert the
        # Mott-Gurney law exactly, whatever the mobility.
        v = np.array([0.5, 1.0, 2.0, 3.0])
        for mu in (2.7e-3, 1.0, 4.1e-7):
            j = mott_gurney_current(v, mu=mu, **self.material)
            mu_eff = effective_mobility(v, j, np.full(v.shape, 0.5), **self.material)
            np.testing.assert_allclose(mu_eff, mu, rtol=1e-12)

    def test_shape_factor_scales_the_mott_gurney_value(self):
        # Away from gamma = 1/2 the result is the Mott-Gurney value times
        # (9/8) / ((1 - g)(2 - g)^2).
        v = np.array([1.0, 2.0, 3.0])
        j = mott_gurney_current(v, mu=2.7e-3, **self.material)
        for g in (0.0, 0.25, 0.75, 1.5, -0.4):
            expected = 2.7e-3 * (9 / 8) / ((1 - g) * (2 - g) ** 2)
            mu_eff = effective_mobility(v, j, np.full(v.shape, g), **self.material)
            np.testing.assert_allclose(mu_eff, expected, rtol=1e-12)

    def test_sub_ohmic_slope_gives_a_negative_mobility(self):
        v, j = np.array([2.0]), np.array([1e-3])
        self.assertLess(effective_mobility(v, j, np.array([1.5]), **self.material)[0], 0.0)

    def test_dimensional_scaling_of_thickness_and_permittivity(self):
        v, j, g = np.array([2.0]), np.array([1e-3]), np.array([0.5])
        base = effective_mobility(v, j, g, **self.material)[0]
        doubled_l = effective_mobility(v, j, g, thickness=2 * self.material['thickness'],
                                       eps_r=self.material['eps_r'])[0]
        doubled_eps = effective_mobility(v, j, g, thickness=self.material['thickness'],
                                         eps_r=2 * self.material['eps_r'])[0]
        self.assertAlmostEqual(doubled_l / base, 8.0, places=9)      # L^3
        self.assertAlmostEqual(doubled_eps / base, 0.5, places=9)    # 1 / eps_r

    def test_units_are_m2_per_volt_second(self):
        # 1 A/m^2 at 1 V across 1 m of vacuum, gamma = 1/2: mu = (8/9)/eps_0.
        mu_eff = effective_mobility(np.array([1.0]), np.array([1.0]), np.array([0.5]),
                                    thickness=1.0, eps_r=1.0)[0]
        self.assertAlmostEqual(mu_eff, (8 / 9) / EPSILON_0, delta=1e-3)

    def test_current_sign_does_not_change_the_mobility(self):
        v = np.array([-3.0, -1.0, 1.0, 3.0])
        j = np.array([-4e-3, -1e-3, 1e-3, 4e-3])
        g = np.full(v.shape, 0.5)
        np.testing.assert_allclose(effective_mobility(v, j, g, **self.material),
                                   effective_mobility(v, -j, g, **self.material),
                                   rtol=1e-12)
        np.testing.assert_allclose(effective_mobility(v, j, g, **self.material)[:2][::-1],
                                   effective_mobility(v, j, g, **self.material)[2:],
                                   rtol=1e-12)

    def test_undefined_points_are_nan(self):
        v = np.array([1.0, 0.0, 1.0, 1.0, 1.0, np.nan])
        j = np.array([1e-3, 1e-3, 1e-3, 1e-3, np.nan, 1e-3])
        g = np.array([1.0, 0.5, 2.0, np.nan, 0.5, 0.5])
        self.assertTrue(np.isnan(effective_mobility(v, j, g, **self.material)).all())

    def test_rejects_mismatched_arrays_and_unphysical_geometry(self):
        v = np.array([1.0, 2.0, 3.0])
        with self.assertRaises(ValueError):
            effective_mobility(v, v[:-1], v, **self.material)
        with self.assertRaises(ValueError):
            effective_mobility(v.reshape(1, 3), v.reshape(1, 3), v.reshape(1, 3), **self.material)
        with self.assertRaises(ValueError):
            effective_mobility(v, v, v, thickness=0.0, eps_r=25.5)
        with self.assertRaises(ValueError):
            effective_mobility(v, v, v, thickness=6e-4, eps_r=-1.0)


if __name__ == "__main__":
    unittest.main()
