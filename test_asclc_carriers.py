import unittest

import numpy as np

from asclc_backend import (EPSILON_0, E_CHARGE, THETA_MODELS, carrier_densities,
                           effective_mobility, equal_density_voltage)

L, EPS_R, MU_0 = 6.0e-4, 25.5, 2.7e-3


def mott_gurney(voltage, *, mobility=MU_0, thickness=L, eps_r=EPS_R):
    """Equation (3), written out independently of the module under test."""
    return 9 / 8 * EPSILON_0 * eps_r * mobility * voltage ** 2 / thickness ** 3


def ohms_law(voltage, p_f, gamma, *, mobility=MU_0, thickness=L):
    """Equation (S14): j = e mu_0 p_f (2 - gamma) U / L."""
    return E_CHARGE * mobility * p_f * (2 - gamma) * voltage / thickness


class CarrierDensityChecks(unittest.TestCase):
    def test_eq5_inverts_ohms_law_for_any_slope(self):
        """Equation (5) returns the p_f that equation (S14) was built from."""
        v = np.geomspace(0.01, 3.0, 24)
        expected = np.geomspace(1e12, 1e18, 24)
        for gamma in (0.0, 0.25, 0.5, 0.9, 1.5, -0.2):
            g = np.full(v.shape, gamma)
            p_f, _, _, _ = carrier_densities(
                v, ohms_law(v, expected, gamma), g,
                thickness=L, eps_r=EPS_R, mu_0=MU_0)
            np.testing.assert_allclose(p_f, expected, rtol=1e-12)

    def test_eq6_matches_an_independent_evaluation_and_units(self):
        v, g = np.array([0.5, 1.0, 2.5]), np.array([0.5, 0.0, 0.25])
        p_f, p_t, p_s, _ = carrier_densities(
            v, mott_gurney(v), g, thickness=L, eps_r=EPS_R, mu_0=MU_0,
            theta_model="absolute_over_total")
        expected = EPSILON_0 * EPS_R * (1 - g) * (2 - g) * v / (E_CHARGE * L ** 2)
        np.testing.assert_allclose(p_t, expected, rtol=1e-12)
        np.testing.assert_allclose(p_s, np.abs(p_f) + np.abs(p_t), rtol=1e-12)
        # Hand-checked magnitude, m^-3: eps_0 eps_r (0.5)(1.5)(0.5) / (e L^2).
        self.assertAlmostEqual(p_t[0] / 1.4685e15, 1.0, places=3)

    def test_mott_gurney_limit_has_no_trapped_charge(self):
        """Exact Eq. (3) data at mu_0: p_f = p_s, p_t = 0, Theta = 1."""
        v = np.geomspace(0.01, 100, 40)
        p_f, p_t, p_s, theta = carrier_densities(
            v, mott_gurney(v), np.full(v.shape, 0.5),
            thickness=L, eps_r=EPS_R, mu_0=MU_0)
        expected = 0.75 * EPSILON_0 * EPS_R * v / (E_CHARGE * L ** 2)
        np.testing.assert_allclose(p_f, expected, rtol=1e-12)
        np.testing.assert_allclose(p_s, expected, rtol=1e-12)
        np.testing.assert_allclose(p_t, 0.0, atol=1e-4 * expected[0])
        np.testing.assert_allclose(theta, 1.0, rtol=1e-12)

    def test_theta_is_the_mobility_ratio_of_equation_4(self):
        """Theta = mu_eff / mu_0 identically, on scattered data."""
        rng = np.random.default_rng(20260913)
        v = np.geomspace(0.05, 3.0, 120)
        j = mott_gurney(v) * rng.lognormal(0.0, 0.4, v.size)
        g = rng.uniform(-0.2, 0.9, v.size)
        mu_eff = effective_mobility(v, j, g, thickness=L, eps_r=EPS_R)
        p_f, _, p_s, theta = carrier_densities(
            v, j, g, thickness=L, eps_r=EPS_R, mu_0=MU_0)
        np.testing.assert_allclose(theta, mu_eff / MU_0, rtol=1e-12)
        np.testing.assert_allclose(theta, p_f / p_s, rtol=1e-12)

    def test_equal_densities_mean_half_the_microscopic_mobility(self):
        """The article's test of mu_0: p_t = p_f at mu_eff = mu_0 / 2."""
        v = np.geomspace(0.05, 3.0, 60)
        j = mott_gurney(v, mobility=MU_0 / 2)
        g = np.full(v.shape, 0.5)
        mu_eff = effective_mobility(v, j, g, thickness=L, eps_r=EPS_R)
        p_f, p_t, _, theta = carrier_densities(
            v, j, g, thickness=L, eps_r=EPS_R, mu_0=MU_0)
        np.testing.assert_allclose(mu_eff, MU_0 / 2, rtol=1e-12)
        np.testing.assert_allclose(p_t, p_f, rtol=1e-12)
        np.testing.assert_allclose(theta, 0.5, rtol=1e-12)

    def test_the_three_readings_of_equation_6(self):
        v, j, g = np.array([2.0]), np.array([1e-2]), np.array([0.25])
        kwargs = dict(thickness=L, eps_r=EPS_R, mu_0=MU_0)
        eq6 = EPSILON_0 * EPS_R * 0.75 * 1.75 * v / (E_CHARGE * L ** 2)
        free = carrier_densities(v, j, g, **kwargs)[0]
        results = {name: carrier_densities(v, j, g, theta_model=name, **kwargs)
                   for name in THETA_MODELS}
        for name, (p_f, _, _, _) in results.items():
            np.testing.assert_allclose(p_f, free, rtol=1e-12, err_msg=name)
        np.testing.assert_allclose(results["free_over_total"][1], eq6 - free)
        np.testing.assert_allclose(results["free_over_total"][3], free / eq6)
        for name in ("absolute_over_total", "mobility_ratio"):
            np.testing.assert_allclose(results[name][1], eq6, err_msg=name)
        np.testing.assert_allclose(results["absolute_over_total"][2], free + eq6)
        np.testing.assert_allclose(results["absolute_over_total"][3],
                                   free / (free + eq6))
        np.testing.assert_allclose(results["mobility_ratio"][2], free + eq6)
        np.testing.assert_allclose(results["mobility_ratio"][3], free / eq6)

    def test_absolute_reading_survives_a_sign_flipping_noise_floor(self):
        v = np.array([0.1, 0.2, 0.3])
        j = np.array([5e-9, -4e-9, 6e-9])
        g = np.array([-0.1, 0.05, 0.0])
        for name in THETA_MODELS:
            p_f, p_t, p_s, theta = carrier_densities(
                v, j, g, thickness=L, eps_r=EPS_R, mu_0=MU_0, theta_model=name)
            for values in (p_f, p_t, p_s, theta):
                self.assertTrue(np.isfinite(values).all(), name)
            self.assertTrue((p_f > 0).all(), name)

    def test_scaling_with_mu_0_thickness_and_permittivity(self):
        v, j, g = np.array([1.5]), np.array([3e-3]), np.array([0.4])
        base = carrier_densities(v, j, g, thickness=L, eps_r=EPS_R, mu_0=MU_0,
                                 theta_model="mobility_ratio")
        doubled = carrier_densities(v, j, g, thickness=L, eps_r=EPS_R,
                                    mu_0=2 * MU_0, theta_model="mobility_ratio")
        np.testing.assert_allclose(doubled[0], base[0] / 2, rtol=1e-12)
        np.testing.assert_allclose(doubled[1], base[1], rtol=1e-12)
        np.testing.assert_allclose(doubled[3], base[3] / 2, rtol=1e-12)
        thicker = carrier_densities(v, j, g, thickness=2 * L, eps_r=EPS_R,
                                    mu_0=MU_0, theta_model="mobility_ratio")
        np.testing.assert_allclose(thicker[0], 2 * base[0], rtol=1e-12)
        np.testing.assert_allclose(thicker[1], base[1] / 4, rtol=1e-12)
        denser = carrier_densities(v, j, g, thickness=L, eps_r=2 * EPS_R,
                                   mu_0=MU_0, theta_model="mobility_ratio")
        np.testing.assert_allclose(denser[0], base[0], rtol=1e-12)
        np.testing.assert_allclose(denser[1], 2 * base[1], rtol=1e-12)

    def test_reverse_branch_maps_onto_the_forward_one(self):
        v = np.geomspace(0.05, 3.0, 30)
        j, g = mott_gurney(v), np.full(v.shape, 0.4)
        forward = carrier_densities(v, j, g, thickness=L, eps_r=EPS_R, mu_0=MU_0)
        reverse = carrier_densities(-v, -j, g, thickness=L, eps_r=EPS_R, mu_0=MU_0)
        for a, b in zip(forward, reverse):
            np.testing.assert_allclose(a, b, rtol=1e-12)

    def test_shape_factor_sign_and_singular_points(self):
        v = np.array([1.0, 1.0, 1.0, 1.0])
        j = np.array([1e-3, 1e-3, 1e-3, 1e-3])
        g = np.array([1.5, 2.0, 1.0, 0.5])
        p_f, p_t, _, _ = carrier_densities(
            v, j, g, thickness=L, eps_r=EPS_R, mu_0=MU_0,
            theta_model="absolute_over_total")
        self.assertLess(p_t[0], 0.0)            # 1 < gamma < 2: Eq. (6) negative
        self.assertTrue(np.isnan(p_f[1]))       # gamma = 2: Eq. (5) singular
        self.assertAlmostEqual(p_t[2], 0.0)     # gamma = 1: no space charge
        self.assertGreater(p_t[3], 0.0)
        theta = carrier_densities(v, j, g, thickness=L, eps_r=EPS_R, mu_0=MU_0)[3]
        self.assertTrue(np.isnan(theta[1]))     # gamma = 2: the whole row goes
        self.assertTrue(np.isnan(theta[2]))     # gamma = 1: theta is 0/0

    def test_nonfinite_inputs_and_zero_voltage(self):
        v = np.array([1.0, 0.0, 2.0, np.nan, 3.0])
        j = np.array([1e-3, 1e-3, np.inf, 1e-3, 1e-3])
        g = np.array([0.5, 0.5, 0.5, 0.5, np.nan])
        for name in THETA_MODELS:
            values = carrier_densities(v, j, g, thickness=L, eps_r=EPS_R,
                                       mu_0=MU_0, theta_model=name)
            for array in values:
                np.testing.assert_array_equal(
                    np.isfinite(array), [True, False, False, False, False],
                    err_msg=name)

    def test_invalid_parameters_and_shapes(self):
        v, j, g = np.array([1.0]), np.array([1e-3]), np.array([0.5])
        for key in ("thickness", "eps_r", "mu_0"):
            for value in (0.0, -1.0, np.inf, np.nan):
                kwargs = dict(thickness=L, eps_r=EPS_R, mu_0=MU_0)
                kwargs[key] = value
                with self.assertRaises(ValueError):
                    carrier_densities(v, j, g, **kwargs)
        with self.assertRaises(ValueError):
            carrier_densities(v, j, g, thickness=L, eps_r=EPS_R, mu_0=MU_0,
                              theta_model="paper")
        for arrays in ((v, np.array([1e-3, 2e-3]), g), (1.0, 1e-3, 0.5),
                       (np.zeros((2, 2)), np.zeros((2, 2)), np.zeros((2, 2)))):
            with self.assertRaises(ValueError):
                carrier_densities(*arrays, thickness=L, eps_r=EPS_R, mu_0=MU_0)


class EqualDensityChecks(unittest.TestCase):
    def test_crossing_is_where_mu_eff_is_half_mu_0(self):
        """A scan through the crossing: Theta = 1/2 at the returned voltage."""
        v = np.geomspace(0.1, 10.0, 200)
        g = np.full(v.shape, 0.5)
        # j rising as V^3 crosses the Mott-Gurney V^2 shape at a known voltage.
        j = mott_gurney(v, mobility=MU_0 / 2) * (v / 4.0)
        p_f, p_t, _, theta = carrier_densities(
            v, j, g, thickness=L, eps_r=EPS_R, mu_0=MU_0)
        crossing = equal_density_voltage(v, p_f, p_t)
        self.assertEqual(crossing.size, 1)
        self.assertAlmostEqual(crossing[0] / 4.0, 1.0, places=3)
        np.testing.assert_allclose(np.interp(4.0, v, theta), 0.5, rtol=1e-3)

    def test_no_crossing_and_exact_rows(self):
        v = np.array([1.0, 2.0, 3.0])
        self.assertEqual(equal_density_voltage(v, [1e14, 2e14, 3e14],
                                               [1e16, 2e16, 3e16]).size, 0)
        np.testing.assert_allclose(
            equal_density_voltage(v, [1e14, 2e14, 3e14], [1e14, 5e14, 9e14]),
            [1.0])

    def test_nonpositive_rows_are_skipped(self):
        v = np.array([1.0, 2.0, 3.0, 4.0])
        p_f = np.array([1e14, np.nan, 1e16, 1e17])
        p_t = np.array([1e16, 1e15, -1e14, 1e16])
        crossing = equal_density_voltage(v, p_f, p_t)
        self.assertEqual(crossing.size, 1)
        self.assertTrue(1.0 < crossing[0] < 4.0)
        with self.assertRaises(ValueError):
            equal_density_voltage(v, p_f[:2], p_t)


if __name__ == "__main__":
    unittest.main()
