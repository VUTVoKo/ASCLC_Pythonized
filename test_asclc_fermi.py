import unittest

import numpy as np

from asclc_backend import fermi_level, effective_dos, K_B, E_CHARGE, H, M_E


class FermiLevelChecks(unittest.TestCase):
    def test_eq7_round_trip_scalar_and_per_point_temperature(self):
        e_v, m_eff_h = -5.58, .305
        for temperature in (299., np.linspace(285., 310., 12)):
            expected = np.linspace(-5.2, -4.6, 12)
            kt = K_B*np.asarray(temperature, float)/E_CHARGE
            n_v = effective_dos(m_eff_h, temperature)
            p_f = n_v*np.exp(-(expected-e_v)/kt)
            actual, returned_n_v, nondegenerate = fermi_level(
                p_f, E_v=e_v, m_eff_h=m_eff_h, temperature=temperature)
            np.testing.assert_allclose(actual, expected, rtol=1e-12)
            np.testing.assert_allclose(returned_n_v, np.broadcast_to(n_v, p_f.shape))
            self.assertTrue(nondegenerate.all())

    def test_per_point_temperature_differs_from_scalar(self):
        p_f = np.full(3, 1.4e12)
        t = np.array([290., 299., 310.])
        per_point = fermi_level(p_f, E_v=-5.58, m_eff_h=.305, temperature=t)
        scalar = fermi_level(p_f, E_v=-5.58, m_eff_h=.305, temperature=299.)
        self.assertAlmostEqual(per_point[0][1], scalar[0][1], places=12)
        self.assertTrue(np.all(np.diff(per_point[0]) > 0))
        self.assertGreater(abs(per_point[0][0]-scalar[0][0]), 1e-3)
        for i, temperature in enumerate(t):
            single = fermi_level(p_f[i:i+1], E_v=-5.58, m_eff_h=.305,
                                 temperature=temperature)
            self.assertAlmostEqual(per_point[0][i], single[0][0], places=12)
            self.assertAlmostEqual(per_point[1][i], single[1][0], places=6)

    def test_effective_dos_matches_its_definition(self):
        for m_eff, temperature in ((.305, 299.), (.32, 300.)):
            expected = 2*(2*np.pi*M_E*m_eff*K_B*temperature/H**2)**1.5
            self.assertAlmostEqual(effective_dos(m_eff, temperature)/expected, 1., 12)
        np.testing.assert_allclose(effective_dos(.305, np.array([299., 598.])),
                                   effective_dos(.305, 299.)*np.array([1., 2**1.5]))

    def test_degeneracy_flag(self):
        e_v, m_eff_h, temperature = -5.58, .305, 299.
        n_v = effective_dos(m_eff_h, temperature)
        kt = K_B*temperature/E_CHARGE
        # The flag turns over at E_F-E_v = 3kT; the boundary itself is a
        # convention, so check either side of it rather than rounding on it.
        p_f = n_v*np.exp(-np.array([3.1, 3.05, 2.9]))
        e_f, _, nondegenerate = fermi_level(p_f, E_v=e_v, m_eff_h=m_eff_h,
                                            temperature=temperature)
        np.testing.assert_allclose(e_f, e_v + kt*np.array([3.1, 3.05, 2.9]))
        np.testing.assert_array_equal(nondegenerate, [True, True, False])

    def test_invalid_points_retained_as_nan(self):
        p_f = np.array([1e12, 0., -1e12, np.nan, np.inf, 1e12])
        t = np.array([299., 299., 299., 299., 299., np.nan])
        e_f, n_v, nondegenerate = fermi_level(p_f, E_v=-5.58, m_eff_h=.305,
                                              temperature=t)
        for values in (e_f, n_v):
            np.testing.assert_array_equal(
                np.isfinite(values), [True, False, False, False, False, False])
        np.testing.assert_array_equal(
            nondegenerate, [True, False, False, False, False, False])

    def test_invalid_parameters_and_shapes(self):
        for kwargs in (dict(E_v=np.nan, m_eff_h=.305), dict(E_v=-5.58, m_eff_h=0),
                       dict(E_v=-5.58, m_eff_h=-1), dict(E_v=-5.58, m_eff_h=np.inf)):
            with self.assertRaises(ValueError):
                fermi_level([1e12], temperature=299., **kwargs)
        for arguments in (([1e12, 2e12], [299.]), ([[1e12]], 299.),
                          ([1e12], [[299., 300.]])):
            with self.assertRaises(ValueError):
                fermi_level(arguments[0], E_v=-5.58, m_eff_h=.305,
                            temperature=arguments[1])


if __name__ == '__main__':
    unittest.main()
