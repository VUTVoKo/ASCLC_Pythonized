import unittest

import numpy as np

from asclc_backend import (valence_band_dos, conduction_band_dos, total_dos,
                           trap_dos, effective_dos, DOS_FLOOR,
                           K_B, E_CHARGE, H, M_E)

# Shipped MAPbBr3 parameters.
E_V, E_C, M_H, M_E_RATIO = -5.58, -3.36, .305, .32
N_T, E_T, T_T = 4.7e16, -4.82, 30.

# The reference constants round h, e and m0 to 6.62607004e-34, 1.602e-19 and 9.109e-31, which
# scales its C = 4 pi (2 m0 m* e)^1.5/h^3 by a constant against CODATA. The bands
# differ from the cached columns by exactly that factor and nothing else.
C_RATIO = ((M_E/9.109e-31)*(E_CHARGE/1.602e-19))**1.5*(6.62607004e-34/H)**3

# Cached g(E) columns F, I and Q at representative grid energies.
CACHED = ((-6.999, 1.366556424849117e+27, 1e-30, 1.366556424849117e+27),
          (-6.000, 7.43465834866779e+26, 1e-30, 7.43465834866779e+26),
          (-5.583, 6.283433135829589e+25, 1e-30, 6.283433135829589e+25),
          (-5.580, 1e-30, 1e-30, 2e-30),
          (-3.360, 1e-30, 1e-30, 2e-30),
          (-3.357, 1e-30, 6.752618378311664e+25, 6.752618378311664e+25),
          (-3.000, 1e-30, 7.39712281605091e+26, 7.39712281605091e+26),
          (-1.998, 1e-30, 1.438799018391913e+27, 1.438799018391913e+27))
ENERGY, CACHED_F, CACHED_I, CACHED_Q = (np.array(c) for c in zip(*CACHED))


class BandDosChecks(unittest.TestCase):
    def test_reproduces_the_cached_band_columns(self):
        for actual, cached in ((valence_band_dos(ENERGY, E_v=E_V, m_eff_h=M_H), CACHED_F),
                               (conduction_band_dos(ENERGY, E_c=E_C,
                                                    m_eff_e=M_E_RATIO), CACHED_I)):
            band = cached > DOS_FLOOR
            self.assertTrue(band.any())
            np.testing.assert_allclose(actual[band], cached[band], rtol=3e-4)
            ratio = actual[band]/cached[band]        # a pure scale, not a shape error
            np.testing.assert_allclose(ratio, C_RATIO, rtol=1e-12)
            np.testing.assert_array_equal(actual[~band], DOS_FLOOR)

    def test_floor_at_and_beyond_each_band_edge(self):
        self.assertEqual(valence_band_dos(E_V, E_v=E_V, m_eff_h=M_H), DOS_FLOOR)
        self.assertEqual(conduction_band_dos(E_C, E_c=E_C, m_eff_e=M_E_RATIO), DOS_FLOOR)
        self.assertGreater(valence_band_dos(E_V - 1e-9, E_v=E_V, m_eff_h=M_H), DOS_FLOOR)
        self.assertGreater(conduction_band_dos(E_C + 1e-9, E_c=E_C, m_eff_e=M_E_RATIO),
                           DOS_FLOOR)
        np.testing.assert_array_equal(
            valence_band_dos([E_V, -5., -3.], E_v=E_V, m_eff_h=M_H), DOS_FLOOR)

    def test_matches_an_independently_written_parabolic_band(self):
        hbar = H/(2*np.pi)
        for depth in (.001, .05, .42, 1.419):
            # (1/2 pi^2)(2m/hbar^2)^1.5 sqrt(E) in joules, converted to per eV
            expected = (E_CHARGE/(2*np.pi**2)*(2*M_E*M_H/hbar**2)**1.5
                        * np.sqrt(depth*E_CHARGE))
            self.assertAlmostEqual(
                float(valence_band_dos(E_V - depth, E_v=E_V, m_eff_h=M_H))/expected,
                1., 12)
            expected_e = (E_CHARGE/(2*np.pi**2)*(2*M_E*M_E_RATIO/hbar**2)**1.5
                          * np.sqrt(depth*E_CHARGE))
            self.assertAlmostEqual(
                float(conduction_band_dos(E_C + depth, E_c=E_C,
                                          m_eff_e=M_E_RATIO))/expected_e, 1., 12)

    def test_boltzmann_integral_returns_the_effective_dos(self):
        temperature = 299.
        kt = K_B*temperature/E_CHARGE
        depth = np.linspace(0., 80*kt, 800001)         # into the band, eV
        for m_eff, band, edge, sign in ((M_H, valence_band_dos, E_V, -1),
                                        (M_E_RATIO, conduction_band_dos, E_C, +1)):
            kwargs = ({'E_v': edge, 'm_eff_h': m_eff} if sign < 0
                      else {'E_c': edge, 'm_eff_e': m_eff})
            g = band(edge + sign*depth, **kwargs)
            integral = np.trapezoid(g*np.exp(-depth/kt), depth)
            self.assertAlmostEqual(integral/effective_dos(m_eff, temperature), 1., 6)

    def test_mirror_symmetry_of_the_two_bands(self):
        depth = np.array([.002, .1, .8, 2.2])
        np.testing.assert_allclose(
            valence_band_dos(E_V - depth, E_v=E_V, m_eff_h=M_H),
            conduction_band_dos(E_C + depth, E_c=E_C, m_eff_e=M_H), rtol=1e-13)
        self.assertAlmostEqual(
            float(conduction_band_dos(E_C + .5, E_c=E_C, m_eff_e=4*M_H))
            / float(valence_band_dos(E_V - .5, E_v=E_V, m_eff_h=M_H)), 8., 10)


class TotalDosChecks(unittest.TestCase):
    PARAMS = dict(E_v=E_V, E_c=E_C, m_eff_h=M_H, m_eff_e=M_E_RATIO,
                  N_t=N_T, E_t=E_T, T_t=T_T)

    def components(self, energy):
        return (valence_band_dos(energy, E_v=E_V, m_eff_h=M_H),
                conduction_band_dos(energy, E_c=E_C, m_eff_e=M_E_RATIO),
                trap_dos(energy, N_t=N_T, E_t=E_T, T_t=T_T))

    def test_is_the_sum_of_the_three_columns(self):
        energy = np.linspace(-7., -2., 501)
        v, c, t = self.components(energy)
        np.testing.assert_allclose(total_dos(energy, **self.PARAMS), v + c + t,
                                   rtol=1e-15)

    def test_switches_drop_each_term_entirely(self):
        energy = np.array([-6., E_V, E_T, E_C, -3.])
        v, c, t = self.components(energy)
        np.testing.assert_allclose(total_dos(energy, valence=False, **self.PARAMS),
                                   c + t, rtol=1e-15)
        np.testing.assert_allclose(total_dos(energy, conduction=False, **self.PARAMS),
                                   v + t, rtol=1e-15)
        np.testing.assert_allclose(total_dos(energy, trap=False, **self.PARAMS),
                                   v + c, rtol=1e-15)
        np.testing.assert_array_equal(                      # 0*F13, floor included
            total_dos(energy, valence=False, conduction=False, trap=False,
                      **self.PARAMS), np.zeros(5))

    def test_reproduces_the_cached_total_where_a_band_dominates(self):
        dominated = CACHED_Q > 1e20
        np.testing.assert_allclose(total_dos(ENERGY, **self.PARAMS)[dominated],
                                   CACHED_Q[dominated], rtol=3e-4)
        edges = ENERGY == E_V
        np.testing.assert_allclose(total_dos(ENERGY, **self.PARAMS)[edges],
                                   CACHED_Q[edges], rtol=1e-12)

    def test_nonfinite_energy_and_shapes(self):
        out = total_dos([E_T, np.nan, np.inf], **self.PARAMS)
        self.assertTrue(np.isfinite(out[0]) and np.isnan(out[1]) and np.isnan(out[2]))
        for band, kwargs in ((valence_band_dos, {'E_v': E_V, 'm_eff_h': M_H}),
                             (conduction_band_dos, {'E_c': E_C,
                                                    'm_eff_e': M_E_RATIO})):
            values = band([np.nan, -np.inf, np.inf], **kwargs)
            self.assertTrue(np.all(np.isnan(values)))
        self.assertEqual(np.shape(total_dos(np.zeros((2, 3)) - 5., **self.PARAMS)),
                         (2, 3))
        self.assertEqual(np.shape(total_dos(-5., **self.PARAMS)), ())

    def test_invalid_parameters(self):
        for kwargs in (dict(E_v=np.nan, m_eff_h=M_H), dict(E_v=E_V, m_eff_h=0.),
                       dict(E_v=E_V, m_eff_h=-1.), dict(E_v=E_V, m_eff_h=np.inf)):
            with self.assertRaises(ValueError):
                valence_band_dos(-6., **kwargs)
        for kwargs in (dict(E_c=np.inf, m_eff_e=M_E_RATIO), dict(E_c=E_C, m_eff_e=0.),
                       dict(E_c=E_C, m_eff_e=np.nan)):
            with self.assertRaises(ValueError):
                conduction_band_dos(-3., **kwargs)
        with self.assertRaises(ValueError):
            total_dos(-5., **{**self.PARAMS, 'm_eff_e': -1.})


if __name__ == '__main__':
    unittest.main()
