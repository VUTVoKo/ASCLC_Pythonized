import math
import unittest

import numpy as np

from asclc_backend import trap_dos, K_B, E_CHARGE

# Shipped MAPbBr3 S2 trap parameters.
N_T, E_T, T_T = 4.7e16, -4.82, 30.

# The reference data rounds k_B and e to 1.38e-23 and 1.602e-19, which makes its
# k_B T_t 0.036 % smaller than CODATA's. Passing the T_t that reproduces the
# reference k_B T_t exactly isolates the profile's form from that difference.
T_T_REFERENCE = (1.38e-23*T_T/1.602e-19)*E_CHARGE/K_B

# Cached values of the reference trap profile, at representative energies.
CACHED = ((-4.821, 4.2263538478369044e+18),   # nearest grid point to the peak
          (-4.818, 3.458334134557805e+18),
          (-4.800, 3959754928233774.5),
          (-4.770, 35981223464.07873),
          (-4.701, 0.09128790905390642),
          (-4.500, 1.5196506214396976e-35),
          (-5.199, 1.8476168167088505e-45))   # the far side of E_t


class TrapDosChecks(unittest.TestCase):
    def test_reproduces_the_cached_reference_column(self):
        energy, cached = (np.array(c) for c in zip(*CACHED))
        actual = trap_dos(energy, N_t=N_T, E_t=E_T, T_t=T_T_REFERENCE)
        np.testing.assert_allclose(actual, cached, rtol=1e-15)

    def test_matches_an_independently_written_sech(self):
        kt_t = K_B*T_T/E_CHARGE
        for energy in (E_T, E_T + .01, E_T - .01, E_T + .2, E_T - .25):
            expected = N_T/(4*kt_t*math.cosh((energy - E_T)/kt_t))
            self.assertAlmostEqual(
                float(trap_dos(energy, N_t=N_T, E_t=E_T, T_t=T_T))/expected, 1., 12)

    def test_integrates_to_pi_over_four_N_t(self):
        energy = np.linspace(E_T - 1.2, E_T + 1.2, 2400001)
        integral = np.trapezoid(trap_dos(energy, N_t=N_T, E_t=E_T, T_t=T_T), energy)
        self.assertAlmostEqual(integral/N_T, np.pi/4, 9)
        self.assertAlmostEqual(integral, 3.6914e16, delta=1e12)

    def test_peak_height_and_symmetry_about_E_t(self):
        kt_t = K_B*T_T/E_CHARGE
        peak = trap_dos(E_T, N_t=N_T, E_t=E_T, T_t=T_T)
        self.assertAlmostEqual(float(peak)/(N_T/(4*kt_t)), 1., 12)
        offsets = np.array([1e-4, 2.5e-3, .01, .1, .7])
        np.testing.assert_allclose(trap_dos(E_T + offsets, N_t=N_T, E_t=E_T, T_t=T_T),
                                   trap_dos(E_T - offsets, N_t=N_T, E_t=E_T, T_t=T_T),
                                   rtol=0, atol=0)
        self.assertTrue(np.all(trap_dos(E_T + offsets, N_t=N_T, E_t=E_T, T_t=T_T) < peak))

    def test_far_tail_underflows_instead_of_overflowing(self):
        kt_t = K_B*T_T/E_CHARGE
        u = np.array([745., 1000., 2000.])           # the reference grid reaches 745
        with np.errstate(over="ignore"):
            self.assertTrue(np.all(np.isinf(np.cosh(u))))   # the naive form dies here
        tail = trap_dos(E_T + u*kt_t, N_t=N_T, E_t=E_T, T_t=T_T)
        self.assertTrue(np.all(np.isfinite(tail)))
        self.assertTrue(np.all(np.diff(tail) <= 0))
        self.assertEqual(tail[-1], 0.)

    def test_scales_with_N_t_and_inversely_with_T_t(self):
        energy = np.linspace(E_T - .02, E_T + .02, 9)
        base = trap_dos(energy, N_t=N_T, E_t=E_T, T_t=T_T)
        np.testing.assert_allclose(trap_dos(energy, N_t=3*N_T, E_t=E_T, T_t=T_T),
                                   3*base, rtol=1e-13)
        self.assertAlmostEqual(
            float(trap_dos(E_T, N_t=N_T, E_t=E_T, T_t=2*T_T))
            / float(trap_dos(E_T, N_t=N_T, E_t=E_T, T_t=T_T)), .5, 12)
        shifted = trap_dos(energy + .3, N_t=N_T, E_t=E_T + .3, T_t=T_T)
        np.testing.assert_allclose(shifted, base, rtol=1e-13)

    def test_nonfinite_energy_and_shape(self):
        out = trap_dos([E_T, np.nan, np.inf, -np.inf], N_t=N_T, E_t=E_T, T_t=T_T)
        self.assertTrue(np.isfinite(out[0]))
        self.assertTrue(np.isnan(out[1]))
        np.testing.assert_array_equal(out[2:], [0., 0.])
        self.assertEqual(np.shape(trap_dos(np.zeros((2, 3)) + E_T,
                                           N_t=N_T, E_t=E_T, T_t=T_T)), (2, 3))
        self.assertEqual(trap_dos(E_T, N_t=0., E_t=E_T, T_t=T_T), 0.)

    def test_invalid_parameters(self):
        for kwargs in (dict(N_t=np.nan, E_t=E_T, T_t=T_T),
                       dict(N_t=-1., E_t=E_T, T_t=T_T),
                       dict(N_t=N_T, E_t=np.inf, T_t=T_T),
                       dict(N_t=N_T, E_t=E_T, T_t=0.),
                       dict(N_t=N_T, E_t=E_T, T_t=-30.),
                       dict(N_t=N_T, E_t=E_T, T_t=np.nan)):
            with self.assertRaises(ValueError):
                trap_dos(E_T, **kwargs)


if __name__ == '__main__':
    unittest.main()
