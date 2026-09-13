"""M2: reference evidence and independent occupation identities, fixed inputs."""
from pathlib import Path
import unittest

import numpy as np

from asclc_backend import model_occupations


class ModelOccupationsChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with np.load(Path(__file__).parent / 'tests/data/m2_reference.npz') as data:
            cls.energy = data['energy'].copy()
            cls.cached = dict(zip(data['reference_keys'], data['reference_values']))
            cls.kwargs = {key: data[key].copy() for key in
                          ('g_total', 'g_conduction', 'g_valence')}
            cls.kwargs.update(E_F0=float(data['E_F0']),
                              thermal_energy=float(data['thermal_energy']))

    def test_cached_reference_populations_and_theta(self):
        rows = [9, 1009, 1396, 1409, 1609, 2009, 3009]
        ef = np.array([self.cached[f'A{i}'] for i in rows])
        model = model_occupations(self.energy, ef, **self.kwargs)
        for field, column in [('n_f', 'E'), ('p_f', 'L')]:
            np.testing.assert_allclose(getattr(model, field),
                [self.cached[f'{column}{i}'] for i in rows], rtol=2e-11)
        # Subtracting ~1e27 backgrounds loses ~1e12 per rounding unit.
        # Independent summation methods accumulate in different orders. Bound that absolute
        # error instead of demanding relative agreement near equilibrium.
        background_error = 32 * np.spacing(max(model.n_s0, model.p_s0))
        for field, column in [('n_s', 'C'), ('p_s', 'J'), ('n_t', 'H'), ('p_t', 'O')]:
            np.testing.assert_allclose(getattr(model, field),
                [self.cached[f'{column}{i}'] for i in rows],
                rtol=2e-11, atol=background_error)
        for field, col, total in [('theta_n', 'G', model.n_s), ('theta_p', 'N', model.p_s)]:
            expected = np.array([self.cached[f'{col}{i}'] for i in rows])
            relative_bound = 2e-11 + 2 * background_error / np.abs(total)
            self.assertTrue(np.all(np.abs(getattr(model, field) / expected - 1) < relative_bound))

    def test_equilibrium_singularity_and_fixed_free_reference(self):
        model = model_occupations(self.energy, [self.kwargs['E_F0']], **self.kwargs)
        self.assertEqual(model.n_s[0], 0)
        self.assertEqual(model.p_s[0], 0)
        self.assertEqual(model.n_t[0], -model.n_f0)
        self.assertEqual(model.p_t[0], -model.p_f0)
        self.assertTrue(np.isinf(model.theta_n[0]))
        self.assertTrue(np.isinf(model.theta_p[0]))

    def test_charge_complementarity_and_monotonicity(self):
        ef = np.array([self.cached[f'A{i}'] for i in (9, 1009, 1409, 2009, 3009)])
        model = model_occupations(self.energy, ef, **self.kwargs)
        error = 16 * np.spacing(max(model.n_s0, model.p_s0))
        np.testing.assert_allclose(model.n_s + model.p_s, 0, atol=error)
        self.assertTrue(np.all(np.diff(model.n_f) >= 0))
        self.assertTrue(np.all(np.diff(model.p_f) <= 0))
        # Full occupied + empty state count follows f + (1-f) = 1.
        available = np.sum(-np.diff(self.energy) * self.kwargs['g_total'][:-1])
        self.assertAlmostEqual((model.n_s0 + model.p_s0) / available, 1, 14)

    def test_grid_validation_and_inputs_preserved(self):
        original = self.energy.copy()
        with self.assertRaises(ValueError):
            model_occupations(self.energy[::-1], [self.kwargs['E_F0']], **self.kwargs)
        model_occupations(self.energy, [], **self.kwargs)
        np.testing.assert_array_equal(self.energy, original)


if __name__ == '__main__':
    unittest.main()
