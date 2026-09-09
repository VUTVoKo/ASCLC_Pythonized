from pathlib import Path
from tempfile import NamedTemporaryFile
import unittest

import numpy as np

from asclc_backend import load_measurements


class StandardCSVChecks(unittest.TestCase):
    def load(self, text):
        with NamedTemporaryFile(dir=Path(__file__).resolve().parent,
                                suffix='.csv', delete=False) as stream:
            path = Path(stream.name)
        try:
            path.write_text(text, encoding='utf-8-sig')
            return load_measurements(path)
        finally:
            path.unlink()

    def test_standard_units_and_preservation(self):
        result = self.load('U(V),I(A),T(C)\n2,-3e-9,25\n1,nan,20\n')
        np.testing.assert_array_equal(result.U, [2, 1])
        self.assertEqual(result.I[0], -3e-9)
        self.assertTrue(np.isnan(result.I[1]))
        np.testing.assert_array_equal(result.T, [25, 20])
        np.testing.assert_array_equal(result.finite, [True, False])
        self.assertEqual(result.units, dict(U='V', I='A', T='C'))

    def test_wrong_units_order_and_column_count_rejected(self):
        for contents in ['U(V),I(A),T(K)\n1,2,300\n',
                         'I(A),U(V),T(C)\n1,2,25\n',
                         'U(V),I(A),T(C)\n1,2\n',
                         'U(V),I(A),T(C)\n1,2,25,4\n']:
            with self.subTest(contents=contents), self.assertRaises(ValueError):
                self.load(contents)


if __name__ == '__main__':
    unittest.main()
