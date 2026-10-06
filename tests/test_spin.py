import unittest

import numpy as np

from stateid.spin import build_s2_ph, determinant_s2, lr_s2, tda_s2


class SpinTests(unittest.TestCase):
    def test_closed_shell(self):
        c = np.eye(3)[:, :2]
        self.assertAlmostEqual(determinant_s2(c, c), 0)

    def test_doublet_and_triplet(self):
        self.assertAlmostEqual(determinant_s2(np.eye(2)[:, :1], np.empty((2, 0))), .75)
        self.assertAlmostEqual(determinant_s2(np.eye(2), np.empty((2, 0))), 2)

    def test_beta_majority(self):
        self.assertAlmostEqual(determinant_s2(np.empty((2, 0)), np.eye(2)), 2)

    def test_broken_symmetry_determinant(self):
        self.assertAlmostEqual(determinant_s2([[1], [0]], [[0], [1]]), 1)
        # Nonorthogonal alpha-beta spatial orbitals: 1-|overlap|² = 3/4.
        self.assertAlmostEqual(determinant_s2([[1], [0]], [[.5], [np.sqrt(.75)]]), .75)

    def test_nonorthogonal_ao_and_complex_orbitals(self):
        a = np.array([[1, .4j], [0, 1.2]], dtype=complex)
        s = a.conj().T @ a
        ca = np.linalg.solve(a, [[1], [0]])
        cb = np.linalg.solve(a, [[0], [1j]])
        self.assertAlmostEqual(determinant_s2(ca, cb, s), 1)

    def test_requires_normalization(self):
        with self.assertRaises(ValueError):
            determinant_s2([[2]], [[1]])

    def test_tda_absolute_s2(self):
        op = np.diag([0, 2])
        self.assertAlmostEqual(tda_s2([1, 1j], op), 1)
        self.assertAlmostEqual(tda_s2([0, 3j], op), 2)

    def test_tda_invalid_inputs(self):
        for x, op in [([0, 0], np.eye(2)), ([1, 0], [[0, 1], [0, 2]]),
                      ([1, np.nan], np.eye(2)), ([[1, 0]], np.eye(2)), ([1], np.eye(2))]:
            with self.subTest(x=x), self.assertRaises(ValueError):
                tda_s2(x, op)

    def test_future_apis_fail_explicitly(self):
        with self.assertRaises(NotImplementedError):
            build_s2_ph(None, None)
        with self.assertRaises(NotImplementedError):
            lr_s2(None, None, None, None)


if __name__ == "__main__":
    unittest.main()
