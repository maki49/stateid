import unittest

import numpy as np
from numpy.testing import assert_allclose

from stateid.symmetry import (
    analyze_c3v, class_characters, complex_to_real_basis, magnetic_order,
    match_characters, orthonormalize, particle_hole_operation,
    project_representation, rotation_in_real_basis,
)


def c3v_ops():
    # Analytic Cartesian A1(z)+E(x,y); test includes all six elements.
    c, s = -0.5, np.sqrt(3) / 2
    r = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    f = np.diag([1, 1, -1])
    return {"E": np.eye(3), "C3": r, "C3^2": r @ r,
            "sv1": f, "sv2": r @ f, "sv3": r @ r @ f}


class SymmetryTests(unittest.TestCase):
    def test_a1_and_e_from_coefficients(self):
        for cols, chars, label in [([0], [1, 1, 1], "A1"), ([1, 2], [2, -1, 0], "E")]:
            with self.subTest(label=label):
                result = analyze_c3v(np.eye(3)[:, cols], c3v_ops(), operator_kind="coefficient")
                assert_allclose(result.characters, chars, atol=1e-12)
                self.assertEqual(result.match.irrep, label)
                self.assertLess(max(rep.closure_error for rep in result.representations.values()), 1e-12)

    def test_random_complex_gauge_invariance(self):
        c = np.eye(3)[:, 1:]
        rng = np.random.default_rng(724)
        u, _ = np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))
        original = analyze_c3v(c, c3v_ops(), operator_kind="coefficient")
        rotated = analyze_c3v(c @ u, c3v_ops(), operator_kind="coefficient")
        assert_allclose(original.characters, rotated.characters, atol=1e-12)
        for name in c3v_ops():
            assert_allclose(rotated.representations[name].matrix,
                            u.conj().T @ original.representations[name].matrix @ u, atol=1e-12)

    def test_nonorthogonal_ao_and_both_operator_conventions(self):
        # A maps nonorthogonal AO coefficients to orthonormal physical coordinates.
        a = np.array([[1, .3j, .2], [0, 1.4, .1j], [0, 0, .7]], dtype=complex)
        s = a.conj().T @ a
        c = np.linalg.solve(a, np.eye(3)[:, 1:]) @ np.array([[2, .4j], [0, .5]])
        ops = {name: np.linalg.solve(a, op @ a) for name, op in c3v_ops().items()}
        for kind in ("coefficient", "matrix_element"):
            supplied = ops if kind == "coefficient" else {k: s @ t for k, t in ops.items()}
            result = analyze_c3v(c, supplied, s, operator_kind=kind)
            assert_allclose(result.characters, [2, -1, 0], atol=1e-12)
            self.assertEqual(result.match.irrep, "E")

    def test_single_e_component_is_not_closed(self):
        with self.assertRaisesRegex(ValueError, "not closed"):
            analyze_c3v(np.eye(3)[:, [1]], c3v_ops(), operator_kind="coefficient")

    def test_single_orbital_is_normalized(self):
        rep = project_representation([[3], [0]], np.eye(2), operator_kind="coefficient")
        self.assertAlmostEqual(rep.character.real, 1)

    def test_rejects_singular_subspace(self):
        with self.assertRaisesRegex(ValueError, "dependent"):
            orthonormalize([[1, 1], [0, 0]])

    def test_rejects_invalid_metric(self):
        for s in ([[1, 2], [0, 1]], [[1, 0], [0, -1]], [[1, 0], [0, 0]]):
            with self.subTest(s=s), self.assertRaises(ValueError):
                orthonormalize(np.eye(2), s)

    def test_check_only_rejects_unnormalized_orbitals(self):
        with self.assertRaisesRegex(ValueError, "S-orthonormal"):
            project_representation([[2]], [[1]], operator_kind="coefficient", orthonormalize_basis=False)

    def test_metric_preservation(self):
        ops = c3v_ops()
        ops["C3"] = 2 * ops["C3"]
        with self.assertRaisesRegex(ValueError, "preserve"):
            analyze_c3v(np.eye(3), ops, operator_kind="coefficient")

    def test_group_laws(self):
        ops = c3v_ops()
        ops["sv2"], ops["sv3"] = ops["sv3"], ops["sv2"]
        with self.assertRaisesRegex(ValueError, "group relations"):
            analyze_c3v(np.eye(3), ops, operator_kind="coefficient")

    def test_match_a2_reducible_and_noise(self):
        self.assertEqual(match_characters([1, 1, -1]).irrep, "A2")
        match = match_characters([3, 0, 1])
        self.assertTrue(match.valid)
        self.assertIsNone(match.irrep)
        self.assertEqual(match.multiplicities, {"A1": 1, "E": 1})
        self.assertEqual(match_characters([2+1e-7, -1-1e-7, 1e-7]).irrep, "E")

    def test_no_forced_match(self):
        for chars in ([2, -.7, 0], [1, -.5, 0], [2, -2, 0], [0, 0, 0], [-1, -1, -1], [2, -1, .1j]):
            with self.subTest(chars=chars):
                self.assertFalse(match_characters(chars).valid)

    def test_invalid_characters(self):
        for chars in ([2, -1], [np.nan, -1, 0]):
            with self.assertRaises(ValueError):
                match_characters(chars)
        with self.assertRaises(ValueError):
            match_characters([2, -1, 0], atol=-1)

    def test_class_consistency_and_missing_operations(self):
        chi = {key: np.trace(value) for key, value in c3v_ops().items()}
        chi["C3"] += .1
        with self.assertRaisesRegex(ValueError, "conjugacy"):
            class_characters(chi)
        del chi["C3"]
        with self.assertRaises(ValueError):
            class_characters(chi)

    def test_harmonic_order_and_unitarity(self):
        self.assertEqual(magnetic_order(2), (0, 1, -1, 2, -2))
        for ell in range(5):
            b = complex_to_real_basis(ell)
            assert_allclose(b.conj().T @ b, np.eye(2*ell+1), atol=1e-14)
        for ell in (-1, 1.5, True):
            with self.assertRaises(ValueError):
                magnetic_order(ell)

    def test_p_harmonic_quarter_turn_convention(self):
        # ABACUS symm_rotation_test.cpp C41 reference, m=(0,+1,-1).
        d = np.diag([1, -1j, 1j])
        assert_allclose(rotation_in_real_basis(d, 1), [[1, 0, 0], [0, 0, -1], [0, 1, 0]], atol=1e-14)

    def test_tda_complete_ph_space(self):
        # A1 occupied -> E virtual on an A2 reference gives total E.
        ops = c3v_ops()
        ph_ops = {name: particle_hole_operation(op[:1, :1], op[1:, 1:],
                  reference_character=(-1 if name.startswith("sv") else 1)) for name, op in ops.items()}
        result = analyze_c3v(np.eye(2), ph_ops, operator_kind="coefficient")
        self.assertEqual(result.match.irrep, "E")

    def test_ph_hole_conjugation_and_reference_phase(self):
        op = particle_hole_operation([[1j]], [[-1j]], reference_character=-1)
        assert_allclose(op, [[1]])


if __name__ == "__main__":
    unittest.main()
