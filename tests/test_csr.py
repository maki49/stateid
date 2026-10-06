from pathlib import Path
import tempfile
import unittest

import numpy as np
from numpy.testing import assert_allclose

from stateid.io import AbacusReader, iter_csr, read_csr, read_wavefunctions
from stateid.realspace import RealSpaceMatrix, k_cartesian_to_fractional


def sample_blocks():
    hopping = np.array([[0, .1+.2j], [0, 0]])
    return {(0, 0, 0): np.array([[1.5, .05j], [-.05j, 2]]),
            (1, 0, 0): hopping, (-1, 0, 0): hopping.conj().T,
            (0, 1, 0): np.zeros((2, 2))}


def csr_text(blocks, *, new=False, step=0, empty_arrays=True):
    """Synthetic writer for the two documented layouts, not production code."""
    size = next(iter(blocks.values())).shape[0]
    if new:
        lines = [f"--- Ionic Step {step} ---", "# print S matrix in real space S(R)",
                 "1 # number of spin directions", "1 # spin index",
                 f"{size} # number of localized basis", f"{len(blocks)} # number of Bravais lattice vector R",
                 "user_defined_lattice", "10.0", "1 0 0", ".2 1 0", "0 .3 1",
                 "X", "1", "Direct", "0 0 0", "# CSR Format"]
    else:
        lines = [f"STEP: {step}", f"Matrix Dimension of S(R): {size}", f"Matrix number of S(R): {len(blocks)}"]
    for r, matrix in blocks.items():
        rows, cols = np.nonzero(matrix)
        values = matrix[rows, cols]
        lines.append(" ".join(map(str, (*r, len(values)))))
        if not len(values) and not (new and empty_arrays):
            continue
        if new:
            lines.append("# CSR values")
        # Deliberately wrap and include spaces inside C++ complex literals.
        lines.extend(f"({value.real:.16e}, {value.imag:.16e})" if np.iscomplexobj(matrix)
                     else f"{value:.16e}" for value in values)
        if new:
            lines.append("# CSR column indices")
        lines.append(" ".join(map(str, cols)))
        if new:
            lines.append("# CSR row pointers")
        pointers = np.r_[0, np.cumsum(np.bincount(rows, minlength=size))]
        lines.append(" ".join(map(str, pointers)))
    return "\n".join(lines) + "\n"


class CSRTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "S.csr"

    def load(self, text, **kwargs):
        self.path.write_text(text)
        return read_csr(self.path, **kwargs)

    def test_old_new_equal_and_metadata(self):
        old = self.load(csr_text(sample_blocks()))
        new = self.load(csr_text(sample_blocks(), new=True, step=3))
        assert_allclose(old.values.toarray(), new.values.toarray())
        assert_allclose(old.translations, new.translations)
        self.assertEqual(old.format, "abacus_legacy")
        self.assertEqual(new.format, "abacus_new")
        self.assertEqual(new.ionic_step, 3)
        self.assertEqual(new.spin_index, 1)
        self.assertEqual(new.lattice_constant_bohr, 10)

    def test_generic_matrix_from_arrays(self):
        data = RealSpaceMatrix([[0, 0, 0], [1, 0, 0]], [[2], [.2j]], 1)
        assert_allclose(data.to_k([.25, 0, 0]), [[1.8]], atol=1e-14)
        for size in (True, 1.5, -1):
            with self.assertRaises(ValueError):
                RealSpaceMatrix([[0, 0, 0]], [[1]], size)

    def test_phase_sign_full_matrix_and_batch(self):
        data = self.load(csr_text(sample_blocks(), new=True))
        ks = np.array([[0, 0, 0], [.25, .13, 0], [.37, 0, 0]])
        actual = data.to_k(ks)
        for k, sk in zip(ks, actual):
            upper = .05j + (.1+.2j) * np.exp(2j*np.pi*k[0])
            assert_allclose(sk, [[1.5, upper], [upper.conjugate(), 2]], atol=1e-14)
            assert_allclose(data.to_k(k), sk, atol=1e-14)
        assert_allclose(data.to_k(ks + [1, -2, 3]), actual, atol=1e-14)

    def test_no_implicit_hermitian_completion(self):
        data = self.load(csr_text({(1, 0, 0): np.array([[0., .1], [0., 0.]])}))
        assert_allclose(data.to_k([.25, 0, 0]), [[0, .1j], [0, 0]], atol=1e-15)

    def test_real_values_and_fortran_exponent(self):
        text = csr_text({(0, 0, 0): np.eye(2)}).replace("e+00", "D+00")
        assert_allclose(self.load(text).to_k([0, 0, 0]), np.eye(2))

    def test_empty_blocks_in_both_new_variants(self):
        for labels in (False, True):
            with self.subTest(labels=labels):
                data = self.load(csr_text(sample_blocks(), new=True, empty_arrays=labels))
                self.assertEqual(data.values.getrow(3).nnz, 0)

    def test_multiframe_requires_selection(self):
        for new in (False, True):
            text = csr_text(sample_blocks(), new=new, step=7) + csr_text(sample_blocks(), new=new, step=11)
            with self.subTest(new=new), self.assertRaisesRegex(ValueError, "multiple CSR"):
                self.load(text)
            self.assertEqual(self.load(text, frame=1).ionic_step, 11)
            self.assertEqual(len(list(iter_csr(self.path))), 2)
            with self.assertRaises(ValueError):
                read_csr(self.path, frame=2)

    def test_no_legacy_step_header(self):
        text = csr_text(sample_blocks()).split("\n", 1)[1]
        self.assertIsNone(self.load(text).ionic_step)

    def test_corrupt_csr_rejected(self):
        valid = csr_text({(0, 0, 0): np.eye(2)})
        variants = [valid.replace("0 1\n", "0 2\n"),  # out-of-range column
                    valid.replace("0 1 2\n", "0 2 1\n"),
                    valid.replace("0 1 2\n", "1 1 2\n"),
                    valid.replace("0 0 0 2\n", "0 0 0 -2\n"),
                    valid.replace("1.0000000000000000e+00", "nan", 1),
                    valid.rsplit("\n", 2)[0], valid + "unexpected trailing data\n",
                    valid.replace("Matrix number of S(R): 1", "Matrix number of S(R): 2")]
        for text in variants:
            with self.subTest(text=text), self.assertRaises(ValueError):
                self.load(text)

    def test_duplicate_r_rejected(self):
        text = csr_text(sample_blocks()).replace("1 0 0 1\n", "0 0 0 1\n", 1)
        with self.assertRaisesRegex(ValueError, "duplicate"):
            self.load(text)

    def test_invalid_k_and_frame(self):
        data = self.load(csr_text(sample_blocks()))
        for k in ([0, 0], [[0, 0, np.nan]], [], [0, 0, 1j]):
            with self.subTest(k=k), self.assertRaises(ValueError):
                data.to_k(k)
        for frame in (-1, True, .5):
            with self.assertRaises(ValueError):
                read_csr(self.path, frame=frame)

    def test_cartesian_fractional_conversion_nonorthogonal_cell(self):
        a = np.array([[1, .2, 0], [0, 1, .3], [.1, 0, 1]])
        kfrac = np.array([[.1, .2, -.3], [.25, 0, .3]])
        kcart = kfrac @ np.linalg.inv(a).T
        assert_allclose(k_cartesian_to_fractional(kcart, a), kfrac, atol=1e-14)
        with self.assertRaises(ValueError):
            k_cartesian_to_fractional([0, 0, 0], np.zeros((3, 3)))

    def test_overlap_api_and_orthonormality_at_multiple_k(self):
        self.load(csr_text(sample_blocks(), new=True))
        ks = [[0, 0, 0], [.25, 0, 0], [.37, 0, 0]]
        sks = AbacusReader().read_overlap(self.path, format="abacus_csr", kpoints=ks)
        self.assertEqual(sks.shape, (3, 2, 2))
        for sk in sks:
            eigenvalues, u = np.linalg.eigh(sk)
            c = u / np.sqrt(eigenvalues)
            assert_allclose(c.conj().T @ sk @ c, np.eye(2), atol=1e-14)
        assert_allclose(AbacusReader().read_overlap(self.path, format="abacus_csr"), sks[0])

    def test_overlap_rejects_nonhermitian_or_indefinite(self):
        for matrix in (np.array([[1, .2], [0, 1]]), np.diag([1, -1])):
            self.load(csr_text({(0, 0, 0): matrix}))
            with self.assertRaises(ValueError):
                AbacusReader().read_overlap(self.path, format="abacus_csr")

    def test_overlap_npy_cannot_be_fourier_transformed(self):
        with self.assertRaises(ValueError):
            AbacusReader().read_overlap("unused.npy", format="stateid_npy", kpoints=[0, 0, 0])

    def test_multik_wavefunctions(self):
        text = ("2 (index of k points)\n0.25 0.1 -0.2\n2 (number of bands)\n"
                "2 (number of orbitals)\n1 (band)\n-0.1 (Ry)\n0.5 (Occupations)\n"
                "1 2 3\n4\n2 (band)\n0.2 (Ry)\n0 (Occupations)\n5 6 7 8\n")
        self.path.write_text(text)
        data = read_wavefunctions(self.path, spin="beta")
        self.assertEqual(data.k_index, 2)
        assert_allclose(data.k_cartesian, [.25, .1, -.2])
        assert_allclose(data.coefficients, [[1+2j, 5+6j], [3+4j, 7+8j]])
        for invalid in (text.replace("5 6 7 8", "5 6 7"), text + text,
                        text.replace("0.25 0.1 -0.2", "0.25 nan -0.2"),
                        text.replace("2 (index of k points)", "0 (index of k points)")):
            self.path.write_text(invalid)
            with self.assertRaises((ValueError, NotImplementedError)):
                read_wavefunctions(self.path)


if __name__ == "__main__":
    unittest.main()
