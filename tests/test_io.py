from pathlib import Path
import tempfile
import unittest

import numpy as np
from numpy.testing import assert_allclose

from stateid.io import AbacusReader, UnsupportedFormatError, read_gamma_wavefunctions, read_lr_npz, read_overlap_npy

FIXTURE = Path(__file__).parent / "fixtures" / "gamma_real.txt"


class IOTests(unittest.TestCase):
    def test_gamma_text(self):
        data = AbacusReader().read_wavefunctions(FIXTURE, spin="beta")
        assert_allclose(data.coefficients, np.eye(3))
        assert_allclose(data.energies_ry, [-1, .2, .2])
        assert_allclose(data.band_numbers, [1, 2, 3])
        self.assertEqual(data.spin, "beta")

    def test_gamma_rejects_corrupt_and_unsupported(self):
        original = FIXTURE.read_text()
        variants = [original + original, "1 (index of k points)\n0 0 0\n" + original,
                    original.replace("3 (band)", "4 (band)"),
                    original.replace("0.0 0.0 1.0", "0.0 0.0"),
                    original.replace("0.0 0.0 1.0", "0.0 nan 1.0"),
                    original.replace("0.0 0.0 1.0", "0.0 0.0 1.0 2.0")]
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "input.txt"
            for text in variants:
                path.write_text(text)
                with self.subTest(text=text[-40:]), self.assertRaises((ValueError, UnsupportedFormatError)):
                    read_gamma_wavefunctions(path)

    def test_overlap_npy(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "s.npy"
            s = np.array([[1, .1j], [-.1j, 1]])
            np.save(path, s)
            assert_allclose(read_overlap_npy(path), s)
            np.save(path, np.zeros((2, 2)))
            with self.assertRaises(ValueError):
                read_overlap_npy(path)

    def test_lr_exchange_and_full_lr(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "lr.npz"
            args = dict(schema_version=1, method="TDA", X=np.eye(2),
                        energies_ev=[2.1, 2.1], transitions=[[0, 1, 125, 1, 126], [0, 1, 125, 1, 127]])
            np.savez(path, **args)
            data = AbacusReader().read_lr_eigenvectors(path, format="stateid_npz")
            self.assertEqual(data.method, "TDA")
            self.assertIsNone(data.y)
            assert_allclose(data.x, np.eye(2))
            args.update(method="full_lr", Y=np.zeros((2, 2)))
            np.savez(path, **args)
            self.assertEqual(read_lr_npz(path).y.shape, (2, 2))

    def test_lr_rejects_bad_schema(self):
        base = dict(schema_version=1, method="TDA", X=np.eye(2), energies_ev=[2.1, 2.1],
                    transitions=[[0, 1, 125, 1, 126], [0, 1, 125, 1, 127]])
        edits = [dict(schema_version=2), dict(schema_version=1.5), dict(method="full_lr"),
                 dict(Y=np.zeros((2, 2))), dict(energies_ev=[2.1]), dict(X=np.ones((2, 2, 1))),
                 dict(transitions=[[0, 1, 125, 1, 126]] * 2),
                 dict(transitions=[[0, 1, 125, 1, 125], [0, 1, 125, 1, 127]])]
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "lr.npz"
            for edit in edits:
                np.savez(path, **(base | edit))
                with self.subTest(edit=edit), self.assertRaises(ValueError):
                    read_lr_npz(path)

    def test_native_placeholders_are_honest(self):
        reader = AbacusReader()
        for call in [lambda: reader.read_output("OUT"),
                     lambda: reader.read_overlap("sr_nao.csr", format="abacus_csr"),
                     lambda: reader.read_lr_eigenvectors("X.dat", format="abacus_lr"),
                     lambda: reader.read_wavefunctions("wfc.dat", format="binary")]:
            with self.assertRaises(UnsupportedFormatError):
                call()


if __name__ == "__main__":
    unittest.main()
