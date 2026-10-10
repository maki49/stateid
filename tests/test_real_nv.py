"""Regression from unmodified local ABACUS non-Gamma NV outputs, not toy C/S."""
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from numpy.testing import assert_allclose

from stateid.io import (read_structure, read_csr, read_wavefunctions,
                        read_ao_labels, read_kpoints_info, read_scf_log)
from stateid.realspace import k_cartesian_to_fractional
from stateid.symmetry import (build_bloch_ao_operation, symmetry_from_structure,
    spgrep_irreps, analyze_little_group, c3v_labels, project_representation,
    project_sewing_matrix)

FIXTURE = Path(__file__).parent/"fixtures/nv_non_gamma"
HAS_SYMMETRY = all(importlib.util.find_spec(p) is not None for p in ("spgrep", "spglib"))


class RealNVTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temp.cleanup)
        cls.directory = Path(cls.temp.name)
        manifest = json.loads((FIXTURE/"provenance.json").read_text())
        for name, checksum in manifest["fixture_raw_sha256"].items():
            packed = FIXTURE/(name+".gz")
            content = gzip.decompress(packed.read_bytes()) if packed.exists() else (FIXTURE/name).read_bytes()
            if hashlib.sha256(content).hexdigest() != checksum:
                raise AssertionError(f"fixture integrity failure: {name}")
            (cls.directory/name).write_bytes(content)
        cls.cell = read_structure(cls.directory/"STRU")
        cls.k = read_kpoints_info(cls.directory/"KPT.info").fractional
        cls.sr = read_csr(cls.directory/"srg1_nao.csr")
        cls.s = cls.sr.to_k(cls.k)
        cls.labels = read_ao_labels(cls.directory/"Orbital")
        cls.data = [read_wavefunctions(cls.directory/f"wfk{i}s2g1_nao.txt", spin="beta") for i in (1, 2, 3)]

    def setup_symmetry(self):
        c = self.cell
        w, t = symmetry_from_structure(c.lattice_angstrom, c.fractional_positions, c.species_ids)
        ref = spgrep_irreps(w, t, self.k[0])
        ops = [build_bloch_ao_operation(c.fractional_positions, c.lattice_angstrom, c.species,
               self.labels, r, self.k[0], translation_fractional=tr)
               for r, tr in zip(ref.group.rotations, ref.group.translations)]
        return ref, ops

    def test_actual_metadata_csr_units_and_norms(self):
        final = read_structure(self.directory/"STRU_FINAL")
        assert_allclose(final.lattice_angstrom, self.cell.lattice_angstrom, atol=1e-9, rtol=0)
        assert_allclose(final.fractional_positions, self.cell.fractional_positions, atol=1e-10, rtol=0)
        self.assertAlmostEqual(self.sr.lattice_constant_angstrom, 1.)
        assert_allclose(self.sr.lattice_vectors*self.sr.lattice_constant_angstrom,
                        self.cell.lattice_angstrom, atol=4e-6, rtol=0)
        for ik, d in enumerate(self.data):
            self.assertEqual(d.k_index, 9+ik)
            assert_allclose(k_cartesian_to_fractional(d.k_cartesian, self.cell.lattice_vectors),
                            self.k[ik], atol=1e-8, rtol=0)
            assert_allclose(d.coefficients.conj().T@self.s[ik]@d.coefficients,
                            np.eye(91), atol=4e-8, rtol=0)
        log = read_scf_log(self.directory/"running_scf.log")
        self.assertTrue(log["converged"])
        self.assertEqual(log["iterations"], 24)

    def test_wrong_fourier_sign_is_detectable(self):
        c = self.data[1].coefficients
        wrong_s = self.sr.to_k(-self.k[1])
        self.assertGreater(np.max(abs(c.conj().T@wrong_s@c-np.eye(91))), .1)

    @unittest.skipUnless(HAS_SYMMETRY, "optional spgrep/spglib")
    def test_actual_a1_e_and_degenerate_gauge(self):
        ref, ops = self.setup_symmetry()
        names = c3v_labels(ref)
        for cols, label in (([13], "A1"), ([14, 15], "E")):
            c = self.data[0].coefficients[:, cols]
            result = analyze_little_group(c, ops, ref, self.s[0], atol=1e-6)
            self.assertEqual(names.labels[result.match.irrep], label)
            self.assertLess(max(r.closure_error for r in result.representations), 1e-7)
            if label == "E":
                # Character table in the actual operation order: E, C3, C3^2, mirrors.
                assert_allclose(result.characters, [2, -1, -1, 0, 0, 0], atol=1e-10)
                rng = np.random.default_rng(724)
                u, _ = np.linalg.qr(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)))
                rotated = analyze_little_group(c@u, ops, ref, self.s[0], atol=1e-6)
                assert_allclose(rotated.characters, result.characters, atol=1e-12)
                with self.assertRaisesRegex(ValueError, "not closed"):
                    analyze_little_group(c[:, :1], ops, ref, self.s[0], atol=1e-6)

    @unittest.skipUnless(HAS_SYMMETRY, "optional spgrep/spglib")
    def test_actual_cross_k_sewing(self):
        ref, _ = self.setup_symmetry()
        w = next(w for w in ref.group.rotations
                 if np.allclose(np.linalg.solve(w.T, self.k[1]), self.k[2]))
        cell = self.cell
        op = build_bloch_ao_operation(cell.fractional_positions, cell.lattice_angstrom,
                                     cell.species, self.labels, w, self.k[1])
        result = project_sewing_matrix(self.data[1].coefficients[:, :16],
                     self.data[2].coefficients[:, :16], op.matrix, self.s[1], self.s[2])
        self.assertLess(result.ao_metric_error, 1e-8)
        self.assertLess(result.closure_error, 1e-8)
        assert_allclose(result.matrix.conj().T@result.matrix, np.eye(16), atol=1e-8)

    @unittest.skipUnless(HAS_SYMMETRY, "optional spgrep/spglib")
    def test_real_data_periodic_image_gauge_requires_bloch_phase(self):
        ref, ops = self.setup_symmetry()
        cell, k, s = self.cell, self.k[0], self.s[0]
        shifts = np.zeros_like(cell.fractional_positions)
        shifts[0] = [1, 0, 0]
        v = np.exp(2j*np.pi*(shifts[[a.atom_index for a in self.labels]]@k))
        c = self.data[0].coefficients[:, 14:16]
        new_c, new_s = v[:, None]*c, v[:, None]*s*v.conj()[None, :]
        new_ops = [build_bloch_ao_operation(cell.fractional_positions+shifts,
                    cell.lattice_angstrom, cell.species, self.labels, w, k,
                    translation_fractional=t)
                    for w, t in zip(ref.group.rotations, ref.group.translations)]
        new = analyze_little_group(new_c, new_ops, ref, new_s, atol=1e-6)
        original = analyze_little_group(c, ops, ref, s, atol=1e-6)
        assert_allclose(new.characters, original.characters, atol=1e-12)
        for old, op in zip(ops, new_ops):
            assert_allclose(op.matrix, v[:, None]*old.matrix*v.conj()[None, :], atol=1e-12)
        wrong = [project_representation(new_c, op.ao_operation.matrix, new_s,
                    operator_kind="coefficient") for op in new_ops]
        self.assertGreater(max(r.closure_error for r in wrong), .5)
        # Finite-k translations cannot borrow the zero-translation label table.
        shifted = spgrep_irreps(ref.group.rotations, ref.group.translations+[1,0,0], k)
        with self.assertRaises(ValueError):
            c3v_labels(shifted)
