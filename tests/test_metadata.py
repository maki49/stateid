from pathlib import Path
import tempfile
import unittest

import numpy as np
from numpy.testing import assert_allclose

from stateid.io import (read_structure, read_input_parameters, read_kpoints_info,
                        read_scf_log, read_multik_wavefunctions, UnsupportedFormatError)
from stateid.io.metadata import ABACUS_BOHR_TO_ANGSTROM

LATTICE = np.array([[2., 0, 0], [.3, 3, 0], [.2, .4, 4]])
POSITIONS = np.array([[.125, .25, .5], [.5, -.25, 1.125]])
STRU = """ATOMIC_SPECIES
C 12 C.upf
N 14 N.upf
NUMERICAL_ORBITAL
C.orb
N.orb
LATTICE_CONSTANT
2
LATTICE_VECTORS
2 0 0
.3 3 0
.2 .4 4
ATOMIC_POSITIONS
{mode}
C
0
1
{c} 1 0 1
N
1
1
{n} m 1 1 1 mag .7 f 0 0 0
"""


class MetadataTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "input"

    def write(self, text):
        self.path.write_text(text)
        return self.path

    def stru(self, mode="Direct"):
        coords = POSITIONS if mode == "Direct" else POSITIONS @ LATTICE
        coords = coords * {"Direct": 1, "Cartesian": 1, "Cartesian_au": 2,
                           "Cartesian_angstrom": 2*ABACUS_BOHR_TO_ANGSTROM}[mode]
        return STRU.format(mode=mode, c=" ".join(map(str, coords[0])),
                           n=" ".join(map(str, coords[1])))

    def test_coordinate_units_atom_order_and_fields(self):
        for mode in ("Direct", "Cartesian", "Cartesian_au", "Cartesian_angstrom"):
            cell = read_structure(self.write("# comment\n" + self.stru(mode)))
            assert_allclose(cell.fractional_positions, POSITIONS, atol=1e-14)
            assert_allclose(cell.lattice_angstrom, 2*ABACUS_BOHR_TO_ANGSTROM*LATTICE)
            self.assertEqual(cell.species, ("C", "N"))
            self.assertEqual(cell.species_ids.tolist(), [1, 2])
            self.assertEqual(cell.atom_fields[1][:4], ("m", "1", "1", "1"))

    def test_structure_rejects_unsupported_or_invalid(self):
        text = self.stru()
        variants = [text.replace("Direct", "Cartesian_angstrom_center_xy"),
                    text.replace("LATTICE_VECTORS", "LATTICE_PARAMETERS"),
                    text.replace(".2 .4 4", ".2 .4 0"),
                    text.replace("0.125", "nan"), text.replace("C 12", "C -12"),
                    text.replace("N\n1\n1", "C\n1\n1"),
                    text[:text.rfind("\nN\n")], text+"\nextra",
                    text+text, text.replace("LATTICE_CONSTANT\n2", "LATTICE_CONSTANT\n0")]
        for value in variants:
            with self.subTest(value=value[-35:]), self.assertRaises((ValueError, UnsupportedFormatError)):
                read_structure(self.write(value))

    def test_input_info_empty_values_and_comments(self):
        values = read_input_parameters(self.write(
            "INPUT_PARAMETERS\n# emitted metadata\nnspin 2 # spins\ncalculation scf\nesolver_type ksdft\nread_file_dir # unset\n"))
        self.assertEqual(values["nspin"], "2")
        self.assertEqual(values["read_file_dir"], "")
        with self.assertRaises(ValueError):
            read_input_parameters(self.write("INPUT_PARAMETERS\nnspin 2\nnspin 1"))
        with self.assertRaises(UnsupportedFormatError):
            read_input_parameters(self.write("nspin 2"))

    def kinfo(self):
        return ("nkstot now = 2\nK-POINTS DIRECT COORDINATES\n"
                "KPOINTS DIRECT_X DIRECT_Y DIRECT_Z WEIGHT\n"
                "1 .25 0 0 .3333\n2 -.25 0 0 .6666\n")

    def test_kpoint_metadata_keeps_raw_rounded_weights(self):
        data = read_kpoints_info(self.write(self.kinfo()))
        assert_allclose(data.fractional, [[.25, 0, 0], [-.25, 0, 0]])
        self.assertAlmostEqual(data.weights.sum(), .9999)  # no implicit normalization
        for text in (self.kinfo().replace("2 -.25", "1 -.25"),
                     self.kinfo().replace(".6666", "-.1"),
                     self.kinfo().replace(".3333", "nan"), self.kinfo()+self.kinfo()):
            with self.assertRaises((ValueError, UnsupportedFormatError)):
                read_kpoints_info(self.write(text))

    def test_scf_convergence_and_multiple_steps(self):
        text = ("ABACUS v3.11.0-beta10\nCommit: abc123\n"
                "#ION MOVE# 1 #ELEC ITER# 24\n#SCF IS CONVERGED#\n"
                "!FINAL_ETOT_IS -10.25 eV\n")
        data = read_scf_log(self.write(text))
        self.assertTrue(data["converged"])
        self.assertEqual(data["iterations"], 24)
        self.assertEqual(data["energy_ev"], -10.25)
        self.assertFalse(read_scf_log(self.write(text.replace("#SCF IS CONVERGED#", "")))["converged"])
        for extra in (text, "#ION MOVE# 2 #ELEC ITER# 1\n"):
            with self.assertRaises(UnsupportedFormatError):
                read_scf_log(self.write(text+extra))

    def test_modern_collection_checks_spin_global_index_and_k(self):
        points = read_kpoints_info(self.write(self.kinfo()))
        directory = Path(self.temp.name)
        for spin in (1, 2):
            for ik, point in enumerate(points.fractional, 1):
                cart = point / 2  # dimensionless lattice = 2 I
                (directory/f"wfk{ik}s{spin}g1_nao.txt").write_text(
                    f"{(spin-1)*2+ik} (index of k points)\n"
                    + " ".join(map(str, cart)) + "\n1 (number of bands)\n1 (number of orbitals)\n"
                    "1 (band)\n-.3 (Ry)\n.5 (Occupations)\n1 0\n")
        kwargs = dict(kpoints=points, lattice_vectors=2*np.eye(3), nspin=2)
        data = read_multik_wavefunctions(directory, **kwargs)
        self.assertEqual(data[1][0].spin, "beta")
        self.assertEqual(data[1][0].k_index, 3)
        with self.assertRaises(ValueError):
            read_multik_wavefunctions(directory, **(kwargs | {"lattice_vectors": np.eye(3)}))
        bad = directory/"wfk1s2g1_nao.txt"
        original = bad.read_text()
        bad.write_text(original.replace("3 (index", "1 (index"))
        with self.assertRaises(ValueError):
            read_multik_wavefunctions(directory, **kwargs)
        bad.write_text(original)
        for nspin in (4, 1.5, True, 1.0):
            with self.assertRaises(UnsupportedFormatError):
                read_multik_wavefunctions(directory, **(kwargs | {"nspin": nspin}))

    def test_scalar_collection_uses_no_spin_filename_suffix(self):
        points = read_kpoints_info(self.write(self.kinfo()))
        directory = Path(self.temp.name)
        for ik, point in enumerate(points.fractional, 1):
            (directory/f"wfk{ik}g1_nao.txt").write_text(
                f"{ik} (index of k points)\n" + " ".join(map(str, point))
                + "\n1 (number of bands)\n1 (number of orbitals)\n"
                  "1 (band)\n-.3 (Ry)\n.25 (Occupations)\n1 0\n")
        data = read_multik_wavefunctions(directory, kpoints=points,
                                        lattice_vectors=np.eye(3), nspin=1)
        self.assertEqual(len(data), 1)
        self.assertEqual(data[0][0].spin, "unspecified")
        self.assertEqual(data[0][1].k_index, 2)
