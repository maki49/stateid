import tempfile
from pathlib import Path
import unittest

import numpy as np
from numpy.testing import assert_allclose

from stateid.io import read_ao_labels
from stateid.symmetry import (AOLabel, build_gamma_ao_operation,
                              real_harmonic_operation, magnetic_order, analyze_c3v)


def labels_for(natoms, shells=(0, 1, 2)):
    return [AOLabel(atom, 'C', ell, zeta, m) for atom in range(natoms)
            for zeta in (1, 2) for ell in shells for m in magnetic_order(ell)]


class AOTests(unittest.TestCase):
    def test_orbital_actual_atom_and_component_indices(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'Orbital'
            path.write_text('#io spec l m z sym\n0 C 1 0 2 pz\n'
                            '0 C 1 1 2 px\n0 C 1 2 2 py\n'
                            '#io = Orbital index in supercell\n')
            result = read_ao_labels(path)
            self.assertEqual([r.atom_index for r in result], [0, 0, 0])
            self.assertEqual([r.m for r in result], [0, 1, -1])
            self.assertEqual([r.zeta for r in result], [2, 2, 2])
            for text in ('', '0 C 1 -1 1 px', '0 C 1 3 1 px',
                         '0 C 1 0 0 px', '0 C 1 0 1', 'bad C 1 0 1 px'):
                path.write_text(text)
                with self.assertRaises(ValueError):
                    read_ao_labels(path)

    def test_active_p_rotation_and_inverse_negative(self):
        r = np.array([[0,-1,0], [1,0,0], [0,0,1]])
        # R(-x)=-y; columns in z,-x,-y order.
        expected = np.array([[1,0,0], [0,0,-1], [0,1,0]])
        assert_allclose(real_harmonic_operation(r, 1), expected, atol=1e-14)
        self.assertGreater(np.linalg.norm(real_harmonic_operation(r.T, 1)-expected), 1)

    def test_direct_d_function_values(self):
        r = np.array([[0,-1,0], [1,0,0], [0,0,1]])
        t = real_harmonic_operation(r, 2)
        # z2 invariant, xz -> yz, yz -> -xz, x2-y2 and xy -> negatives.
        expected = np.array([[1,0,0,0,0], [0,0,-1,0,0], [0,1,0,0,0],
                             [0,0,0,-1,0], [0,0,0,0,-1]])
        assert_allclose(t, expected, atol=1e-14)
        p = np.array([.23, -.74, .91])
        def d(v):
            x,y,z = v
            return np.array([(2*z*z-x*x-y*y)/2, -np.sqrt(3)*x*z,
                             -np.sqrt(3)*y*z, np.sqrt(3)*(x*x-y*y)/2,
                             np.sqrt(3)*x*y])
        assert_allclose(d(p) @ t, d(r.T @ p), atol=1e-14)

    def test_noncommuting_products_and_improper_parity(self):
        r = np.array([[0,-1,0], [1,0,0], [0,0,1]])
        f = np.diag([-1,1,1])
        for ell in (0,1,2):
            a, b = (real_harmonic_operation(op, ell) for op in (r,f))
            assert_allclose(real_harmonic_operation(r @ f, ell), a @ b, atol=1e-14)
            assert_allclose(real_harmonic_operation(-r, ell), (-1)**ell*a, atol=1e-14)
            assert_allclose(a.T @ a, np.eye(2*ell+1), atol=1e-14)
        self.assertGreater(np.linalg.norm(real_harmonic_operation(r @ f,1)
                                          - real_harmonic_operation(f @ r,1)), 1)

    def test_periodic_mapping_multiple_shells_arbitrary_ao_order(self):
        positions = [[.25,0,0], [.75,0,0]]
        labels = labels_for(2)[::-1]
        result = build_gamma_ao_operation(positions, np.eye(3), ['C','C'], labels,
                                          np.diag([-1,1,1]))
        assert_allclose(result.atom_mapping, [1,0])
        assert_allclose(result.lattice_shifts, [[-1,0,0], [-1,0,0]])
        assert_allclose(result.matrix @ result.matrix, np.eye(len(labels)), atol=1e-14)
        for col, label in enumerate(labels):
            nonzero = np.flatnonzero(abs(result.matrix[:, col]) > 1e-10)
            for row in nonzero:
                other = labels[row]
                self.assertEqual(other.atom_index, 1-label.atom_index)
                self.assertEqual((other.ell, other.zeta), (label.ell, label.zeta))
        # Nonorthogonal, invariant radial overlap; check full AO metric directly.
        s = np.eye(len(labels))
        for i, a in enumerate(labels):
            for j, b in enumerate(labels):
                if (a.atom_index, a.ell, a.m) == (b.atom_index, b.ell, b.m) and a.zeta != b.zeta:
                    s[i,j] = .2
        assert_allclose(result.matrix.T @ s @ result.matrix, s, atol=1e-14)

    def test_translation_and_reported_approximation(self):
        result = build_gamma_ao_operation([[0,0,0]], np.eye(3), ['C'], labels_for(1),
                                          np.eye(3), translation=[1+1e-7,0,0])
        assert_allclose(result.lattice_shifts, [[1,0,0]])
        self.assertAlmostEqual(result.atom_errors[0], 1e-7)
        self.assertEqual(result.symprec, 1e-6)
        with self.assertRaises(ValueError):
            build_gamma_ao_operation([[0,0,0]], np.eye(3), ['C'], labels_for(1),
                                     np.eye(3), translation=[1+1e-7,0,0], symprec=1e-8)

    def test_full_c3v_nonorthogonal_metric(self):
        c,s = -.5, np.sqrt(3)/2
        r = np.array([[c,-s,0], [s,c,0], [0,0,1]])
        f = np.diag([1,-1,1])
        lattice = np.array([[1,0,0], [-.5,np.sqrt(3)/2,0], [0,0,2]])
        labels = labels_for(1, (1,))
        ops = {name: build_gamma_ao_operation([[0,0,0]], lattice, ['C'], labels, op).matrix
               for name,op in {'E': np.eye(3), 'C3':r, 'C3^2':r@r,
                               'sv1':f, 'sv2':r@f, 'sv3':r@r@f}.items()}
        metric = np.kron([[1,.2],[.2,1.4]], np.eye(3))
        result = analyze_c3v(np.eye(6)[:,[1,2]], ops, metric, operator_kind='coefficient')
        self.assertEqual(result.match.irrep, 'E')
        self.assertLess(result.group_error, 1e-12)

    def test_invalid_mapping_shells_and_geometry(self):
        def build(pos=[[0,0,0]], species=['C'], labels=None, rotation=None, lattice=None):
            return build_gamma_ao_operation(pos, np.eye(3) if lattice is None else lattice,
                       species, labels_for(1) if labels is None else labels,
                       np.eye(3) if rotation is None else rotation)
        for labels in ([], labels_for(1)[:-1], labels_for(1)+labels_for(1)[:1],
                       [AOLabel(0,'N',0,1,0)], [AOLabel(0,'C',1,1,True)]):
            with self.subTest(labels=labels), self.assertRaises(ValueError):
                build(labels=labels)
        with self.assertRaisesRegex(ValueError, 'ambiguous'):
            build(pos=[[0,0,0],[1,0,0]], species=['C','C'], labels=labels_for(2))
        with self.assertRaises(ValueError):
            build(pos=[[.2,0,0]], rotation=-np.eye(3))
        with self.assertRaises(ValueError):
            build(rotation=np.eye(3)*2)
        with self.assertRaises(ValueError):
            build(lattice=np.zeros((3,3)))
        with self.assertRaisesRegex(ValueError, 'periodic lattice'):
            build(rotation=np.array([[.6,-.8,0],[.8,.6,0],[0,0,1]]))
        with self.assertRaises(NotImplementedError):
            real_harmonic_operation(np.eye(3), 3)
        with self.assertRaisesRegex(ValueError, 'different radial shells'):
            build(pos=[[.25,0,0],[.75,0,0]], species=['C','C'],
                  labels=labels_for(1)+[AOLabel(1,'C',0,1,0)], rotation=-np.eye(3))


if __name__ == '__main__':
    unittest.main()
