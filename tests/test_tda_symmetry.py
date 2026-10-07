import unittest
from unittest.mock import patch
import numpy as np
from numpy.testing import assert_allclose
from stateid.symmetry import analyze_tda_c3v, analyze_c3v, particle_hole_operation


def setup():
    c,s = -.5,np.sqrt(3)/2
    r = np.array([[c,-s],[s,c]])
    f = np.diag([1,-1])
    e = {'E':np.eye(2), 'C3':r, 'C3^2':r@r,
         'sv1':f, 'sv2':r@f, 'sv3':r@r@f}
    # Reference alpha occupied E gives A2 determinant, beta occupied A1.
    ref = {name:[op, np.eye(1)] for name,op in e.items()}
    # Unequal spin blocks (2x2 and 1x2); beta roots give the E pair.
    ops = {name:[(op,op), (np.eye(1),op)] for name,op in e.items()}
    x = np.zeros((6,2), dtype=complex)
    x[4:] = np.eye(2)
    return x,ops,ref


class TDATests(unittest.TestCase):
    def test_dense_equivalence_and_reference_a2(self):
        x,ops,ref = setup()
        dense = {}
        for name,blocks in ops.items():
            phase = np.linalg.det(ref[name][0])
            a,b = [particle_hole_operation(o,v,reference_character=phase) for o,v in blocks]
            dense[name] = np.block([[a,np.zeros((4,2))],[np.zeros((2,4)),b]])
        result = analyze_tda_c3v(x,ops,reference_occupied_operations=ref)
        expected = analyze_c3v(x,dense,operator_kind='coefficient')
        self.assertEqual(result.symmetry.match.irrep,'E')
        for name in ops:
            assert_allclose(result.symmetry.representations[name].matrix,
                            expected.representations[name].matrix,atol=1e-14)
            self.assertAlmostEqual(result.reference_characters[name].real,
                                   -1 if name.startswith('sv') else 1)

    def test_complex_gauges_and_two_spin_blocks(self):
        x,ops,ref = setup()
        rng = np.random.default_rng(1729)
        def unitary(n):
            return np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))[0]
        # Complete ph space ensures the nontrivial first spin contraction is
        # tested too; complex orbital gauges require hole conjugation.
        full = np.eye(6,dtype=complex)
        a,b,c,d = unitary(2),unitary(2),unitary(1),unitary(2)
        ops = {name:[(a.conj().T@pairs[0][0]@a,b.conj().T@pairs[0][1]@b),
                     (c.conj().T@pairs[1][0]@c,d.conj().T@pairs[1][1]@d)]
               for name,pairs in ops.items()}
        result = analyze_tda_c3v(full,ops,reference_occupied_operations=ref)
        for name,pairs in ops.items():
            phase = result.reference_characters[name]
            a,b = [particle_hole_operation(o,v,reference_character=phase) for o,v in pairs]
            expected = np.block([[a,np.zeros((4,2))],[np.zeros((2,4)),b]])
            assert_allclose(result.symmetry.representations[name].matrix,expected,atol=1e-14)
        x,ops,ref = setup()
        u = unitary(2)
        rotated = analyze_tda_c3v(x@u,ops,reference_occupied_operations=ref)
        assert_allclose(rotated.symmetry.characters,[2,-1,0],atol=1e-14)
        self.assertLess(max(r.closure_error for r in rotated.symmetry.representations.values()),1e-14)

    def test_does_not_construct_kron_or_ph_identity(self):
        x,ops,ref = setup()
        with patch('numpy.kron',side_effect=AssertionError('dense ph operator')):
            analyze_tda_c3v(x,ops,reference_occupied_operations=ref)

    def test_reference_phase_changes_total_irrep(self):
        # An A1 excitation on an A2 reference is total A2, not A1.
        _,_,ref = setup()
        ops = {name:[(np.eye(1),np.eye(1))] for name in ref}
        result = analyze_tda_c3v([[1]],ops,reference_occupied_operations=ref)
        self.assertEqual(result.symmetry.match.irrep,'A2')
        empty_beta = {name:[blocks[0],np.empty((0,0))] for name,blocks in ref.items()}
        self.assertEqual(analyze_tda_c3v([[1]],ops,reference_occupied_operations=empty_beta).symmetry.match.irrep,'A2')

    def test_rejects_single_e_root_and_bad_norm_layout_group(self):
        x,ops,ref = setup()
        for roots,pattern in ((x[:,:1],'not closed'),(2*x,'not orthonormal'),(x[:-1],'layout')):
            with self.assertRaisesRegex(ValueError,pattern):
                analyze_tda_c3v(roots,ops,reference_occupied_operations=ref)
        bad = dict(ops)
        bad['sv2'],bad['sv3'] = bad['sv3'],bad['sv2']
        with self.assertRaisesRegex(ValueError,'group relations'):
            analyze_tda_c3v(x,bad,reference_occupied_operations=ref)
        bad = {name:[np.eye(1)*2,np.eye(1)] for name in ref}
        with self.assertRaises(ValueError):
            analyze_tda_c3v(x,ops,reference_occupied_operations=bad)
        with self.assertRaises(ValueError):
            analyze_tda_c3v(x,ops,reference_occupied_operations={})
        with self.assertRaises(ValueError):
            analyze_tda_c3v(x,ops,reference_occupied_operations={name:[np.eye(1)] for name in ref})


if __name__ == '__main__':
    unittest.main()
