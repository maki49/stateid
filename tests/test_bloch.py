import unittest
import numpy as np
from numpy.testing import assert_allclose
from stateid.symmetry import (
    AOLabel, build_bloch_ao_operation, bloch_from_ao_operation,
    build_gamma_ao_operation, find_little_group, transform_kpoint,
    project_representation, project_sewing_matrix,
)


def c3v_fractional():
    r = np.array([[0,-1,0], [1,-1,0], [0,0,1]])
    f = np.array([[1,-1,0], [0,-1,0], [0,0,1]])
    return np.array([np.eye(3,dtype=int), r, r@r, f, r@f, r@r@f])


def glide(k):
    return build_bloch_ao_operation(
        [[0,.17,0], [.5,-.17,0]], np.eye(3), ['C','C'],
        [AOLabel(i,'C',0,1,0) for i in range(2)],
        np.diag([1,-1,1]), k, translation_fractional=[.5,0,0])


class BlochTests(unittest.TestCase):
    def test_glide_analytic_columns_and_target_k(self):
        result = glide([.25,.21,0])
        assert_allclose(result.mapped_kpoint, [.25,-.21,0])
        assert_allclose(result.matrix, [[0,-1j],[1,0]], atol=1e-14)
        assert_allclose(result.ao_operation.lattice_shifts, [[0,0,0],[1,0,0]])
        # Independent real-space Bloch-sum relabelling, R' = W R + L_source.
        w = np.diag([1,-1,1])
        k = np.array([.25,.21,0])
        for source in range(2):
            target = result.ao_operation.atom_mapping[source]
            shift = result.ao_operation.lattice_shifts[source]
            for r in ([0,0,0], [2,-3,1], [-1,4,0]):
                lhs = np.exp(2j*np.pi*np.dot(k,r))
                rp = w @ r + shift
                rhs = result.matrix[target,source]*np.exp(
                    2j*np.pi*np.dot(result.mapped_kpoint,rp))
                assert_allclose(lhs, rhs, atol=1e-14)

    def test_glide_squared_is_translation_not_identity(self):
        b = glide([.37,0,0]).matrix
        assert_allclose(b @ b, np.exp(-2j*np.pi*.37)*np.eye(2), atol=1e-14)
        assert_allclose(glide([1.37,0,0]).matrix, b, atol=1e-14)

    def test_gamma_limit_and_reuse(self):
        result = glide([0,0,0])
        assert_allclose(result.matrix, result.ao_operation.matrix, atol=1e-14)
        cached = bloch_from_ao_operation(result.ao_operation, [0,1],
                                        [.25,.1,0], np.diag([1,-1,1]))
        assert_allclose(cached.matrix, glide([.25,.1,0]).matrix)

    def test_full_translation_phase(self):
        b = build_bloch_ao_operation([[.2,.1,.3]], np.eye(3), ['C'],
            [AOLabel(0,'C',0,1,0)], np.eye(3), [.2,.3,.4],
            translation_fractional=[1,-2,0])
        assert_allclose(b.matrix, [[np.exp(-2j*np.pi*(.2-.6))]])

    def test_nonorthogonal_lattice_reciprocal_action_and_harmonics(self):
        a = np.array([[1,0,0],[-.5,np.sqrt(3)/2,0],[0,0,2]])
        w = c3v_fractional()[1]
        k = np.array([.23,.14,.09])
        kp = transform_kpoint(k,w)
        self.assertGreater(np.linalg.norm(kp-w@k), .1)
        for r in ([2,1,0],[-1,2,1]):
            self.assertAlmostEqual(k @ r, kp @ (w @ r))
        labels = [AOLabel(0,'C',1,1,m) for m in (0,1,-1)]
        b = build_bloch_ao_operation([[0,0,0]],a,['C'],labels,w,k)
        cart = a.T @ w @ np.linalg.inv(a.T)
        gamma = build_gamma_ao_operation([[0,0,0]],a,['C'],labels,cart)
        assert_allclose(b.matrix, gamma.matrix)

    def test_little_group_generic_axis_and_boundary(self):
        w = c3v_fractional()
        self.assertEqual(find_little_group(w,np.zeros((6,3)),[0,0,.27]).order,6)
        self.assertEqual(find_little_group(w,np.zeros((6,3)),[.123,.217,.31]).order,1)
        inv = np.array([np.eye(3),-np.eye(3)],int)
        little = find_little_group(inv,np.zeros((2,3)),[.5,0,0])
        assert_allclose(little.reciprocal_shifts, [[0,0,0],[-1,0,0]])
        self.assertEqual(find_little_group(inv,np.zeros((2,3)),[.2,0,0]).order,1)

    def test_glide_factor_system_and_nonzero_identity_representative(self):
        w = np.array([np.eye(3),np.diag([1,-1,1])],int)
        t = np.array([[0,0,0],[.5,0,0]])
        g = find_little_group(w,t,[.5,0,0])
        assert_allclose(g.multiplication,[[0,1],[1,0]])
        assert_allclose(g.factor_system,[[1,1],[1,-1]],atol=1e-14)
        t[0]=[1,0,0]
        g = find_little_group(w,t,[.5,0,0])
        assert_allclose(g.factor_system,[[-1,-1],[-1,1]],atol=1e-14)

    def test_bad_geometry_k_and_operation_sets(self):
        w = c3v_fractional()
        for ops, trans in [(w[:2],np.zeros((2,3))),
                           (np.array([np.eye(3),np.eye(3)]),[[0,0,0],[.5,0,0]]),
                           (w,np.full((6,3),.123))]:
            with self.assertRaises(ValueError):
                find_little_group(ops,trans,[0,0,0])
        for k in ([np.nan,0,0],[0,0],[0,0,1j]):
            with self.assertRaises(ValueError):
                glide(k)
        with self.assertRaises(ValueError):
            transform_kpoint([0,0,0],np.zeros((3,3)))
        with self.assertRaises(ValueError):
            bloch_from_ao_operation(glide([0,0,0]).ao_operation,[0,2],
                                    [0,0,0],np.eye(3))
        # Antiunitarity is not inferred by a sign on k or a spin flag.

    def test_target_k_phase_for_inversion(self):
        k=np.array([.21,.13,.07])
        b=build_bloch_ao_operation(
            [[.2,.1,.3],[.8,.9,.7]],np.eye(3),['C','C'],
            [AOLabel(i,'C',0,1,0) for i in range(2)],-np.eye(3),k)
        phase=np.exp(-2j*np.pi*k.sum())
        assert_allclose(b.matrix,phase*np.array([[0,1],[1,0]]),atol=1e-14)
        self.assertGreater(abs(phase-phase.conjugate()),.5)

    def test_noncommuting_cross_k_composition(self):
        r=np.array([[0,-1,0],[1,0,0],[0,0,1]])
        f=np.diag([-1,1,1])
        labels=[AOLabel(0,'C',ell,1,m)
                for ell,ms in [(0,[0]),(1,[0,1,-1]),(2,[0,1,-1,2,-2])]
                for m in ms]
        def build(w,k):
            return build_bloch_ao_operation([[.5,.5,0]],np.eye(3),['C'],labels,w,k)
        k=[.173,.291,.07]
        right=build(f,k)
        left=build(r,right.mapped_kpoint)
        product=build(r@f,k)
        assert_allclose(left.matrix@right.matrix,product.matrix,atol=1e-13)
        self.assertGreater(np.linalg.norm(build(f@r,k).matrix-product.matrix),1)

    def test_metric_covariance_from_independent_periodic_gaussian_overlaps(self):
        from itertools import product
        positions=np.array([[0,.17,0],[.5,-.17,0]])
        cells=np.array(list(product(range(-4,5),repeat=3)))
        # Direct Gaussian AO overlaps <mu,0|nu,R>, not generated from B.
        sr=np.exp(-4*np.sum((cells[:,None,None,:]
            +positions[None,None,:,:]-positions[None,:,None,:])**2,axis=-1))
        def overlap(k):
            return np.einsum('r,rij->ij',np.exp(2j*np.pi*(cells@k)),sr)
        k=np.array([.23,.31,.17])
        b=glide(k)
        assert_allclose(b.matrix.conj().T@overlap(b.mapped_kpoint)@b.matrix,
                        overlap(k),atol=1e-13)

    def test_two_metric_sewing_and_independent_gauges(self):
        a = np.array([[1,.2j],[0,1.3]],complex)
        at = np.array([[.8,.1],[0,1.2]],complex)
        u = np.array([[0,1j],[1,0]])
        s, st = a.conj().T@a, at.conj().T@at
        c, ct = np.linalg.inv(a), np.linalg.inv(at)
        b = np.linalg.solve(at,u@a)
        for kind,op in [('coefficient',b),('matrix_element',st@b)]:
            result=project_sewing_matrix(c,ct,op,s,st,operator_kind=kind)
            assert_allclose(result.matrix,u,atol=1e-13)
            self.assertLess(result.ao_metric_error,1e-13)
            self.assertLess(result.closure_error,1e-13)
            self.assertFalse(hasattr(result,'character'))
        us=np.diag([1j,1]); ut=np.array([[0,1],[1,0]])
        gauged=project_sewing_matrix(c@us,ct@ut,b,s,st)
        assert_allclose(gauged.matrix,ut.conj().T@u@us,atol=1e-13)
        self.assertGreater(abs(np.trace(gauged.matrix)-np.trace(u)),.1)
        same=project_sewing_matrix(c,c,np.eye(2),s,s)
        old=project_representation(c,np.eye(2),s,operator_kind='coefficient')
        assert_allclose(same.matrix,old.matrix)
        bad=project_sewing_matrix(c,ct,2*b,s,st)
        self.assertGreater(bad.ao_metric_error,1)
        with self.assertRaises(ValueError):
            project_sewing_matrix(c,ct[:,:1],b,s,st)


if __name__ == '__main__':
    unittest.main()
