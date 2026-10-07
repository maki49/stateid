import importlib.util
import unittest
from unittest.mock import patch
import numpy as np
from numpy.testing import assert_allclose
from stateid.symmetry import (
    AOLabel, spgrep_irreps, build_bloch_ao_operation, analyze_little_group,
    decompose_characters, c3v_labels, match_irrep_labels,
    labels_from_character_table, irrep_labels, symmetry_from_structure,
)
from test_bloch import c3v_fractional, glide

HAS_SPGREP = importlib.util.find_spec('spgrep') is not None
HAS_IRREP = (importlib.util.find_spec('irrep') is not None
             and importlib.util.find_spec('irreptables') is not None)


class CharacterTests(unittest.TestCase):
    def test_decomposition_complex_noise_and_invalid_rows(self):
        rows={'plus':[1,-1j], 'minus':[1,1j]}
        result=decompose_characters([3,-1j],rows)
        self.assertEqual(result.multiplicities,{'plus':2,'minus':1})
        self.assertTrue(decompose_characters([1+1e-7,1j],rows).valid)
        for chi in ([1,.8j],[0,0],[-1,-1j],[1,np.nan]):
            if np.all(np.isfinite(chi)):
                self.assertFalse(decompose_characters(chi,rows).valid)
            else:
                with self.assertRaises(ValueError):
                    decompose_characters(chi,rows)
        with self.assertRaisesRegex(ValueError,'orthonormal'):
            decompose_characters([1,1],{'a':[1,1],'b':[1,1]})
        with self.assertRaises(ValueError):
            decompose_characters([1,1],{'a':[1,1]},weights=[1,-1])


@unittest.skipUnless(HAS_SPGREP,'optional stateid[symmetry] not installed')
class SpgrepTests(unittest.TestCase):
    def c3v(self, shuffle=False):
        w=c3v_fractional()
        if shuffle:
            w=w[[4,2,0,5,1,3]]
        ref=spgrep_irreps(w,np.zeros((6,3)),[0,0,0])
        a=np.array([[1,0,0],[-.5,np.sqrt(3)/2,0],[0,0,2]])
        labels=[AOLabel(0,'C',1,1,m) for m in (0,1,-1)]
        ops=[build_bloch_ao_operation([[0,0,0]],a,['C'],labels,r,[0,0,0])
             for r in ref.group.rotations]
        return ref,ops

    def test_c3v_e_a1_reducible_and_shuffled_operations(self):
        for shuffle in (False,True):
            ref,ops=self.c3v(shuffle)
            labels=c3v_labels(ref)
            for cols,expected in [([0],{'A1':1}),([1,2],{'E':1}),
                                   ([0,1,2],{'A1':1,'E':1})]:
                result=analyze_little_group(np.eye(3)[:,cols],ops,ref)
                named={labels.labels[key]:n for key,n in result.match.multiplicities.items()}
                self.assertEqual(named,expected)
                self.assertLess(result.group_error,1e-12)

    def test_random_complex_gauge_nonorthogonal_metric(self):
        ref,ops=self.c3v()
        a=np.array([[1,.2j,0],[0,1.3,.1j],[0,0,.8]],complex)
        s=a.conj().T@a
        transformed=[np.linalg.solve(a,op.matrix@a) for op in ops]
        c=np.linalg.solve(a,np.eye(3)[:,1:])
        rng=np.random.default_rng(621)
        u,_=np.linalg.qr(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)))
        first=analyze_little_group(c,transformed,ref,s)
        second=analyze_little_group(c@u,transformed,ref,s)
        assert_allclose(first.characters,second.characters,atol=1e-13)
        self.assertEqual(first.match.multiplicities,second.match.multiplicities)

    def test_closure_metric_group_and_k_failures(self):
        ref,ops=self.c3v()
        with self.assertRaisesRegex(ValueError,'closed'):
            analyze_little_group(np.eye(3)[:,[1]],ops,ref)
        wrong=[op.matrix.copy() for op in ops]
        wrong[1]*=2
        with self.assertRaisesRegex(ValueError,'metric'):
            analyze_little_group(np.eye(3),wrong,ref)
        wrong=[op.matrix for op in ops]
        wrong[4],wrong[5]=wrong[5],wrong[4]
        with self.assertRaisesRegex(ValueError,'relations'):
            analyze_little_group(np.eye(3),wrong,ref)
        with self.assertRaises(ValueError):
            analyze_little_group(np.eye(3),ops,ref,operator_kind='matrix_element')
        with self.assertRaises(ValueError):
            analyze_little_group(np.eye(3),ops[:-1],ref)

    def test_glide_projective_irreps_and_band_assignment(self):
        w=np.array([np.eye(3),np.diag([1,-1,1])],int)
        t=np.array([[0,0,0],[.5,0,0]])
        ref=spgrep_irreps(w,t,[.5,0,0])
        self.assertEqual(len(ref.matrices),2)
        assert_allclose(sorted(np.round([v[1].imag for v in ref.characters.values()])),[-1,1])
        b=glide([.5,0,0]).matrix
        # B [i,1]^T = i [i,1]^T, independent analytic solution.
        c=np.array([[1j],[1]])/np.sqrt(2)
        result=analyze_little_group(c,[np.eye(2),b],ref)
        self.assertTrue(result.match.valid)
        assert_allclose(result.characters,[1,1j],atol=1e-13)
        # Omitting the Bloch phase keeps closure for some vectors but breaks
        # the projective relation. Full space isolates that failure.
        with self.assertRaisesRegex(ValueError,'relations'):
            analyze_little_group(np.eye(2),[np.eye(2),[[0,1],[1,0]]],ref)

    def test_spgrep_preserves_integer_parts_including_identity(self):
        w=np.array([np.eye(3),np.diag([1,-1,1])],int)
        standard=spgrep_irreps(w,[[0,0,0],[.5,0,0]],[.25,0,0])
        shifted=spgrep_irreps(w,[[1,0,0],[1.5,0,0]],[.25,0,0])
        # Every representative gained a_x, hence every matrix gains -i.
        for a,b in zip(standard.matrices,shifted.matrices):
            assert_allclose(b,-1j*a,atol=1e-13)
        assert_allclose([chi[0] for chi in shifted.characters.values()],[-1j,-1j])
        pure=spgrep_irreps([np.eye(3)],[[1,0,0]],[.2,0,0])
        assert_allclose(pure.matrices[0],[[[np.exp(-.4j*np.pi)]]],atol=1e-13)

    def test_label_alignment_integer_translation_phase(self):
        w=np.array([np.eye(3),np.diag([1,-1,1])],int)
        ref=spgrep_irreps(w,[[0,0,0],[.5,0,0]],[.5,0,0])
        result=labels_from_character_table(ref,rotations=w[::-1],
            translations=[[1.5,0,0],[0,0,0]],kpoint=[.5,0,0],
            characters={'positive_branch':[-1j,1],'negative_branch':[1j,1]},
            source='analytic shifted glide table')
        positive=next(key for key,chi in ref.characters.items() if chi[1].imag>0)
        self.assertEqual(result.labels[positive],'positive_branch')
        self.assertLess(max(result.residuals.values()),1e-12)
        with self.assertRaisesRegex(ValueError,'origin'):
            labels_from_character_table(ref,rotations=w,
                translations=[[0,0,0],[.6,0,0]],kpoint=[.5,0,0],
                characters={'a':[1,1j],'b':[1,-1j]},source='wrong origin')
        with self.assertRaises(ValueError):
            labels_from_character_table(ref,rotations=w,
                translations=[[0,0,0],[.5,0,0]],kpoint=[-.5,0,0],
                characters={'a':[1,1j],'b':[1,-1j]},source='different branch')

    def test_no_labels_by_enumeration_or_ambiguous_match(self):
        ref,ops=self.c3v()
        chars=ref.characters
        bad={str(i):next(iter(chars.values())) for i in range(3)}
        with self.assertRaisesRegex(ValueError,'missing or ambiguous'):
            match_irrep_labels(ref,bad,source='bad')
        with self.assertRaises(ValueError):
            match_irrep_labels(ref,chars,source='')
        # Changing the irrep enumeration only changes local IDs, not names.
        from stateid.symmetry import IrrepSet
        flipped=IrrepSet(ref.group,ref.matrices[::-1],ref.backend)
        names=c3v_labels(flipped)
        for key,chi in flipped.characters.items():
            if abs(chi[0]-2)<1e-10:
                self.assertEqual(names.labels[key],'E')

    def test_structure_discovery_and_generic_k(self):
        lattice=np.array([[1.1,0,0],[.2,1.3,0],[.3,.4,1.7]])
        w,t=symmetry_from_structure(lattice,[[0,0,0]],[1])
        self.assertEqual(len(w),2)  # generic triclinic P-1
        ref=spgrep_irreps(w,t,[.123,.217,.319])
        self.assertEqual(ref.group.order,1)

    def test_missing_optional_dependency_is_explicit(self):
        import builtins
        original=builtins.__import__
        def guarded(name,*args,**kwargs):
            if name=='spgrep':
                raise ImportError('blocked for test')
            return original(name,*args,**kwargs)
        with patch('builtins.__import__',side_effect=guarded):
            with self.assertRaisesRegex(ImportError,'stateid\\[symmetry\\]'):
                self.c3v()


@unittest.skipUnless(HAS_SPGREP and HAS_IRREP,'optional stateid[labels] not installed')
class IrrepIntegrationTests(unittest.TestCase):
    def test_actual_nonsymmorphic_irrep_labels_at_boundary(self):
        a=np.array([[1,0,0],[0,1.3,0],[.2,0,1.7]])
        p=[[.13,.17,.23],[.63,-.17,.23],[.27,.32,.41],[.77,-.32,.41]]
        numbers=[1,1,2,2]
        w,t=symmetry_from_structure(a,p,numbers)
        self.assertEqual(len(w),2)  # Pc in an a-glide setting
        ref=spgrep_irreps(w,t,[-.5,0,0])
        named=irrep_labels(ref,a,p,numbers,'B')
        self.assertEqual(set(named.labels.values()),{'B1','B2'})
        self.assertIn('BCS SG 7',named.source)
        self.assertLess(max(named.residuals.values()),1e-12)
        assert_allclose(sorted(round(x[1].imag) for x in ref.characters.values()),[-1,1])

    def test_actual_irrep_tables_inversion_with_shifted_origin(self):
        lattice=np.array([[1.1,0,0],[.2,1.3,0],[.3,.4,1.7]])
        positions=[[.13,.21,.34]]
        w,t=symmetry_from_structure(lattice,positions,[1])
        ref=spgrep_irreps(w,t,[0,0,0])
        named=irrep_labels(ref,lattice,positions,[1],'GM')
        self.assertEqual(set(named.labels.values()),{'GM1+','GM1-'})
        self.assertIn('IrRep 3.3.0',named.source)
        for key,chi in ref.characters.items():
            inv=np.flatnonzero(np.all(w==-np.eye(3),axis=(1,2)))[0]
            self.assertEqual(named.labels[key], 'GM1+' if chi[inv].real>0 else 'GM1-')
        with self.assertRaisesRegex(ValueError,'lookup failed'):
            irrep_labels(ref,lattice,positions,[1],'not_a_kpoint')


if __name__ == '__main__':
    unittest.main()
