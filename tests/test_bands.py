import unittest
import numpy as np
from numpy.testing import assert_array_equal
from stateid.bands import select_band_edges


class BandEdgeTests(unittest.TestCase):
    def test_unsorted_degenerate_edges_and_explicit_weight_scale(self):
        e = [2, -.2, -3, -.2+2e-8, 2+3e-8, 4]
        for scale in (2, 1, .25, .125):
            result = select_band_edges(e, np.array([0,1,1,1,0,0])*scale,
                                       full_occupation=scale)
            assert_array_equal(result.valence_indices, [1,3])
            assert_array_equal(result.conduction_indices, [0,4])
            self.assertAlmostEqual(result.gap, 2.2-2e-8)

    def test_cluster_uses_distance_from_edge_without_chain_extension(self):
        result = select_band_edges([-1.0000015,-1.0000008,-1,1], [1,1,1,0],
                                   full_occupation=1)
        assert_array_equal(result.valence_indices, [1,2])

    def test_partial_occupation_and_wrong_scale_rejected(self):
        for o,scale in (([1,.1,0],1),([.25,.25,0],1),([1,-.1,0],1)):
            with self.assertRaisesRegex(ValueError, 'occupations'):
                select_band_edges([-1,0,1],o,full_occupation=scale)

    def test_missing_edges_and_unresolved_gap_rejected(self):
        for e,o in (([-1,1],[1,1]),([-1,1],[0,0]),([0,0],[1,0]),
                    ([1,-1],[1,0]),([0,1e-8],[1,0])):
            with self.assertRaises(ValueError):
                select_band_edges(e,o,full_occupation=1)

    def test_finite_vector_and_tolerance_validation(self):
        for e,o in (([],[]),([1,2],[1]),([np.nan,2],[1,0]),([1,2],[1,np.inf]),([[1,2]],[[1,0]])):
            with self.assertRaises(ValueError):
                select_band_edges(e,o,full_occupation=1)
        for kwargs in ({'full_occupation':0},{'full_occupation':np.nan},
                       {'full_occupation':1,'occupation_tol':.5},
                       {'full_occupation':1,'energy_tol':0}):
            with self.assertRaises(ValueError):
                select_band_edges([-1,1],[1,0],**kwargs)
