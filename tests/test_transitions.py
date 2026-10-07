import tempfile
from pathlib import Path
import unittest
from stateid.io.abacus_transitions import read_transition_analysis

HEADER = '''State Excitation Energy (Ry, eV) Transition dipole x, y, z Oscillator strength
0 .1 1.36 1 0 0 .1
1 .1 1.36 0 1 0 .1
State Occupied orbital Virtual orbital Excitation amplitude Excitation rate k-point
'''
class TransitionTests(unittest.TestCase):
    def read(self, text):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'analysis.dat'; p.write_text(text)
            return read_transition_analysis(p)
    def test_spin_continuation_and_missing_weight(self):
        states=self.read(HEADER+'0 126b 127b .8 .64 1\n126b 128b .5 .25 1\n1 126b 128b .9 .81 1\n')
        self.assertEqual(states[0].transitions[1].occupied,126)
        self.assertEqual(states[0].transitions[0].spin,'beta')
        self.assertAlmostEqual(states[0].summary()['reported_weight'],.89)
        self.assertEqual(states[1].index,1)
    def test_spin_flip_rejected(self):
        with self.assertRaises(ValueError): self.read(HEADER+'0 126a 127b .8 .64 1\n')
    def test_bad_weight_rejected(self):
        with self.assertRaises(ValueError): self.read(HEADER+'0 126b 127b .8 .2 1\n')
    def test_empty_state_allowed(self):
        self.assertEqual(self.read(HEADER)[0].summary()['reported_weight'],0)
    def test_incomplete_file_rejected(self):
        with self.assertRaises(ValueError): self.read('not an analysis file')

class RealNVTests(unittest.TestCase):
    def test_endpoints(self):
        fixture=Path(__file__).parent/'fixtures'
        avg=read_transition_analysis(fixture/'nv_avg_tda.dat')
        jt=read_transition_analysis(fixture/'nv_jt_tda.dat')
        self.assertEqual(len(avg),6)
        self.assertAlmostEqual(avg[0].summary()['reported_weight'],.963884)
        self.assertAlmostEqual(avg[1].summary()['reported_weight'],.964171)
        self.assertEqual(jt[0].transitions[0].virtual,127)
        self.assertEqual(jt[1].transitions[0].virtual,128)
        self.assertAlmostEqual(jt[0].summary()['reported_weight'],.966687)

class PolarizationTests(unittest.TestCase):
    def test_known_direction_and_scale(self):
        from stateid.io.abacus_transitions import TransitionState
        state=TransitionState(0,.1,1.36,(1.,-1.,0.),.1,())
        self.assertAlmostEqual(state.longitudinal_fraction([1,1,1]),0.)
        self.assertAlmostEqual(state.longitudinal_fraction([2,-2,0]),1.)
        with self.assertRaises(ValueError): state.longitudinal_fraction([0,0,0])
    def test_dark_state_rejected(self):
        from stateid.io.abacus_transitions import TransitionState
        state=TransitionState(0,.1,1.36,(0.,0.,0.),0.,())
        with self.assertRaises(ValueError): state.longitudinal_fraction([1,1,1])
