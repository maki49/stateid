"""Read ABACUS's thresholded transition analysis, not raw LR eigenvectors.

State indices are zero-based; orbital and k-point labels remain one-based.
Printed weights are not renormalized: omitted small amplitudes are real data.
Spin a/b denotes alpha/beta; no spin multiplicity or orbital irrep is inferred.
"""
from dataclasses import asdict, dataclass
from pathlib import Path
import argparse
import json
import re
import numpy as np

@dataclass(frozen=True)
class Transition:
    occupied: int
    virtual: int
    spin: str
    kpoint: int
    amplitude: float
    weight: float

@dataclass(frozen=True)
class TransitionState:
    index: int
    energy_ry: float
    energy_ev: float
    dipole_au: tuple[float, float, float]
    oscillator_strength: float
    transitions: tuple[Transition, ...]

    def longitudinal_fraction(self, axis):
        """Fraction |d dot axis_hat|² / |d|²; undefined for a dark state."""
        axis = np.asarray(axis, dtype=float)
        if axis.shape != (3,) or not np.all(np.isfinite(axis)) or np.linalg.norm(axis) == 0:
            raise ValueError('axis must be a finite nonzero three-vector')
        d = np.asarray(self.dipole_au)
        if np.dot(d, d) == 0:
            raise ValueError('polarization fraction is undefined for a zero dipole')
        return float(abs(np.dot(d, axis / np.linalg.norm(axis)))**2 / np.dot(d, d))

    def summary(self):
        weights = {}
        for t in self.transitions:
            weights[t.spin] = weights.get(t.spin, 0.0) + t.weight
        return {**asdict(self), 'reported_weight': sum(weights.values()),
                'reported_spin_weights': weights,
                'thresholded': True, 'normalized_by_parser': False}

def read_transition_analysis(path):
    """Read real-TDA trans_analysis_*.dat; reject unsupported complex/LR data.

    The caller must establish TDA from the calculation settings. ABACUS's
    full-LR analysis currently prints X only, so that file is not a CI weight.
    """
    metadata, transitions = {}, {}
    section, current = None, None
    for lineno, line in enumerate(Path(path).read_text().splitlines(), 1):
        if 'State' in line and 'Excitation Energy' in line:
            section = 'energies'; continue
        if 'State' in line and 'Occupied orbital' in line:
            section = 'transitions'; continue
        tokens = line.split()
        if not tokens or not re.match(r'^\d', tokens[0]):
            continue
        try:
            if section == 'energies':
                if len(tokens) != 7: raise ValueError('expected seven energy fields')
                index = int(tokens[0])
                if index in metadata: raise ValueError('duplicate state')
                vals = tuple(float(v) for v in tokens[1:])
                if not np.all(np.isfinite(vals)): raise ValueError('nonfinite metadata')
                metadata[index] = vals
                transitions[index] = []
            elif section == 'transitions':
                if len(tokens) == 6:
                    current = int(tokens.pop(0))
                if len(tokens) != 5 or current not in metadata:
                    raise ValueError('invalid transition row or missing state')
                occupied, virtual = [re.fullmatch(r'(\d+)([ab]?)', t) for t in tokens[:2]]
                if occupied is None or virtual is None or occupied[2] != virtual[2]:
                    raise ValueError('invalid orbital labels or spin-flip transition')
                spin = {'a': 'alpha', 'b': 'beta', '': 'unspecified'}[occupied[2]]
                amplitude, weight = map(float, tokens[2:4])
                kpoint = int(tokens[4])
                if min(int(occupied[1]), int(virtual[1]), kpoint) < 1:
                    raise ValueError('band and k-point labels must be positive')
                if not np.isfinite(amplitude) or not np.isfinite(weight) or weight < 0:
                    raise ValueError('invalid amplitude or weight')
                if not np.isclose(weight, amplitude**2, rtol=2e-5, atol=2e-6):
                    raise ValueError('printed weight inconsistent with real TDA amplitude')
                t = Transition(int(occupied[1]), int(virtual[1]), spin, kpoint, amplitude, weight)
                if any((p.occupied,p.virtual,p.spin,p.kpoint)==(t.occupied,t.virtual,t.spin,t.kpoint)
                       for p in transitions[current]): raise ValueError('duplicate transition')
                transitions[current].append(t)
        except (ValueError, KeyError) as exc:
            raise ValueError(f'{path}:{lineno}: {exc}') from exc
    if not metadata or section != 'transitions':
        raise ValueError('missing energy/transition table')
    return tuple(TransitionState(i, v[0], v[1], v[2:5], v[5], tuple(transitions[i]))
                 for i, v in sorted(metadata.items()))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('paths', nargs='+')
    parser.add_argument('--states', type=int, nargs='+', default=[0, 1])
    parser.add_argument('--axis', type=float, nargs=3, help='NV or other reference axis in Cartesian coordinates')
    args = parser.parse_args()
    result = {p: [s.summary() for s in read_transition_analysis(p) if s.index in args.states]
              for p in args.paths}
    if args.axis is not None:
        for p, summaries in result.items():
            states = {s.index: s for s in read_transition_analysis(p)}
            for summary in summaries:
                state = states[summary['index']]
                summary['longitudinal_fraction'] = state.longitudinal_fraction(args.axis)
    print(json.dumps(result, indent=2, ensure_ascii=False))

if __name__ == '__main__':
    main()
