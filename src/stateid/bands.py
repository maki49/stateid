"""Select complete candidate band-edge subspaces at one sampled k point."""
from dataclasses import dataclass
import numpy as np
from ._validation import positive_tolerance


@dataclass(frozen=True)
class BandEdges:
    valence_indices: np.ndarray  # zero-based coefficient columns
    conduction_indices: np.ndarray
    valence_energy: float
    conduction_energy: float

    @property
    def gap(self):
        return self.conduction_energy - self.valence_energy


def select_band_edges(energies, occupations, *, full_occupation,
                      energy_tol=1e-6, occupation_tol=1e-6):
    """Select local valence/conduction edges from one k/spin block.

    energy_tol uses the SAME unit as energies. Return every occupied/unoccupied
    level within energy_tol of its respective edge, without assuming band order.
    These are candidates; symmetry closure must be checked independently.

    full_occupation explicitly gives the raw value for a fully occupied band:
    for ABACUS wg it is spin degeneracy times this k weight (2*w for scalar
    nspin=1, w per channel for nspin=2). No scale is inferred from the data.
    Fractionally occupied or unresolved-gap cases fail rather than assigning
    insulator band edges by an arbitrary threshold. A zero-weight k point
    requires independently supplied occupations; wg alone carries no evidence.
    """
    positive_tolerance(energy_tol, 'energy_tol')
    positive_tolerance(occupation_tol, 'occupation_tol')
    if occupation_tol >= .5:
        raise ValueError('occupation_tol must be less than 0.5')
    scale = positive_tolerance(full_occupation, 'full_occupation')
    e = np.asarray(energies, dtype=float)
    o = np.asarray(occupations, dtype=float)
    if e.ndim != 1 or not len(e) or o.shape != e.shape or not np.all(np.isfinite(e)) or not np.all(np.isfinite(o)):
        raise ValueError('energies and occupations must be matching nonempty finite vectors')
    fractions = o / scale
    occupied = np.abs(fractions - 1) <= occupation_tol
    empty = np.abs(fractions) <= occupation_tol
    if not np.all(occupied | empty):
        raise ValueError('partial or invalid occupations: insulating band edges are ambiguous')
    if not np.any(occupied) or not np.any(empty):
        raise ValueError('both occupied and unoccupied bands are required')
    ev, ec = float(np.max(e[occupied])), float(np.min(e[empty]))
    if ec - ev <= energy_tol:
        raise ValueError('occupied/unoccupied gap is nonpositive or unresolved within energy_tol')
    return BandEdges(np.flatnonzero(occupied & (ev-e <= energy_tol)),
                     np.flatnonzero(empty & (e-ec <= energy_tol)), ev, ec)
