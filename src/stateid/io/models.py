"""Backend-neutral data contracts. Arrays are dense in the first milestone."""

from dataclasses import dataclass
from typing import Any, Protocol

import numpy as np


@dataclass(frozen=True)
class OrbitalData:
    coefficients: np.ndarray  # (nao, nband), columns are orbitals
    energies_ry: np.ndarray
    occupations: np.ndarray  # raw ABACUS occupation/weight field, not reinterpreted
    band_numbers: np.ndarray  # one-based file labels; NumPy column index is label-1
    spin: str  # alpha, beta, or unspecified; never inferred from filename
    source: str


@dataclass(frozen=True)
class LRData:
    x: np.ndarray  # (nph, nstate)
    y: np.ndarray | None
    energies_ev: np.ndarray
    transitions: np.ndarray  # (nph,5): k, spin_hole, occupied, spin_particle, virtual
    method: str  # TDA or full_lr
    source: str


class OutputReader(Protocol):
    """Future backends implement explicit methods instead of guessing formats."""

    def read_output(self, path) -> dict[str, Any]: ...
    def read_wavefunctions(self, path, *, format: str, spin: str) -> OrbitalData: ...
    def read_overlap(self, path, *, format: str) -> np.ndarray: ...
    def read_lr_eigenvectors(self, path, *, format: str) -> LRData: ...
