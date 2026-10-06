"""Software-independent sparse R-space matrices and cell-gauge Fourier sums."""

from dataclasses import dataclass
import operator

import numpy as np
from scipy.sparse import csr_matrix


def _vectors(value, name):
    a = np.asarray(value)
    if a.dtype.kind not in "fiu" or a.ndim not in (1, 2) or a.shape[-1] != 3 or not a.size or not np.all(np.isfinite(a)):
        raise ValueError(f"{name} must be finite real coordinates of shape (3,) or (n,3)")
    return a.astype(float)


def k_cartesian_to_fractional(k_cartesian, lattice_vectors):
    """ABACUS k_cartesian (units 2pi/lat0) -> reciprocal fractional coordinates.

    lattice_vectors contains dimensionless DIRECT lattice vectors as ROWS;
    do not pass a cell in Angstrom without first dividing it by lat0.
    """
    k = _vectors(k_cartesian, "k_cartesian")
    a = _vectors(lattice_vectors, "lattice_vectors")
    if a.shape != (3, 3) or np.linalg.matrix_rank(a) != 3:
        raise ValueError("lattice_vectors must be a nonsingular 3x3 matrix")
    return k @ a.T


@dataclass(frozen=True)
class RealSpaceMatrix:
    """Full X(R), stored sparsely as (nR, nao*nao), AO columns flattened in C order.

    X(R) need not be Hermitian individually. All entries are retained, and
    neither absent R blocks nor lower triangles are filled by conjugation.
    Units are unchanged from the input (S dimensionless, ABACUS H normally Ry).
    """

    translations: np.ndarray
    values: csr_matrix
    basis_size: int
    ionic_step: int | None = None
    format: str = "arrays"
    spin_index: int | None = None
    lattice_vectors: np.ndarray | None = None
    lattice_constant_bohr: float | None = None
    source: str = ""

    def __post_init__(self):
        r = np.asarray(self.translations)
        try:
            size = operator.index(self.basis_size)
        except TypeError as exc:
            raise ValueError("basis_size must be a positive integer") from exc
        values = csr_matrix(self.values, dtype=complex, copy=True)
        if r.ndim != 2 or r.shape[1] != 3 or r.dtype.kind not in "iu" or len(np.unique(r, axis=0)) != len(r):
            raise ValueError("translations must be unique integer R vectors of shape (nR,3)")
        if isinstance(self.basis_size, bool) or size < 1 or values.shape != (len(r), size**2):
            raise ValueError("sparse matrix shape must be (nR,basis_size**2)")
        values.check_format(full_check=True)
        if not np.all(np.isfinite(values.data)):
            raise ValueError("matrix values must be finite")
        object.__setattr__(self, "translations", r.copy())
        object.__setattr__(self, "values", values)

    def to_k(self, kpoints):
        """X(k)=sum_R exp(+2pi i k.R) X(R), without weights or 1/nR.

        kpoints are reciprocal FRACTIONAL coordinates. A (3,) input returns
        (nao,nao); (nk,3) returns (nk,nao,nao). Only the result is dense.
        For large k meshes call this method on manageable batches.
        """
        k = _vectors(kpoints, "kpoints")
        phases = np.exp(2j * np.pi * (np.atleast_2d(k) @ self.translations.T))
        result = np.asarray(phases @ self.values).reshape(-1, self.basis_size, self.basis_size)
        return result[0] if k.ndim == 1 else result
