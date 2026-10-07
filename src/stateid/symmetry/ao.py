"""Explicit Gamma-point AO operations in the ABACUS real s/p/d basis.

Lattice vectors are rows; Cartesian rotations act on columns. No symmetry
finder, k phase, spinor, or approximate-symmetry repair is implicit here.
"""
from dataclasses import dataclass
import operator

import numpy as np

from stateid._validation import matrix, positive_tolerance
from .harmonics import magnetic_order


@dataclass(frozen=True)
class AOLabel:
    """One coefficient row: zero-based atom, signed m, one-based zeta."""
    atom_index: int
    species: str
    ell: int
    zeta: int
    m: int


@dataclass(frozen=True)
class AOOperation:
    matrix: np.ndarray
    atom_mapping: np.ndarray  # source atom -> target atom
    lattice_shifts: np.ndarray  # transformed source = target + shift @ lattice
    atom_errors: np.ndarray  # Cartesian distances, in the supplied lattice unit
    symprec: float


def _integer(value, name, minimum=0):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be an integer")
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be an integer") from exc
    if result < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    return result


def _values(points, ell):
    x, y, z = points.T
    if ell == 0:
        return np.ones((len(points), 1))
    if ell == 1:
        return np.column_stack((z, -x, -y))
    return np.column_stack(((2*z*z-x*x-y*y)/2, -np.sqrt(3)*x*z,
                            -np.sqrt(3)*y*z, np.sqrt(3)*(x*x-y*y)/2,
                            np.sqrt(3)*x*y))


def real_harmonic_operation(rotation, ell, *, atol=1e-10):
    """T for (R f)(r)=f(R^-1 r), including improper parity, for ell<=2.

    Independent Cartesian homogeneous polynomials implement the documented
    complex_to_real_basis convention; common normalization cancels in T.
    """
    positive_tolerance(atol)
    order = magnetic_order(ell)
    ell = (len(order)-1)//2
    if ell > 2:
        raise NotImplementedError("AO rotations currently support only s/p/d shells")
    r = matrix(rotation, "rotation")
    if r.shape != (3, 3) or np.max(np.abs(r.imag)) > atol:
        raise ValueError("rotation must be a real 3x3 orthogonal Cartesian matrix")
    r = r.real
    if not np.allclose(r.T @ r, np.eye(3), atol=atol, rtol=0):
        raise ValueError("rotation must be orthogonal")
    # Overdetermined exact polynomial interpolation, with a fixed full-rank set.
    points = np.array([[1,0,0], [0,1,0], [0,0,1], [1,1,0],
                       [1,0,1], [0,1,1], [1,-1,1]], dtype=float)
    return np.linalg.lstsq(_values(points, ell), _values(points @ r, ell), rcond=None)[0]


def build_gamma_ao_operation(fractional_positions, lattice_vectors, species, labels,
                             rotation, *, translation=None, symprec=1e-6, atol=1e-10):
    """Build columns R|phi_nu> from explicit geometry and AO row labels.

    translation is Cartesian in the lattice unit, and rotation is about the
    Cartesian origin. AO labels must enumerate every complete radial shell;
    shell ordering may be arbitrary. Matching is periodic, species-preserving,
    bijective, and unique within symprec. No C/S or electronic invariance is
    inferred: subsequently validate T†ST and orbital closure with analyze_c3v.
    """
    symprec = positive_tolerance(symprec, "symprec")
    positive_tolerance(atol)
    positions = matrix(fractional_positions, "fractional_positions")
    lattice = matrix(lattice_vectors, "lattice_vectors")
    if positions.shape[1] != 3 or not len(positions) or lattice.shape != (3, 3):
        raise ValueError("positions must be (nat,3) and lattice must be (3,3)")
    if np.any(positions.imag) or np.any(lattice.imag):
        raise ValueError("geometry must be real")
    positions, lattice = positions.real, lattice.real
    if np.linalg.matrix_rank(lattice) != 3:
        raise ValueError("lattice must be nonsingular")
    species, labels = tuple(species), tuple(labels)
    if len(species) != len(positions) or any(not isinstance(s, str) or not s for s in species):
        raise ValueError("species must contain one nonempty label per atom")
    r = matrix(rotation, "rotation")
    real_harmonic_operation(r, 0, atol=atol)  # validates even for all-s shells
    r = r.real
    # A crystal operation must also preserve the periodic translation lattice.
    lattice_map = lattice @ r.T @ np.linalg.inv(lattice)
    if not np.allclose(lattice_map, np.rint(lattice_map), atol=atol, rtol=0):
        raise ValueError("rotation does not preserve the periodic lattice")
    shift = np.zeros(3) if translation is None else np.asarray(translation, dtype=float)
    if shift.shape != (3,) or not np.all(np.isfinite(shift)):
        raise ValueError("translation must be a finite Cartesian 3-vector")
    transformed = (positions @ lattice @ r.T + shift) @ np.linalg.inv(lattice)
    mapping, shifts, errors = [], [], []
    for atom, pos in enumerate(transformed):
        delta = pos - positions
        integers = np.rint(delta)
        distances = np.linalg.norm((delta-integers) @ lattice, axis=1)
        candidates = [j for j in range(len(species))
                      if species[j] == species[atom] and distances[j] <= symprec]
        if len(candidates) != 1:
            raise ValueError(f"atom {atom}: periodic mapping missing or ambiguous within symprec")
        target = candidates[0]
        mapping.append(target)
        shifts.append(integers[target].astype(int))
        errors.append(distances[target])
    if len(set(mapping)) != len(mapping):
        raise ValueError("atom mapping is not bijective")
    shells = {}
    for index, label in enumerate(labels):
        atom = _integer(label.atom_index, "atom_index")
        ell = _integer(label.ell, "ell")
        zeta = _integer(label.zeta, "zeta", 1)
        if atom >= len(species) or label.species != species[atom]:
            raise ValueError("AO species/atom does not match geometry")
        _integer(abs(label.m), "abs(m)")
        if isinstance(label.m, (bool, np.bool_)) or label.m not in magnetic_order(ell):
            raise ValueError("invalid signed magnetic index")
        shell = shells.setdefault((atom, ell, zeta), {})
        if label.m in shell:
            raise ValueError("duplicate AO within a radial shell")
        shell[label.m] = index
    if not shells or set(key[0] for key in shells) != set(range(len(species))):
        raise ValueError("AO labels must cover every atom")
    for (atom, ell, zeta), shell in shells.items():
        if set(shell) != set(magnetic_order(ell)):
            raise ValueError("incomplete AO radial shell")
    result = np.zeros((len(labels), len(labels)))
    for (atom, ell, zeta), shell in shells.items():
        target = shells.get((mapping[atom], ell, zeta))
        if target is None:
            raise ValueError("mapped atoms have different radial shells")
        order = magnetic_order(ell)
        result[np.ix_([target[m] for m in order], [shell[m] for m in order])] = real_harmonic_operation(r, ell, atol=atol)
    return AOOperation(result, np.array(mapping), np.array(shifts), np.array(errors), symprec)
