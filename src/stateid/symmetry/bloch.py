"""Bloch AO coefficient maps in the cell gauge, without an AO-position phase."""
from dataclasses import dataclass
import numpy as np
from stateid._validation import matrix, positive_tolerance
from .ao import AOOperation, build_gamma_ao_operation
from .groups import fractional_rotation, real_array, transform_kpoint


@dataclass(frozen=True)
class BlochAOOperation:
    matrix: np.ndarray
    kpoint: np.ndarray
    mapped_kpoint: np.ndarray
    ao_operation: AOOperation


def bloch_from_ao_operation(ao_operation, atom_indices, kpoint, rotation_fractional,
                          *, atol=1e-8):
    """Add exp(-2 pi i k'.L_source) to each source AO column of Gamma T.

    atom_indices maps coefficient rows to zero-based atoms. Cache the expensive
    AOOperation and reuse this function at many k points. The supplied W must
    be the operation used to build that AOOperation.
    """
    positive_tolerance(atol)
    t = matrix(ao_operation.matrix, "AO operation")
    if t.shape[0] != t.shape[1] or not len(t):
        raise ValueError("AO operation must be nonempty and square")
    shifts = np.asarray(ao_operation.lattice_shifts)
    if shifts.ndim != 2 or shifts.shape[1:] != (3,):
        raise ValueError("lattice_shifts must be (nat,3)")
    shifts = real_array(shifts, shifts.shape, "lattice_shifts")
    if not np.allclose(shifts, np.rint(shifts), atol=atol, rtol=0):
        raise ValueError("lattice_shifts must be integer")
    atoms = real_array(atom_indices, (len(t),), "atom_indices")
    if np.any(atoms != np.rint(atoms)) or np.any(atoms < 0) or np.any(atoms >= len(shifts)):
        raise ValueError("atom_indices must index the lattice-shift rows")
    k = real_array(kpoint, (3,), "kpoint")
    kp = transform_kpoint(k, rotation_fractional, atol=atol)
    phase = np.exp(-2j * np.pi * (shifts[atoms.astype(int)] @ kp))
    return BlochAOOperation(t * phase[None, :], k, kp, ao_operation)


def build_bloch_ao_operation(fractional_positions, lattice_vectors, species, labels,
                             rotation_fractional, kpoint, *,
                             translation_fractional=None, symprec=1e-6, atol=1e-8):
    """Build B(g,k) for scalar s/p/d AO columns, using fractional Seitz inputs.

    Unlike build_gamma_ao_operation, W and t here are FRACTIONAL. Lattice
    vectors are rows: R_cart=A.T W A^{-T}, t_cart=t @ A.
    The basis is |mu,k> = sum_R exp(+2 pi i k.R)|mu,R>.
    W tau_source+t=tau_target+L gives B_target,source =
    exp(-2 pi i (W^{-T}k).L) T_target,source.
    """
    labels = tuple(labels)
    w = fractional_rotation(rotation_fractional, atol=atol)
    a = real_array(lattice_vectors, (3, 3), "lattice_vectors")
    if np.linalg.matrix_rank(a) != 3:
        raise ValueError("lattice must be nonsingular")
    t = np.zeros(3) if translation_fractional is None else real_array(
        translation_fractional, (3,), "translation_fractional")
    cart = a.T @ w @ np.linalg.inv(a.T)
    ao = build_gamma_ao_operation(fractional_positions, a, species, labels, cart,
                                  translation=t @ a, symprec=symprec, atol=atol)
    return bloch_from_ao_operation(ao, [label.atom_index for label in labels],
                                   kpoint, w, atol=atol)
