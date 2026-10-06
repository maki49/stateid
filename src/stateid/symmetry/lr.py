"""A first TDA symmetry building block for a complete particle-hole product."""

import numpy as np

from stateid._validation import matrix, positive_tolerance


def particle_hole_operation(occupied_operation, virtual_operation, *, reference_character, atol=1e-8):
    """Return chi_ref(R) * (D_occ(R)* tensor D_virt(R)).

    Basis: |Phi_i^a>=a_a†a_i|Phi_ref>, i outer / a inner (virtual index fast).
    Requires invariant, orthonormal occupied/virtual spaces in ONE spin block,
    a one-dimensional symmetry-pure reference, and the full Cartesian ph set.
    Spin-adapted or truncated ABACUS LR layouts require a separate adapter.
    reference_character=1 returns the excitation-operator transformation when
    the physical reference phase is intentionally omitted; do not then label
    the total many-electron symmetry as if that phase had been included.
    """
    positive_tolerance(atol)
    occ, vir = matrix(occupied_operation, "occupied_operation"), matrix(virtual_operation, "virtual_operation")
    for op in (occ, vir):
        if not op.size or op.shape[0] != op.shape[1] or not np.allclose(op.conj().T @ op, np.eye(op.shape[0]), atol=atol, rtol=0):
            raise ValueError("occupied and virtual operations must be nonempty square unitary matrices")
    phase = complex(reference_character)
    if not np.isfinite(phase) or abs(abs(phase) - 1) > atol:
        raise ValueError("reference_character must be a unit-modulus scalar for a 1-D reference")
    return phase * np.kron(occ.conj(), vir)
