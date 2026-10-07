"""A first TDA symmetry building block for a complete particle-hole product."""

from dataclasses import dataclass

import numpy as np

from stateid._validation import matrix, positive_tolerance
from .analysis import SymmetryAnalysis, analyze_c3v, _match_c3v_representations
from .representation import Representation
from .tables import C3V_CLASSES


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


@dataclass(frozen=True)
class TDASymmetryAnalysis:
    symmetry: SymmetryAnalysis
    reference_characters: dict[str, complex]
    root_orthonormality_error: float


def _checked_orbital_blocks(operations, atol):
    required = {name for members in C3V_CLASSES.values() for name in members}
    if set(operations) != required:
        raise ValueError(f"expected exactly these operations: {sorted(required)}")
    count = len(operations['E'])
    if not count or any(len(blocks) != count for blocks in operations.values()):
        raise ValueError("operation block counts must be nonzero and consistent")
    result = []
    for block in range(count):
        ops = {name: matrix(values[block], 'orbital operation') for name, values in operations.items()}
        size = ops['E'].shape[0]
        if any(op.shape != (size, size) for op in ops.values()):
            raise ValueError("orbital operation dimensions must be consistent and square")
        if size == 0:
            # Empty beta space of a fully polarized reference is legitimate.
            result.append(ops)
            continue
        analyze_c3v(np.eye(size), ops, operator_kind='coefficient', atol=atol)
        result.append(ops)
    return result


def analyze_tda_c3v(amplitudes, operations, *, reference_occupied_operations, atol=1e-5):
    """Analyze total spatial irrep of orthonormal Gamma TDA roots without kron.

    operations[name] is a sequence of (D_occ, D_virt) pairs, one per spin
    excitation block. Rows of amplitudes concatenate those complete Cartesian
    ph products, occupied outer / virtual fast. No sparse/truncated ph layout,
    multi-k, spin-flip or full LR X/Y semantics are supported.

    reference_occupied_operations[name] contains TWO matrices for the COMPLETE
    occupied alpha/beta spaces of an integer-occupied Slater reference. These
    can differ from the excitation occupied windows; their determinant product
    supplies the single common reference phase. Matrices must come from closed,
    S-orthonormal orbital spaces. Fractional occupations or a non-determinant
    reference cannot be represented by this contract.
    """
    atol = positive_tolerance(atol)
    x = matrix(amplitudes, 'amplitudes')
    if not x.shape[0] or not x.shape[1]:
        raise ValueError("select at least one TDA root and one particle-hole row")
    required = {name for members in C3V_CLASSES.values() for name in members}
    if set(operations) != required:
        raise ValueError(f"expected exactly these operations: {sorted(required)}")
    if any(not isinstance(pair, (tuple, list)) or len(pair) != 2
           for blocks in operations.values() for pair in blocks):
        raise ValueError("each excitation block must specify occupied and virtual operations")
    if any(len(blocks) != 2 for blocks in reference_occupied_operations.values()):
        raise ValueError("reference must specify complete alpha and beta occupied operations")
    reference = _checked_orbital_blocks(reference_occupied_operations, atol)
    occ = _checked_orbital_blocks({name: [pair[0] for pair in blocks]
                                   for name, blocks in operations.items()}, atol)
    vir = _checked_orbital_blocks({name: [pair[1] for pair in blocks]
                                   for name, blocks in operations.items()}, atol)
    if len(occ) not in (1, 2) or len(vir) != len(occ):
        raise ValueError("expected one or two spin-conserving excitation blocks")
    dimensions = [(o['E'].shape[0], v['E'].shape[0]) for o,v in zip(occ,vir)]
    if any(no == 0 or nv == 0 for no,nv in dimensions):
        raise ValueError("excitation blocks must be nonempty")
    if sum(no*nv for no,nv in dimensions) != x.shape[0]:
        raise ValueError("amplitude rows do not match the complete particle-hole block layout")
    orth_error = float(np.linalg.norm(x.conj().T @ x-np.eye(x.shape[1])))
    if orth_error > atol:
        raise ValueError(f"TDA roots are not orthonormal (error={orth_error:.3g})")
    phases, reps = {}, {}
    for name in required:
        phase = complex(np.prod([np.linalg.det(block[name]) for block in reference]))
        if abs(abs(phase)-1) > atol:
            raise ValueError("reference determinant character does not have unit modulus")
        phases[name] = phase
        transformed = np.empty_like(x)
        offset = 0
        for o,v,(no,nv) in zip(occ,vir,dimensions):
            end = offset+no*nv
            tensor = x[offset:end].reshape(no,nv,x.shape[1])
            # For each root: D_occ.conj() @ X @ D_virt.T.
            transformed[offset:end] = (phase*np.einsum('ij,jbs,ab->ias',
                o[name].conj(), tensor, v[name], optimize=True)).reshape(no*nv,x.shape[1])
            offset = end
        d = x.conj().T @ transformed
        closure = float(np.linalg.norm(transformed-x @ d)/np.sqrt(x.shape[1]))
        if closure > atol:
            raise ValueError(f"{name}: TDA root subspace is not closed (error={closure:.3g})")
        # Orbital block unitarity was checked before contraction; no ph metric
        # matrix is materialized. This diagnostic bounds all product blocks.
        error = max(np.linalg.norm(block[name].conj().T @ block[name]-np.eye(block[name].shape[0]))
                    for block in occ+vir)
        reps[name] = Representation(d, complex(np.trace(d)), closure, float(error))
    return TDASymmetryAnalysis(_match_c3v_representations(reps, atol=atol), phases, orth_error)
