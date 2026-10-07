"""Optional spgrep backend and checked scalar little-group projection."""
from dataclasses import dataclass
from importlib.metadata import version
import numpy as np
from stateid._validation import positive_tolerance
from .groups import LittleGroup, find_little_group, representation_errors
from .representation import Representation, project_representation
from .tables import CharacterMatch, decompose_characters


@dataclass(frozen=True)
class IrrepSet:
    group: LittleGroup
    matrices: tuple[np.ndarray, ...]
    backend: str

    @property
    def characters(self):
        # IDs are local to this result, NOT conventional irrep labels.
        return {f"irrep_{i}": np.trace(d, axis1=1, axis2=2)
                for i, d in enumerate(self.matrices)}


def spgrep_irreps(rotations, translations, kpoint, *, atol=1e-8):
    """Generate complex scalar small irreps, aligned to input little-group order.

    Primitive cell, unitary spatial operations only. Returned matrices already
    include nonsymmorphic phases; do not multiply a second translation phase.
    No standard labels, SOC, magnetic co-reps or automatic cell changes.
    """
    positive_tolerance(atol)
    group = find_little_group(rotations, translations, kpoint, atol=atol)
    try:
        import spgrep
    except ImportError as exc:
        raise ImportError("install stateid[symmetry] to generate irreps with spgrep") from exc
    # spgrep requires its identity representative to have zero translation.
    # Reduce only at the backend boundary, then restore the caller's Seitz
    # representatives by their lattice-translation phases.
    integer_shifts = np.floor(group.translations + atol).astype(int)
    canonical_t = group.translations - integer_shifts
    canonical_t[np.abs(canonical_t) < atol] = 0
    reps, mapping = spgrep.get_spacegroup_irreps_from_primitive_symmetry(
        group.rotations, canonical_t, group.kpoint,
        real=False, rtol=0, atol=atol)
    mapping = np.asarray(mapping)
    if sorted(mapping.tolist()) != list(range(group.order)):
        raise RuntimeError("spgrep and stateid disagree about little-group membership")
    order = np.argsort(mapping)
    phases = np.exp(-2j * np.pi * (integer_shifts @ group.kpoint))
    reps = tuple(np.asarray(rep, dtype=complex)[order] * phases[:, None, None]
                 for rep in reps)
    if sum(d.shape[1] ** 2 for d in reps) != group.order:
        raise RuntimeError("spgrep returned an incomplete set of irreducible representations")
    for d in reps:
        if max(representation_errors(d, group)) > 10 * atol:
            raise RuntimeError("spgrep representations violate the chosen Seitz phase convention")
    chars = np.array([np.trace(d, axis1=1, axis2=2) for d in reps])
    if not np.allclose(chars.conj() @ chars.T / group.order,
                       np.eye(len(reps)), atol=10*atol, rtol=0):
        raise RuntimeError("spgrep irreducible characters are not orthonormal")
    return IrrepSet(group, reps, f"spgrep {version('spgrep')}")


@dataclass(frozen=True)
class LittleGroupAnalysis:
    representations: tuple[Representation, ...]
    characters: np.ndarray
    match: CharacterMatch  # run-local IDs; assign conventional labels separately
    group_error: float
    unitarity_error: float
    reference: IrrepSet


def analyze_little_group(coefficients, operations, reference, overlap=None, *,
                         operator_kind="coefficient", atol=1e-5):
    """Project and validate before decomposing against an IrrepSet.

    operations contains ONLY the little-group operations in reference.group
    order. Supply B(g,k), not Gamma T. Matrices may be raw arrays or
    BlochAOOperation objects (the latter require operator_kind='coefficient').
    """
    from .bloch import BlochAOOperation
    positive_tolerance(atol)
    ops = tuple(operations)
    if len(ops) != reference.group.order:
        raise ValueError("supply one AO operation per reference little-group operation")
    reps = []
    for i, op in enumerate(ops):
        if isinstance(op, BlochAOOperation):
            if operator_kind != "coefficient":
                raise ValueError("BlochAOOperation contains coefficients B, not S B")
            k = reference.group.kpoint
            kp = np.linalg.solve(reference.group.rotations[i].T, k)
            if (not np.allclose(op.kpoint, k, atol=atol, rtol=0)
                    or not np.allclose(op.mapped_kpoint, kp, atol=atol, rtol=0)):
                raise ValueError("Bloch operation k points do not match the little group")
            op = op.matrix
        rep = project_representation(coefficients, op, overlap,
                                     operator_kind=operator_kind, atol=atol)
        if rep.ao_metric_error > atol:
            raise ValueError(f"operation {i}: AO map does not preserve the overlap metric")
        if rep.closure_error > atol:
            raise ValueError(f"operation {i}: selected subspace is not closed")
        reps.append(rep)
    unitary, law = representation_errors([r.matrix for r in reps], reference.group)
    if max(unitary, law) > atol:
        raise ValueError("projected operations violate little-group relations or unitarity")
    chi = np.array([r.character for r in reps])
    match = decompose_characters(chi, reference.characters, atol=atol)
    return LittleGroupAnalysis(tuple(reps), chi, match, law, unitary, reference)
