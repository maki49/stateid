"""A checked C3v workflow over all six unitary operations."""

from dataclasses import dataclass

import numpy as np

from .representation import Representation, project_representation
from .tables import C3V_CLASSES, CharacterMatch, class_characters, match_characters


@dataclass(frozen=True)
class SymmetryAnalysis:
    representations: dict[str, Representation]
    characters: np.ndarray
    match: CharacterMatch
    group_error: float


def analyze_c3v(coefficients, operations, overlap=None, *, operator_kind, atol=1e-5):
    """Require closure, metric preservation and C3v relations before matching.

    Operation names: E, C3, C3^2, sv1, sv2=C3*sv1, sv3=C3^2*sv1.
    Products use active column-vector operators (rightmost acts first).
    """
    required = {name for members in C3V_CLASSES.values() for name in members}
    if set(operations) != required:
        raise ValueError(f"expected exactly these operations: {sorted(required)}")
    reps = {name: project_representation(coefficients, op, overlap,
            operator_kind=operator_kind, atol=atol) for name, op in operations.items()}
    for name, rep in reps.items():
        if rep.ao_metric_error > atol:
            raise ValueError(f"{name}: AO operation does not preserve the overlap metric")
        if rep.closure_error > atol:
            raise ValueError(f"{name}: subspace is not closed (error={rep.closure_error:.3g})")
    d = {name: rep.matrix for name, rep in reps.items()}
    r, f = d["C3"], d["sv1"]
    eye = np.eye(r.shape[0])
    differences = [d["E"] - eye, r @ r @ r - eye, f @ f - eye,
                   f @ r @ f - r @ r, d["C3^2"] - r @ r,
                   d["sv2"] - r @ f, d["sv3"] - r @ r @ f]
    group_error = float(max(np.linalg.norm(delta) for delta in differences))
    if group_error > atol:
        raise ValueError(f"projected operations violate C3v group relations ({group_error:.3g})")
    chi = class_characters({name: rep.character for name, rep in reps.items()}, atol=atol)
    return SymmetryAnalysis(reps, chi, match_characters(chi, atol=atol), group_error)
