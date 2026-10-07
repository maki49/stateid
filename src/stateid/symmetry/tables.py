"""Class characters, decomposition, and conservative irrep matching."""

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np

from stateid._validation import positive_tolerance


@dataclass(frozen=True)
class CharacterTable:
    name: str
    classes: tuple[str, ...]
    class_sizes: tuple[int, ...]
    irreps: Mapping[str, tuple[complex, ...]]


C3V = CharacterTable(
    "C3v", ("E", "2C3", "3sigma_v"), (1, 2, 3),
    MappingProxyType({"A1": (1, 1, 1), "A2": (1, 1, -1), "E": (2, -1, 0)}),
)
C3V_CLASSES = MappingProxyType({
    "E": ("E",), "2C3": ("C3", "C3^2"), "3sigma_v": ("sv1", "sv2", "sv3"),
})


@dataclass(frozen=True)
class CharacterMatch:
    irrep: str | None
    multiplicities: dict[str, int]
    raw_multiplicities: dict[str, complex]
    residual: float
    valid: bool

    @property
    def label(self):
        if not self.valid:
            return "unmatched"
        return " + ".join(name if n == 1 else f"{n}{name}" for name, n in self.multiplicities.items())


def decompose_characters(characters, reference_characters, *, weights=None, atol=1e-5):
    """Decompose against named reference characters, without forcing a match.

    Default: one value per operation. Optional positive integer weights allow
    class representatives. All rows must describe the same ordered operations,
    same Seitz representatives, k and factor system. Character compatibility
    alone does not certify matrix group laws or closure.
    """
    positive_tolerance(atol)
    chi = np.asarray(characters, dtype=complex)
    if chi.ndim != 1 or not len(chi) or not np.all(np.isfinite(chi)):
        raise ValueError("characters must be a finite nonempty vector")
    names = list(reference_characters)
    if not names or any(not isinstance(n, str) or not n for n in names):
        raise ValueError("reference characters need nonempty string identifiers")
    rows = np.asarray(list(reference_characters.values()), dtype=complex)
    if rows.shape != (len(names), len(chi)) or not np.all(np.isfinite(rows)):
        raise ValueError("reference character rows must match the operation count")
    w = np.ones(len(chi)) if weights is None else np.asarray(weights)
    if (w.shape != chi.shape or not np.isrealobj(w) or not np.all(np.isfinite(w))
            or np.any(w <= 0) or np.any(w != np.rint(w))):
        raise ValueError("weights must be positive integer class sizes")
    w = w / w.sum()
    if not np.allclose((rows.conj() * w) @ rows.T, np.eye(len(rows)),
                       atol=atol, rtol=0):
        raise ValueError("reference characters must be mutually orthonormal")
    raw = rows.conj() @ (w * chi)
    rounded = np.rint(raw.real)
    residual = float(np.max(np.abs(chi - rounded @ rows)))
    valid = bool(np.all(np.abs(raw - rounded) <= atol)
                 and np.all(rounded >= 0) and np.sum(rounded) > 0 and residual <= atol)
    counts = {name: int(n) for name, n in zip(names, rounded) if n > 0} if valid else {}
    irrep = next(iter(counts)) if valid and sum(counts.values()) == 1 else None
    return CharacterMatch(irrep, counts, dict(zip(names, map(complex, raw))), residual, valid)


def match_characters(characters, table=C3V, *, atol=1e-5):
    """Match one character per class (not class sums), using class-size weights."""
    return decompose_characters(characters, table.irreps,
                                weights=table.class_sizes, atol=atol)


def class_characters(operation_characters, classes=C3V_CLASSES, *, atol=1e-5):
    """Average all group elements within classes, rejecting class inconsistency."""
    positive_tolerance(atol)
    required = {name for members in classes.values() for name in members}
    if set(operation_characters) != required:
        raise ValueError(f"expected exactly these operations: {sorted(required)}")
    means = []
    for members in classes.values():
        values = np.array([operation_characters[name] for name in members], dtype=complex)
        if not np.all(np.isfinite(values)) or np.max(np.abs(values - values.mean())) > atol:
            raise ValueError("characters are inconsistent within a conjugacy class")
        means.append(values.mean())
    return np.asarray(means)
