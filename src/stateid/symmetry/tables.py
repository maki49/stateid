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


def match_characters(characters, table=C3V, *, atol=1e-5):
    """Match one character PER CLASS (not summed over class members).

    Decompose using class-size weights, then require nonnegative integral
    multiplicities AND reconstruction within atol. Never force a nearest label.
    A valid reducible representation has irrep=None and nonempty multiplicities.
    This function alone cannot certify group laws or subspace closure.
    """
    positive_tolerance(atol)
    chi = np.asarray(characters, dtype=complex)
    if chi.shape != (len(table.classes),) or not np.all(np.isfinite(chi)):
        raise ValueError("supply a finite character for each class in table order")
    rows = np.asarray(list(table.irreps.values()), dtype=complex)
    raw = rows.conj() @ (np.asarray(table.class_sizes) * chi) / sum(table.class_sizes)
    rounded = np.rint(raw.real)
    reconstructed = rounded @ rows
    residual = float(np.max(np.abs(chi - reconstructed)))
    valid = bool(
        np.all(np.abs(raw - rounded) <= atol)
        and np.all(rounded >= 0) and np.sum(rounded) > 0 and residual <= atol
    )
    names = list(table.irreps)
    counts = {name: int(n) for name, n in zip(names, rounded) if n > 0} if valid else {}
    irrep = next(iter(counts)) if valid and sum(counts.values()) == 1 else None
    return CharacterMatch(irrep, counts, dict(zip(names, map(complex, raw))), residual, valid)


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
