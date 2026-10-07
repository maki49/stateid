"""Unitary Seitz operations in a primitive cell; fractional column coordinates.

x' = W x + t; reciprocal fractional k' = W^{-T} k. No time reversal,
spin rotations, cell standardization, or symmetrization is inferred.
"""
from dataclasses import dataclass
import numpy as np
from stateid._validation import positive_tolerance


def real_array(value, shape, name):
    a = np.asarray(value)
    if a.shape != shape or not np.isrealobj(a) or not np.all(np.isfinite(a)):
        raise ValueError(f"{name} must be a finite real array of shape {shape}")
    return a.astype(float)


def fractional_rotation(value, *, atol=1e-8):
    positive_tolerance(atol)
    w = real_array(value, (3, 3), "fractional rotation")
    if not np.allclose(w, np.rint(w), atol=atol, rtol=0):
        raise ValueError("fractional rotation must be integer")
    w = np.rint(w).astype(int)
    if not np.isclose(abs(np.linalg.det(w)), 1, atol=atol, rtol=0):
        raise ValueError("fractional rotation must be unimodular")
    return w


def transform_kpoint(kpoint, rotation, *, atol=1e-8):
    """Return W^{-T} k without folding; W acts on direct fractional columns."""
    k = real_array(kpoint, (3,), "kpoint")
    return np.linalg.solve(fractional_rotation(rotation, atol=atol).T, k)


def _operations(rotations, translations, atol):
    values = np.asarray(rotations)
    if values.ndim != 3 or values.shape[1:] != (3, 3) or len(values) == 0:
        raise ValueError("rotations must be a nonempty (n,3,3) array")
    w = np.array([fractional_rotation(v, atol=atol) for v in values])
    t = real_array(translations, (len(w), 3), "translations")
    if len({tuple(v.flat) for v in w}) != len(w):
        raise ValueError("use primitive-cell coset representatives: repeated rotations "
                         "(centering translations or duplicate operations) are unsupported")
    return w, t


def _products(rotations, translations, atol):
    index = {tuple(w.flat): i for i, w in enumerate(rotations)}
    size = len(rotations)
    table = np.empty((size, size), dtype=int)
    shifts = np.empty((size, size, 3), dtype=int)
    for a, wa in enumerate(rotations):
        for b, wb in enumerate(rotations):
            c = index.get(tuple((wa @ wb).flat))
            if c is None:
                raise ValueError("operations are not closed under multiplication")
            delta = translations[a] + wa @ translations[b] - translations[c]
            if not np.allclose(delta, np.rint(delta), atol=atol, rtol=0):
                raise ValueError("Seitz products are not closed modulo lattice translations")
            table[a, b] = c
            shifts[a, b] = np.rint(delta).astype(int)
    return table, shifts


@dataclass(frozen=True)
class LittleGroup:
    indices: np.ndarray  # indices in the caller's full operation list
    rotations: np.ndarray
    translations: np.ndarray
    kpoint: np.ndarray
    reciprocal_shifts: np.ndarray  # W^{-T} k - k, integer for selected operations
    multiplication: np.ndarray
    lattice_shifts: np.ndarray  # g_a g_b = {I|L_ab} g_c
    factor_system: np.ndarray  # exp(-2 pi i k.L_ab)
    atol: float

    @property
    def order(self):
        return len(self.indices)


def find_little_group(rotations, translations, kpoint, *, atol=1e-8):
    """Select W^{-T}k = k + G and construct the finite projective group law.

    Input is a complete group modulo primitive translations, with one Seitz
    representative per rotation. Integer parts of translations are retained.
    Geometry and electronic invariance are separate checks.
    """
    positive_tolerance(atol)
    w, t = _operations(rotations, translations, atol)
    _products(w, t, atol)  # reject an incomplete input group, even at generic k
    k = real_array(kpoint, (3,), "kpoint")
    dk = np.array([np.linalg.solve(op.T, k) - k for op in w])
    selected = np.flatnonzero(np.max(abs(dk - np.rint(dk)), axis=1) <= atol)
    if not len(selected):
        raise ValueError("no identity operation in the little group")
    wl, tl = w[selected], t[selected]
    product, shifts = _products(wl, tl, atol)
    factor = np.exp(-2j * np.pi * (shifts @ k))
    return LittleGroup(selected, wl, tl, k.copy(), np.rint(dk[selected]).astype(int),
                       product, shifts, factor, atol)


def symmetry_from_structure(lattice, positions, numbers, *, symprec=1e-6):
    """Discover geometric operations with optional spglib; no electronic claim.

    Returned operations use the supplied cell. Centered/nonprimitive output
    must be converted explicitly before find_little_group; AO rows are never
    silently moved into a new cell.
    """
    positive_tolerance(symprec, "symprec")
    try:
        import spglib
    except ImportError as exc:
        raise ImportError("install stateid[symmetry] for symmetry discovery") from exc
    lattice = real_array(lattice, (3, 3), "lattice")
    p = np.asarray(positions)
    if p.ndim != 2 or p.shape[1:] != (3,) or not len(p):
        raise ValueError("positions must be a nonempty (nat,3) array")
    p = real_array(p, p.shape, "positions")
    n = real_array(numbers, (len(p),), "numbers")
    if np.any(n <= 0) or np.any(n != np.rint(n)) or np.linalg.matrix_rank(lattice) != 3:
        raise ValueError("require positive integer species numbers and nonsingular lattice")
    result = spglib.get_symmetry((lattice, p, n.astype(int)), symprec=symprec)
    if result is None:
        raise ValueError("spglib could not identify geometric symmetry")
    return result["rotations"], result["translations"]


def representation_errors(matrices, group):
    """Return (unitarity, projective multiplication) max normalized residuals."""
    d = np.asarray(matrices, dtype=complex)
    if (d.ndim != 3 or d.shape[0] != group.order or d.shape[1] == 0
            or d.shape[1] != d.shape[2] or not np.all(np.isfinite(d))):
        raise ValueError("matrices must be finite (little_group_order, dim, dim)")
    dim = d.shape[1]
    unitary = max(np.linalg.norm(x.conj().T @ x - np.eye(dim)) for x in d) / np.sqrt(dim)
    law = max(np.linalg.norm(d[a] @ d[b] - group.factor_system[a, b]
                             * d[group.multiplication[a, b]])
              for a in range(group.order) for b in range(group.order)) / np.sqrt(dim)
    return float(unitary), float(law)
