"""ABACUS real-harmonic basis conversion, independent of Euler-angle choices.

This is a mathematical basis-change helper, not an ABACUS AO rotation builder.
See docs/abacus-conventions.md for active/passive and storage-layout caveats.
"""

import operator

import numpy as np

from stateid._validation import matrix


def magnetic_order(ell):
    """Return m=(0,+1,-1,...,+ell,-ell), the local ABACUS shell order."""
    if isinstance(ell, bool):
        raise ValueError("ell must be a nonnegative integer")
    try:
        ell = operator.index(ell)
    except TypeError as exc:
        raise ValueError("ell must be a nonnegative integer") from exc
    if ell < 0:
        raise ValueError("ell must be a nonnegative integer")
    return (0,) + tuple(m for n in range(1, ell + 1) for m in (n, -n))


def complex_to_real_basis(ell):
    """Return B_mm'=<Y_l^m|S_l^m'>; both axes use magnetic_order(ell).

    For m>0: S_m=(Y_m+(-1)^m Y_-m)/sqrt(2),
    S_-m=(-i Y_m+i(-1)^m Y_-m)/sqrt(2).
    """
    order = magnetic_order(ell)
    b = np.zeros((len(order), len(order)), dtype=complex)
    b[0, 0] = 1
    for m in range(1, ell + 1):
        p, n = 2 * m - 1, 2 * m
        b[p, p], b[n, p] = 1 / np.sqrt(2), (-1)**m / np.sqrt(2)
        b[p, n], b[n, n] = -1j / np.sqrt(2), 1j * (-1)**m / np.sqrt(2)
    return b


def rotation_in_real_basis(complex_rotation, ell):
    """Change basis B†D_Y B, preserving the supplied operation convention."""
    b = complex_to_real_basis(ell)
    d = matrix(complex_rotation, "complex_rotation")
    if d.shape != b.shape:
        raise ValueError("rotation dimension does not match ell")
    return b.conj().T @ d @ b
