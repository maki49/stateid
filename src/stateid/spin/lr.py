"""Explicit boundaries between determinant, TDA, and full-LR spin analysis."""

import numpy as np

from stateid._validation import matrix, positive_tolerance


def tda_s2(x, s2_ph, *, atol=1e-8):
    """Contract X† S²_ph X/(X†X) in an ORTHONORMAL determinant/CSF basis.

    Caller supplies the full absolute S² operator (including reference
    contributions), with exactly the same ph ordering and determinant phases.
    This routine does not construct that operator or interpret raw TDDFT X.
    One root only; units hbar². Scaling X leaves the result unchanged.
    """
    positive_tolerance(atol)
    x = np.asarray(x, dtype=complex)
    if x.ndim != 1 or not x.size or not np.all(np.isfinite(x)):
        raise ValueError("x must be a finite, nonempty vector for one TDA root")
    op = matrix(s2_ph, "s2_ph")
    if op.shape != (x.size, x.size):
        raise ValueError("s2_ph and x dimensions do not match")
    if not np.allclose(op, op.conj().T, atol=atol, rtol=0):
        raise ValueError("s2_ph must be Hermitian")
    norm = np.vdot(x, x).real
    if norm <= 0:
        raise ValueError("x has zero norm")
    return float((np.vdot(x, op @ x) / norm).real)


def build_s2_ph(reference, transitions):
    """Reserved: construct S² between spin-resolved excited determinants."""
    raise NotImplementedError("S²_ph construction needs occupied alpha/beta overlaps, ph order and fermion phases")


def lr_s2(x, y, reference, transitions):
    """Reserved: full non-TDA LR requires a documented response-state prescription."""
    raise NotImplementedError("Full LR <S²> cannot be obtained by treating (X,Y) as a normalized CI vector")
