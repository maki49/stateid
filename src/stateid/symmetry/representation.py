"""Represent unitary spatial operations in an orbital subspace."""

from dataclasses import dataclass
from typing import Literal

import numpy as np

from stateid._validation import matrix, metric, positive_tolerance


@dataclass(frozen=True)
class Representation:
    matrix: np.ndarray
    character: complex
    closure_error: float
    ao_metric_error: float


def orthonormalize(coefficients, overlap=None, *, rank_tol=1e-10):
    """Return Q=C(C†SC)^(-1/2), without silently dropping dependent columns.

    Columns are orbital coefficients. A singular AO metric or selected subspace
    raises ValueError; explicit basis truncation belongs to the caller.
    """
    positive_tolerance(rank_tol, "rank_tol")
    c = matrix(coefficients, "coefficients")
    if c.shape[1] == 0:
        raise ValueError("select at least one orbital")
    s = metric(overlap, c.shape[0], rank_tol)
    g = c.conj().T @ s @ c
    eig, u = np.linalg.eigh((g + g.conj().T) / 2)
    if eig[0] <= rank_tol * max(1.0, eig[-1]):
        raise ValueError("selected orbitals are linearly dependent or ill-conditioned")
    return c @ ((u * (1 / np.sqrt(eig))) @ u.conj().T)


def project_representation(
    coefficients,
    operation,
    overlap=None,
    *,
    operator_kind: Literal["coefficient", "matrix_element"],
    orthonormalize_basis=True,
    atol=1e-8,
):
    """Compute D=Q†M Q, chi=Tr(D), and subspace-closure diagnostics.

    coefficient: operation=T, R|phi_nu>=sum_mu |phi_mu>T_mu,nu, M=S T.
    matrix_element: operation=M, M_mu,nu=<phi_mu|R|phi_nu>, T=solve(S,M).
    Q is S-orthonormal. Antiunitary operations and k-changing maps are excluded.
    closure_error = ||TQ-QD||_S / sqrt(number of selected orbitals).
    """
    positive_tolerance(atol)
    c = matrix(coefficients, "coefficients")
    if c.shape[1] == 0:
        raise ValueError("select at least one orbital")
    s = metric(overlap, c.shape[0])
    op = matrix(operation, "operation")
    if op.shape != s.shape:
        raise ValueError("operation and overlap shapes must match")
    if operator_kind == "coefficient":
        t, m = op, s @ op
    elif operator_kind == "matrix_element":
        m, t = op, np.linalg.solve(s, op)
    else:
        raise ValueError("operator_kind must be 'coefficient' or 'matrix_element'")
    if orthonormalize_basis:
        q = orthonormalize(c, s)
    else:
        q = c
        if not np.allclose(q.conj().T @ s @ q, np.eye(q.shape[1]), atol=atol, rtol=0):
            raise ValueError("coefficients are not S-orthonormal")
    d = q.conj().T @ m @ q
    residual = t @ q - q @ d
    closure = np.sqrt(max(0.0, np.trace(residual.conj().T @ s @ residual).real) / q.shape[1])
    # Measure unitarity in orthonormal AO coordinates, independent of AO scaling.
    l = np.linalg.cholesky(s)
    a = l.conj().T
    t_orth = np.linalg.solve(a.T, (a @ t).T).T
    error = np.linalg.norm(t_orth.conj().T @ t_orth - np.eye(s.shape[0])) / np.sqrt(s.shape[0])
    return Representation(d, complex(np.trace(d)), float(closure), float(error))
