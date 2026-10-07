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
    result = project_sewing_matrix(
        coefficients, coefficients, operation, overlap, overlap,
        operator_kind=operator_kind, orthonormalize_basis=orthonormalize_basis, atol=atol)
    return Representation(result.matrix, complex(np.trace(result.matrix)),
                          result.closure_error, result.ao_metric_error)


@dataclass(frozen=True)
class SewingMatrix:
    """A map between independently gauged subspaces; no invariant character."""
    matrix: np.ndarray
    closure_error: float
    ao_metric_error: float


def project_sewing_matrix(source_coefficients, target_coefficients, operation,
                          source_overlap=None, target_overlap=None, *,
                          operator_kind="coefficient", orthonormalize_basis=True, atol=1e-8):
    """D_(k' <- k)=Q_target† S_target B(g,k) Q_source.

    matrix_element means M=S_target B. Checks B†S_target B=S_source in
    orthonormal AO coordinates and closure in the target metric. It does not
    identify k points, infer an irrep, or expose Tr(D) as a character.
    """
    positive_tolerance(atol)
    c = matrix(source_coefficients, "source_coefficients")
    same_basis = (source_coefficients is target_coefficients
                  and source_overlap is target_overlap)
    ct = c if same_basis else matrix(target_coefficients, "target_coefficients")
    if c.shape[1] == 0 or ct.shape[1] != c.shape[1]:
        raise ValueError("source/target subspaces must have the same nonzero dimension")
    s = metric(source_overlap, c.shape[0])
    st = s if same_basis else metric(target_overlap, ct.shape[0])
    op = matrix(operation, "operation")
    if op.shape != (ct.shape[0], c.shape[0]):
        raise ValueError("operation shape must match target/source AO dimensions")
    if operator_kind == "coefficient":
        b, m = op, st @ op
    elif operator_kind == "matrix_element":
        m, b = op, np.linalg.solve(st, op)
    else:
        raise ValueError("operator_kind must be 'coefficient' or 'matrix_element'")
    if orthonormalize_basis:
        q = orthonormalize(c, s)
        qt = q if same_basis else orthonormalize(ct, st)
    else:
        q, qt = c, ct
        for basis, overlap in ((q, s), (qt, st)):
            if not np.allclose(basis.conj().T @ overlap @ basis,
                               np.eye(basis.shape[1]), atol=atol, rtol=0):
                raise ValueError("coefficients are not S-orthonormal")
    d = qt.conj().T @ m @ q
    residual = b @ q - qt @ d
    closure = np.sqrt(max(0.0, np.trace(residual.conj().T @ st @ residual).real)
                      / q.shape[1])
    a = np.linalg.cholesky(s).conj().T
    at = np.linalg.cholesky(st).conj().T
    b_orth = np.linalg.solve(a.T, (at @ b).T).T
    error = np.linalg.norm(b_orth.conj().T @ b_orth - np.eye(s.shape[0])) / np.sqrt(s.shape[0])
    return SewingMatrix(d, float(closure), float(error))
