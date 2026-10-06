"""Spin diagnostics for collinear, integer-occupied Slater determinants."""

import numpy as np

from stateid._validation import matrix, metric, positive_tolerance


def determinant_s2(occupied_alpha, occupied_beta, overlap=None, *, atol=1e-8):
    """Return <S²>/hbar² = Ms² + (Na+Nb)/2 - ||Ca† S Cb||_F².

    Each column is one occupied spatial orbital (one electron of that spin).
    Both spin blocks must be individually S-orthonormal in the SAME AO basis.
    Empty spin blocks use shape (nao,0). Fractional occupations, spinors,
    k-weighted ensembles and excited-state LR amplitudes are not accepted here.
    UKS use diagnoses its auxiliary determinant, not exact interacting <S²>.
    """
    positive_tolerance(atol)
    ca, cb = matrix(occupied_alpha, "occupied_alpha"), matrix(occupied_beta, "occupied_beta")
    if ca.shape[0] != cb.shape[0]:
        raise ValueError("alpha and beta orbitals must use the same AO basis")
    s = metric(overlap, ca.shape[0])
    for c in (ca, cb):
        if not np.allclose(c.conj().T @ s @ c, np.eye(c.shape[1]), atol=atol, rtol=0):
            raise ValueError("occupied orbitals must be S-orthonormal within each spin")
    na, nb = ca.shape[1], cb.shape[1]
    ms = (na - nb) / 2
    return float(ms**2 + (na + nb) / 2 - np.linalg.norm(ca.conj().T @ s @ cb)**2)
