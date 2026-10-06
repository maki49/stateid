"""Shared dense-array validation; no silent repairs of physical input."""

import numpy as np


def positive_tolerance(value, name="atol"):
    if not np.isscalar(value) or not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return float(value)


def matrix(value, name):
    result = np.asarray(value, dtype=complex)
    if result.ndim != 2 or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be a finite 2-D numeric array")
    return result


def metric(value, size, atol=1e-10):
    positive_tolerance(atol)
    result = np.eye(size, dtype=complex) if value is None else matrix(value, "overlap")
    if result.shape != (size, size) or size < 1:
        raise ValueError("overlap must be square and match the AO dimension")
    if not np.allclose(result, result.conj().T, atol=atol, rtol=0):
        raise ValueError("overlap must be Hermitian")
    result = (result + result.conj().T) / 2
    eig = np.linalg.eigvalsh(result)
    if eig[0] <= atol * max(1.0, eig[-1]):
        raise ValueError("overlap is singular, indefinite, or ill-conditioned; reduce the AO basis explicitly")
    return result
