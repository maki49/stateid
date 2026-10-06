"""Versioned, pickle-free exchange format; NOT a native ABACUS format."""

import numpy as np

from stateid._validation import matrix, metric
from .models import LRData


def read_overlap_npy(path):
    """Read a dense S array saved with numpy.save; validate Hermitian positivity."""
    s = matrix(np.load(path, allow_pickle=False), "overlap")
    return metric(s, s.shape[0])


def read_lr_npz(path):
    """Read stateid schema 1 (see docs/io-formats.md).

    File must explicitly identify method and global ph ordering. No reshaping,
    normalization or merging of MPI-local files is attempted.
    """
    with np.load(path, allow_pickle=False) as data:
        required = {"schema_version", "method", "X", "energies_ev", "transitions"}
        if not required <= set(data.files):
            raise ValueError(f"missing LR keys: {sorted(required - set(data.files))}")
        if set(data.files) - required - {"Y"}:
            raise ValueError("unknown LR keys; schema version 1 uses explicit fields only")
        version = data["schema_version"]
        if version.shape != () or version.dtype.kind not in "iu" or version.item() != 1:
            raise ValueError("unsupported LR schema_version; expected integer scalar 1")
        method_array = data["method"]
        if method_array.shape != () or method_array.dtype.kind != "U":
            raise ValueError("method must be a scalar Unicode string")
        method = str(method_array.item())
        if method not in {"TDA", "full_lr"}:
            raise ValueError("method must be TDA or full_lr")
        x = matrix(data["X"], "X").copy()
        if 0 in x.shape:
            raise ValueError("X must contain at least one transition and root")
        e = np.asarray(data["energies_ev"])
        if e.dtype.kind not in "fiu" or e.shape != (x.shape[1],) or not np.all(np.isfinite(e)):
            raise ValueError("energies_ev must contain one finite real energy per root")
        tr = np.asarray(data["transitions"])
        if tr.shape != (x.shape[0], 5) or tr.dtype.kind not in "iu" or np.any(tr < 0):
            raise ValueError("transitions must be nonnegative integer rows (k,spin_h,i,spin_p,a)")
        if np.any(tr[:, [1, 3]] > 1) or len(np.unique(tr, axis=0)) != len(tr):
            raise ValueError("transition spin codes must be 0/1 and rows must be unique")
        if np.any((tr[:, 1] == tr[:, 3]) & (tr[:, 2] == tr[:, 4])):
            raise ValueError("a transition cannot excite an orbital into itself")
        y = matrix(data["Y"], "Y").copy() if "Y" in data.files else None
        if method == "TDA" and y is not None:
            raise ValueError("TDA schema must omit Y")
        if method == "full_lr" and (y is None or y.shape != x.shape):
            raise ValueError("full_lr requires Y with the same shape as X")
        return LRData(x, y, e.astype(float).copy(), tr.copy(), method, str(path))
