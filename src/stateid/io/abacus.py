"""Small, strict ABACUS adapter. Unsupported readers fail explicitly."""

from pathlib import Path
import re

import numpy as np

from .arrays import read_lr_npz, read_overlap_npy
from .models import OrbitalData


class UnsupportedFormatError(NotImplementedError):
    """The input requires an adapter that has not yet been implemented."""


def read_gamma_wavefunctions(path, *, spin="unspecified"):
    """Read one real Gamma LCAO text frame from wfc_nao_write2file.

    Header: 'N (number of bands)', 'M (number of orbitals)'. Each band has
    '(band)', '(Ry)', '(Occupations)', then exactly M real coefficients.
    Multi-k, binary, spinors, append-mode trajectories and MPI shards are outside
    this reader. Spin is supplied by the caller, not guessed from a filename.
    """
    if spin not in {"alpha", "beta", "unspecified"}:
        raise ValueError("spin must be alpha, beta, or unspecified")
    try:
        lines = [line.strip() for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip()]
    except UnicodeError as exc:
        raise UnsupportedFormatError("binary wavefunctions are not supported") from exc
    if not lines or not re.fullmatch(r"\d+\s+\(number of bands\)", lines[0]):
        raise UnsupportedFormatError("expected real Gamma LCAO text header; multi-k/binary formats are not supported")
    pos = 0

    def field(label, integer=False):
        nonlocal pos
        if pos >= len(lines):
            raise ValueError(f"truncated file: expected ({label})")
        match = re.fullmatch(r"(\S+)\s+\(" + re.escape(label) + r"\)", lines[pos])
        if match is None:
            raise ValueError(f"line {pos + 1}: expected ({label})")
        pos += 1
        token = match[1]
        value = int(token) if integer else float(token.replace("D", "E").replace("d", "e"))
        if not np.isfinite(value):
            raise ValueError("non-finite value in wavefunction file")
        return value

    nb = field("number of bands", True)
    nao = field("number of orbitals", True)
    if min(nb, nao) < 1:
        raise ValueError("band and AO dimensions must be positive")
    # Accumulate only available data; do not allocate from an unverified header.
    columns, energies, occupations = [], [], []
    for band in range(1, nb + 1):
        if field("band", True) != band:
            raise ValueError("band labels must be consecutive and one-based")
        energies.append(field("Ry"))
        occupations.append(field("Occupations"))
        column = []
        while len(column) < nao:
            if pos >= len(lines):
                raise ValueError("truncated AO coefficients")
            try:
                column.extend(float(t.replace("D", "E").replace("d", "e")) for t in lines[pos].split())
            except ValueError as exc:
                raise ValueError(f"line {pos + 1}: expected real AO coefficients") from exc
            pos += 1
            if len(column) > nao or not np.all(np.isfinite(column)):
                raise ValueError("wrong coefficient count or non-finite coefficient")
        columns.append(column)
    if pos != len(lines):
        raise UnsupportedFormatError("extra data after one frame; append-mode/multiple frames are not supported")
    return OrbitalData(np.array(columns).T, np.array(energies), np.array(occupations),
                       np.arange(1, nb + 1), spin, str(path))


class AbacusReader:
    """ABACUS-first adapter with explicit format selection."""

    def read_output(self, path):
        raise UnsupportedFormatError("running_scf/running_lr metadata parsing is planned; supply metadata explicitly")

    def read_wavefunctions(self, path, *, format="gamma_text", spin="unspecified"):
        if format != "gamma_text":
            raise UnsupportedFormatError("only one-frame real Gamma LCAO text is currently supported")
        return read_gamma_wavefunctions(path, spin=spin)

    def read_overlap(self, path, *, format):
        if format == "stateid_npy":
            return read_overlap_npy(path)
        raise UnsupportedFormatError("native ABACUS S(R)/S(k) reading is planned; use a validated dense stateid_npy S")

    def read_lr_eigenvectors(self, path, *, format):
        if format == "stateid_npz":
            return read_lr_npz(path)
        raise UnsupportedFormatError("native ABACUS LR amplitudes need MPI distribution and ph-order metadata; adapter pending")
