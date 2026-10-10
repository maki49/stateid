"""ABACUS adapter: common Gamma/multi-k text parser and CSR overlap IO."""

from pathlib import Path
import re

import numpy as np

from stateid._validation import metric
from .arrays import read_lr_npz, read_overlap_npy
from .csr import read_csr
from .models import OrbitalData
from .errors import UnsupportedFormatError
from .metadata import read_structure, read_scf_log


def read_wavefunctions(path, *, format="auto", spin="unspecified"):
    """Read one real-Gamma or complex-k LCAO text file, one frame only.

    Multi-k files prepend '(index of k points)' and three Cartesian k values
    (2pi/lat0 units); each AO coefficient is a real/imaginary pair. Common band
    metadata and coefficient handling are shared with the real Gamma reader.
    Spin and scalar-vs-spinor setup must be established by the caller.
    """
    if format not in {"auto", "gamma_text", "k_text"}:
        raise UnsupportedFormatError("supported wavefunction formats: auto, gamma_text, k_text (text only)")
    if spin not in {"alpha", "beta", "unspecified"}:
        raise ValueError("spin must be alpha, beta, or unspecified")
    try:
        lines = [line.strip() for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip()]
    except UnicodeError as exc:
        raise UnsupportedFormatError("binary wavefunctions are not supported") from exc
    is_complex = bool(lines and re.fullmatch(r"\d+\s+\(index of k points\)", lines[0]))
    is_gamma = bool(lines and re.fullmatch(r"\d+\s+\(number of bands\)", lines[0]))
    if (not (is_complex or is_gamma) or (format == "gamma_text" and not is_gamma)
            or (format == "k_text" and not is_complex)):
        raise UnsupportedFormatError("wavefunction text header does not match the requested format")
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

    def numbers(count):
        nonlocal pos
        values = []
        while len(values) < count:
            if pos >= len(lines):
                raise ValueError("truncated wavefunction numeric data")
            try:
                values.extend(float(t.replace("D", "E").replace("d", "e")) for t in lines[pos].split())
            except ValueError as exc:
                raise ValueError(f"line {pos + 1}: expected numeric data") from exc
            pos += 1
        if len(values) != count or not np.all(np.isfinite(values)):
            raise ValueError("wrong numeric count or non-finite value")
        return np.asarray(values)

    k_index, k_cartesian = None, None
    if is_complex:
        k_index = field("index of k points", True)
        if k_index < 1:
            raise ValueError("k index must be positive and one-based")
        k_cartesian = numbers(3)
    nb, nao = field("number of bands", True), field("number of orbitals", True)
    if min(nb, nao) < 1:
        raise ValueError("band and AO dimensions must be positive")
    columns, energies, occupations = [], [], []
    for band in range(1, nb + 1):
        if field("band", True) != band:
            raise ValueError("band labels must be consecutive and one-based")
        energies.append(field("Ry"))
        occupations.append(field("Occupations"))
        values = numbers(nao * (2 if is_complex else 1))
        columns.append(values[::2] + 1j * values[1::2] if is_complex else values)
    if pos != len(lines):
        raise UnsupportedFormatError("extra data after one frame; append-mode/multiple frames are not supported")
    return OrbitalData(np.array(columns).T, np.array(energies), np.array(occupations),
                       np.arange(1, nb + 1), spin, str(path), k_index, k_cartesian)


def read_gamma_wavefunctions(path, *, spin="unspecified"):
    """Compatibility wrapper requiring a one-frame real Gamma text file."""
    return read_wavefunctions(path, format="gamma_text", spin=spin)


class AbacusReader:
    """ABACUS-first adapter with explicit overlap format selection."""

    def read_output(self, path):
        return read_scf_log(path)

    def read_structure(self, path, **kwargs):
        return read_structure(path, **kwargs)

    def read_wavefunctions(self, path, *, format="auto", spin="unspecified"):
        return read_wavefunctions(path, format=format, spin=spin)

    def read_overlap(self, path, *, format, kpoints=None, frame=None, atol=1e-8):
        """Return S(k) for one or many reciprocal fractional k points.

        abacus_csr autodetects old/new TEXT CSR. Omitted kpoints means Gamma.
        stateid_npy is already a single S(k), so kpoints/frame must be omitted.
        Checks Hermiticity and positive definiteness at every requested k.
        """
        if format == "stateid_npy":
            if kpoints is not None or frame is not None:
                raise ValueError("stateid_npy already contains one S(k); kpoints/frame are inapplicable")
            return read_overlap_npy(path)
        if format == "abacus_csr":
            realspace = read_csr(path, frame=frame)
            s = realspace.to_k([0, 0, 0] if kpoints is None else kpoints)
            if s.ndim == 2:
                return metric(s, realspace.basis_size, atol)
            return np.stack([metric(sk, realspace.basis_size, atol) for sk in s])
        raise UnsupportedFormatError("supported overlap formats: stateid_npy, abacus_csr (old/new text)")

    def read_lr_eigenvectors(self, path, *, format):
        if format == "stateid_npz":
            return read_lr_npz(path)
        raise UnsupportedFormatError("native ABACUS LR amplitudes need MPI distribution and ph-order metadata; adapter pending")


def read_multik_wavefunctions(directory, *, kpoints, lattice_vectors, nspin,
                             ionic_step=1, atol=1e-7):
    """Read modern complex WFC files as (spin, spatial k).

    nspin=1 uses wfk{k}g{step}_nao.txt; nspin=2 adds s{spin}.

    nspin must come from this run's INPUT/INPUT.info, never from filenames.
    Only scalar nspin=1/2, complex text, one ionic step. kpoints is KPointData;
    lattice_vectors must be the ORIGINAL input's dimensionless lattice, since
    WFC kvec_c uses that lat0. Cross-check both global header indices and k.
    Legacy names remain available through read_wavefunctions individually.
    """
    from stateid.realspace import k_cartesian_to_fractional
    from stateid._validation import positive_tolerance
    positive_tolerance(atol)
    if isinstance(nspin, bool) or not isinstance(nspin, (int, np.integer)) or nspin not in (1, 2):
        raise UnsupportedFormatError("multi-k collection requires scalar nspin=1 or 2")
    if isinstance(ionic_step, bool) or not isinstance(ionic_step, (int, np.integer)) or ionic_step < 1:
        raise ValueError("ionic_step must be a positive one-based integer")
    directory = Path(directory)
    nk = len(kpoints.indices)
    result = []
    shape = None
    for ispin in range(nspin):
        spin = ("alpha", "beta")[ispin] if nspin == 2 else "unspecified"
        row = []
        for ik, point in enumerate(kpoints.fractional, 1):
            spin_suffix = f"s{ispin+1}" if nspin == 2 else ""
            data = read_wavefunctions(
                directory / f"wfk{ik}{spin_suffix}g{ionic_step}_nao.txt",
                format="k_text", spin=spin)
            actual_k = k_cartesian_to_fractional(data.k_cartesian, lattice_vectors)
            if data.k_index != ispin * nk + ik or not np.allclose(actual_k, point, atol=atol, rtol=0):
                raise ValueError(f"{data.source}: k index/coordinate disagrees with metadata")
            shape = data.coefficients.shape if shape is None else shape
            if data.coefficients.shape != shape:
                raise ValueError("wavefunction dimensions differ across spin/k blocks")
            row.append(data)
        result.append(tuple(row))
    return tuple(result)
