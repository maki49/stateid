"""Small, explicit ABACUS metadata readers; no inferred spin or cell settings."""
from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np

from .errors import UnsupportedFormatError

# Match source_base/constants.h in the validated ABACUS writer, not CODATA.
ABACUS_BOHR_TO_ANGSTROM = 0.5291770


def _lines(path):
    return [clean for line in Path(path).read_text(encoding="utf-8").splitlines()
            if (clean := re.split(r"[#!]", line, maxsplit=1)[0].strip())]


def _floats(tokens):
    values = np.array([float(t.replace("D", "E").replace("d", "e")) for t in tokens])
    if not np.all(np.isfinite(values)):
        raise ValueError("metadata contains non-finite numbers")
    return values


@dataclass(frozen=True)
class AbacusStructure:
    lattice_constant_bohr: float
    lattice_vectors: np.ndarray  # dimensionless, ROWS, in this file's lat0
    fractional_positions: np.ndarray  # unwrapped, original atom order
    species: tuple[str, ...]  # one label per atom, not inferred atomic numbers
    coordinate_mode: str
    atom_fields: tuple[tuple[str, ...], ...]  # uninterpreted trailing m/v/f/mag...
    source: str
    bohr_to_angstrom: float = ABACUS_BOHR_TO_ANGSTROM

    @property
    def lattice_angstrom(self):
        return self.lattice_vectors * self.lattice_constant_bohr * self.bohr_to_angstrom

    @property
    def species_ids(self):
        """Positive type IDs for spglib; these are NOT atomic numbers."""
        types = {s: i + 1 for i, s in enumerate(dict.fromkeys(self.species))}
        return np.array([types[s] for s in self.species])


def read_structure(path, *, bohr_to_angstrom=ABACUS_BOHR_TO_ANGSTROM):
    """Read one explicit-cell STRU/STRU_FINAL, retaining its original AO order.

    Direct, Cartesian (lat0 units), Cartesian_au and Cartesian_angstrom are
    supported. Per-atom trailing fields are retained, not used as a magnetic
    symmetry model. Centered coordinates and implicit Bravais cells fail.
    STRU_FINAL may rescale lat0: use the ORIGINAL STRU for WFC kvec_c conversion.
    """
    if not np.isfinite(bohr_to_angstrom) or bohr_to_angstrom <= 0:
        raise ValueError("bohr_to_angstrom must be positive")
    sections, current = {}, None
    headers = {"ATOMIC_SPECIES", "NUMERICAL_ORBITAL", "NUMERICAL_DESCRIPTOR",
               "LATTICE_CONSTANT", "LATTICE_VECTORS", "LATTICE_PARAMETERS",
               "ATOMIC_POSITIONS", "PAW_FILES"}
    for line in _lines(path):
        if line in headers:
            if line in sections:
                raise UnsupportedFormatError("multiple STRU sections/frames are not supported")
            sections[line] = []
            current = line
        elif current is None:
            raise ValueError("expected a STRU section header")
        else:
            sections[current].append(line)
    required = {"ATOMIC_SPECIES", "LATTICE_CONSTANT", "LATTICE_VECTORS", "ATOMIC_POSITIONS"}
    if not required <= sections.keys() or "LATTICE_PARAMETERS" in sections:
        raise UnsupportedFormatError("STRU requires an explicit lattice and atom positions")
    records = [line.split() for line in sections["ATOMIC_SPECIES"]]
    if not records or any(len(row) < 2 for row in records):
        raise ValueError("malformed ATOMIC_SPECIES")
    types = [row[0] for row in records]
    if len(set(types)) != len(types):
        raise ValueError("duplicate species label")
    if np.any(_floats([row[1] for row in records]) <= 0):
        raise ValueError("atomic masses must be positive")
    lat0 = _floats(" ".join(sections["LATTICE_CONSTANT"]).split())
    lattice_rows = [line.split() for line in sections["LATTICE_VECTORS"]]
    if len(lat0) != 1 or lat0[0] <= 0 or len(lattice_rows) != 3 or any(len(r) != 3 for r in lattice_rows):
        raise ValueError("invalid lattice metadata")
    lattice = np.array([_floats(r) for r in lattice_rows])
    if np.linalg.matrix_rank(lattice) != 3:
        raise ValueError("singular lattice")
    rows = iter(sections["ATOMIC_POSITIONS"])
    positions, species, fields = [], [], []
    try:
        mode = next(rows)
        if mode not in {"Direct", "Cartesian", "Cartesian_au", "Cartesian_angstrom"}:
            raise UnsupportedFormatError(f"unsupported STRU coordinates: {mode}")
        for label in types:
            if next(rows) != label:
                raise ValueError("atom blocks must follow ATOMIC_SPECIES order")
            mag = _floats(next(rows).split())
            if mag.size != 1:
                raise ValueError("species magnetization must be one scalar")
            count = int(next(rows))
            if count < 1:
                raise ValueError("each species must contain atoms")
            for _ in range(count):
                tokens = next(rows).split()
                if len(tokens) < 3:
                    raise ValueError("atom position needs three coordinates")
                positions.append(_floats(tokens[:3]))
                species.append(label)
                fields.append(tuple(tokens[3:]))
    except StopIteration as exc:
        raise ValueError("truncated ATOMIC_POSITIONS") from exc
    if next(rows, None) is not None:
        raise ValueError("extra data after ATOMIC_POSITIONS")
    positions = np.array(positions)
    if mode != "Direct":
        scale = {"Cartesian": 1, "Cartesian_au": lat0[0],
                 "Cartesian_angstrom": lat0[0] * bohr_to_angstrom}[mode]
        positions = np.linalg.solve(lattice.T, (positions / scale).T).T
    return AbacusStructure(float(lat0[0]), lattice, positions, tuple(species),
                           mode, tuple(fields), str(path), bohr_to_angstrom)


@dataclass(frozen=True)
class KPointData:
    indices: np.ndarray  # one-based spatial k labels, BEFORE spin duplication
    fractional: np.ndarray
    weights: np.ndarray  # raw printed values; no normalization/occupancy inference
    source: str


def read_kpoints_info(path):
    """Read one KPT.info DIRECT table, not a KPT mesh specification.

    Rounded weights are retained as printed; do not demand that they sum
    exactly to one. SOC/spin duplication is not encoded in this table.
    """
    lines = [line.strip() for line in Path(path).read_text().splitlines() if line.strip()]
    if len(lines) < 4 or not re.fullmatch(r"nkstot now\s*=\s*\d+", lines[0]):
        raise UnsupportedFormatError("expected one ABACUS KPT.info table")
    count = int(lines[0].split("=")[1])
    if (lines[1] != "K-POINTS DIRECT COORDINATES"
            or lines[2].split() != ["KPOINTS", "DIRECT_X", "DIRECT_Y", "DIRECT_Z", "WEIGHT"]):
        raise UnsupportedFormatError("only the KPT.info DIRECT table is supported")
    rows = [line.split() for line in lines[3:]]
    if count < 1 or len(rows) != count or any(len(r) != 5 for r in rows):
        raise ValueError("KPT.info row count/shape mismatch")
    indices = np.array([int(r[0]) for r in rows])
    data = np.array([_floats(r[1:]) for r in rows])
    if not np.array_equal(indices, np.arange(1, count + 1)) or np.any(data[:, 3] < 0):
        raise ValueError("invalid KPT.info index/weight")
    return KPointData(indices, data[:, :3], data[:, 3], str(path))


def read_input_parameters(path):
    """Read INPUT/INPUT.info as raw string values; no undocumented defaults."""
    rows = _lines(path)
    if not rows or rows[0] != "INPUT_PARAMETERS":
        raise UnsupportedFormatError("expected INPUT_PARAMETERS")
    values = {}
    for line in rows[1:]:
        tokens = line.split(maxsplit=1)
        if not re.fullmatch(r"[a-z][a-z0-9_]*", tokens[0]) or tokens[0] in values:
            raise ValueError("invalid or duplicate parameter")
        values[tokens[0]] = tokens[1] if len(tokens) == 2 else ""
    return values


def read_scf_log(path):
    """Read status/provenance of one ABACUS running_scf.log SCF calculation.

    Absence of a convergence marker means False (not assumed success).
    Multiple ionic steps/runs are rejected. This is not an LR log reader.
    """
    text = Path(path).read_text()
    versions = re.findall(r"^\s*ABACUS\s+(\S+)\s*$", text, re.M)
    finals = re.findall(r"!FINAL_ETOT_IS\s+(\S+)\s+eV", text)
    iterations = re.findall(r"#ION MOVE#\s+(\d+)\s+#ELEC ITER#\s+(\d+)", text)
    if len(versions) != 1 or not iterations:
        raise UnsupportedFormatError("expected a single running_scf.log with SCF iterations")
    steps = {int(ion) for ion, _ in iterations}
    if len(steps) != 1 or len(finals) > 1 or text.count("#SCF IS CONVERGED#") > 1:
        raise UnsupportedFormatError("multiple ionic steps/runs require an explicit frame parser")
    commits = re.findall(r"^\s*Commit:\s+(\S+)", text, re.M)
    return {"version": versions[0], "commit": commits[0] if commits else None,
            "converged": "#SCF IS CONVERGED#" in text,
            "iterations": int(iterations[-1][1]),
            "energy_ev": float(_floats(finals)[0]) if finals else None,
            "source": str(path)}
