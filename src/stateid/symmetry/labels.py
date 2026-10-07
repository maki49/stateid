"""Conventional names are matched to characters, never to enumeration order."""
from dataclasses import dataclass
from importlib.metadata import version
import numpy as np
from stateid._validation import positive_tolerance
from .groups import _operations, real_array
from .tables import C3V


@dataclass(frozen=True)
class IrrepLabels:
    labels: dict[str, str]  # run-local irrep ID -> conventional label
    source: str
    residuals: dict[str, float]


def match_irrep_labels(reference, label_characters, *, source, atol=1e-5):
    """Match a complete named table in exactly the reference operation order.

    This is a low-level adapter boundary: the caller certifies the same cell,
    origin, k, Seitz representatives and phase convention. Use
    labels_from_character_table when operation alignment is needed.
    """
    positive_tolerance(atol)
    if not isinstance(source, str) or not source.strip():
        raise ValueError("label provenance must be a nonempty source string")
    names = list(label_characters)
    rows = np.asarray(list(label_characters.values()), dtype=complex)
    if (len(names) != len(reference.matrices)
            or rows.shape != (len(names), reference.group.order)
            or not np.all(np.isfinite(rows))
            or any(not isinstance(n, str) or not n for n in names)):
        raise ValueError("provide a complete finite named character table")
    labels, residuals = {}, {}
    for identifier, chi in reference.characters.items():
        errors = np.max(abs(rows - chi), axis=1)
        matches = np.flatnonzero(errors <= atol)
        if len(matches) != 1:
            raise ValueError(f"{identifier}: conventional label is missing or ambiguous")
        j = matches[0]
        labels[identifier], residuals[identifier] = names[j], float(errors[j])
    if len(set(labels.values())) != len(labels):
        raise ValueError("label mapping is not one-to-one")
    return IrrepLabels(labels, source, residuals)


def labels_from_character_table(reference, *, rotations, translations, kpoint,
                                characters, source, atol=1e-5):
    """Align an external per-operation table in the SAME cell and origin.

    Reorders operations and corrects integer representative differences:
    t_stateid=t_table+L => chi_stateid=exp(-2 pi i k.L) chi_table.
    No undocumented cell/origin or k-star conversion is attempted.
    """
    w, t = _operations(rotations, translations, atol)
    k = real_array(kpoint, (3,), "table kpoint")
    g = reference.group
    if not np.allclose(k, g.kpoint, atol=atol, rtol=0):
        raise ValueError("table kpoint must match exactly in the same reciprocal basis")
    if len(w) != g.order:
        raise ValueError("table and reference must describe the same little group")
    indices, phases = [], []
    for wg, tg in zip(g.rotations, g.translations):
        matches = np.flatnonzero(np.all(w == wg, axis=(1, 2)))
        if len(matches) != 1:
            raise ValueError("table and reference rotations differ")
        j = matches[0]
        delta = tg - t[j]
        if not np.allclose(delta, np.rint(delta), atol=atol, rtol=0):
            raise ValueError("table translations differ beyond lattice representatives; check origin")
        indices.append(j)
        phases.append(np.exp(-2j * np.pi * np.dot(k, np.rint(delta))))
    aligned = {}
    for name, chi in characters.items():
        chi = np.asarray(chi, dtype=complex)
        if chi.shape != (len(w),):
            raise ValueError("one character per table operation is required")
        aligned[name] = chi[indices] * phases
    return match_irrep_labels(reference, aligned, source=source, atol=atol)


def c3v_labels(reference, *, atol=1e-5):
    """Mulliken A1/A2/E at Gamma, classified by rotations, independent of axes."""
    positive_tolerance(atol)
    g = reference.group
    if g.order != 6 or not np.allclose(g.kpoint, 0, atol=atol, rtol=0):
        raise ValueError("built-in C3v labels require its six operations at Gamma")
    classes = []
    for w in g.rotations:
        if np.array_equal(w, np.eye(3)):
            classes.append(0)
        elif (round(np.linalg.det(w)) == 1 and np.trace(w) == 0
              and np.array_equal(w @ w @ w, np.eye(3))):
            classes.append(1)
        elif (round(np.linalg.det(w)) == -1 and np.trace(w) == 1
              and np.array_equal(w @ w, np.eye(3))):
            classes.append(2)
        else:
            raise ValueError("rotations do not form C3v")
    if [classes.count(i) for i in range(3)] != [1, 2, 3]:
        raise ValueError("rotations do not have C3v classes")
    chars = {name: np.asarray(row)[classes] for name, row in C3V.irreps.items()}
    return match_irrep_labels(reference, chars, source="stateid C3v Mulliken table", atol=atol)


def irrep_labels(reference, lattice, positions, numbers, kpoint_name, *,
                 symprec=1e-6, atol=1e-5):
    """Optional IrRep 3.3.0 / irreptables 3.1.0 standard-label adapter.

    Supply the SAME cell/origin as reference; IrRep handles conversion to its
    BCS reference setting. Only scalar, unitary space groups and table-covered
    maximal k points are supported. Missing names/cells fail explicitly.
    This adapter uses a version-guarded internal IrRep API, isolated here.
    """
    positive_tolerance(symprec, "symprec")
    positive_tolerance(atol)
    try:
        from irrep.spacegroup_irreps import SpaceGroupIrreps
        from irreptables.irreps import IrrepTable  # ensure the optional tables exist
    except ImportError as exc:
        raise ImportError("install stateid[labels] for IrRep standard labels") from exc
    versions = (version("irrep"), version("irreptables"))
    if versions != ("3.3.0", "3.1.0"):
        raise RuntimeError("IrRep label adapter validated for irrep==3.3.0, "
                           "irreptables==3.1.0; install stateid[labels]")
    if not isinstance(kpoint_name, str) or not kpoint_name:
        raise ValueError("provide an explicit IrRep table k-point name, e.g. GM")
    try:
        sg = SpaceGroupIrreps.from_cell(
            cell=(lattice, positions, numbers), spinor=False, include_TR=False,
            search_cell=True, symprec=symprec)
        table = sg.get_irreps_from_table(kpoint_name, reference.group.kpoint)
    except (ValueError, RuntimeError, KeyError, AssertionError) as exc:
        raise ValueError(f"IrRep cell/table lookup failed for {kpoint_name}: {exc}") from exc
    keys = set(next(iter(table.values())))
    if any(set(row) != keys for row in table.values()):
        raise ValueError("IrRep character rows cover different operation sets")
    ops = [op for op in sg.symmetries if op.ind in keys]
    if len(ops) != len(keys):
        raise ValueError("IrRep operation indices do not align with its character table")
    return labels_from_character_table(
        reference, rotations=[op.rotation for op in ops],
        translations=[op.translation for op in ops], kpoint=reference.group.kpoint,
        characters={name: [row[op.ind] for op in ops] for name, row in table.items()},
        source=f"IrRep {versions[0]} / irreptables {versions[1]}, "
               f"BCS SG {sg.number_str}, k={kpoint_name}", atol=atol)
