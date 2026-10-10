"""Audit a one-step scalar NV ABACUS calculation from its actual outputs.

OPENBLAS_NUM_THREADS=1 python examples/abacus_non_gamma.py RUN --output report.json
Requires stateid[symmetry]. No ABACUS execution or writes to its source tree.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from stateid.io import (read_structure, read_kpoints_info, read_input_parameters,
                        read_scf_log, read_multik_wavefunctions, read_ao_labels, read_csr)
from stateid.symmetry import (symmetry_from_structure, spgrep_irreps,
    build_bloch_ao_operation, analyze_little_group, c3v_labels,
    project_representation, project_sewing_matrix)


def maxabs(a):
    return float(np.max(np.abs(a)))


def audit(run, *, suffix="stateid", atol=1e-5, energy_tol_ry=1e-6):
    run = Path(run)
    out = run / f"OUT.{suffix}"
    metadata = read_scf_log(out / "running_scf.log")
    if not metadata["converged"]:
        raise ValueError("SCF did not converge")
    inp = read_input_parameters(out / "INPUT.info")
    if inp["basis_type"] != "lcao" or inp["calculation"] != "scf" or int(inp["nspin"]) != 2:
        raise ValueError("this NV audit example expects a scalar nspin=2 LCAO SCF")
    if inp.get("gamma_only", "0") not in {"0", "false"}:
        raise ValueError("complex multi-k output is required")
    cell = read_structure(run / "STRU")
    final = read_structure(out / "STRU_FINAL")
    if (cell.species != final.species
            or maxabs(cell.lattice_angstrom-final.lattice_angstrom) > 1e-7
            or maxabs(cell.fractional_positions-final.fractional_positions) > 1e-7):
        raise ValueError("input/final structure mismatch")
    points = read_kpoints_info(out / "KPT.info")
    wfc = read_multik_wavefunctions(out / "WFC", kpoints=points,
             lattice_vectors=cell.lattice_vectors, nspin=int(inp["nspin"]))
    labels = read_ao_labels(out / "Orbital")
    sr = read_csr(out / "srg1_nao.csr")
    overlaps = sr.to_k(points.fractional)
    hamiltonians = [read_csr(out / f"hrs{i}g1_nao.csr").to_k(points.fractional) for i in (1, 2)]
    if sr.basis_size != len(labels):
        raise ValueError("CSR/AO label dimensions differ")
    # CSR cell header has only six significant digits; never use it to convert kvec_c.
    cell_error = maxabs(sr.lattice_vectors * sr.lattice_constant_angstrom - cell.lattice_angstrom)
    if cell_error > 1e-5:
        raise ValueError("CSR/STRU physical cell mismatch beyond printed precision")
    rotations, translations = symmetry_from_structure(
        cell.lattice_angstrom, cell.fractional_positions, cell.species_ids)
    refs = [spgrep_irreps(rotations, translations, k) for k in points.fractional]
    operations = [[build_bloch_ao_operation(
        cell.fractional_positions, cell.lattice_angstrom, cell.species, labels, w, k,
        translation_fractional=t) for w, t in zip(rotations, translations)]
        for k in points.fractional]
    errors = dict(normalization=0., eigen_equation_ry=0., ao_metric=0.,
                  subspace_closure=0., group_law=0., character_reconstruction=0.,
                  cross_k_closure=0., cross_k_metric=0., hamiltonian_covariance_ry=0.)
    analyses, cross_count = [], 0
    for spin, spin_wfc in enumerate(wfc):
        for ik, data in enumerate(spin_wfc):
            c, sk, hk = data.coefficients, overlaps[ik], hamiltonians[spin][ik]
            if c.shape[0] != len(labels):
                raise ValueError("WFC/AO label dimensions differ")
            errors["normalization"] = max(errors["normalization"], maxabs(c.conj().T @ sk @ c-np.eye(c.shape[1])))
            errors["eigen_equation_ry"] = max(errors["eigen_equation_ry"],
                    maxabs(hk @ c - (sk @ c) * data.energies_ry))
            ref = refs[ik]
            ops = [operations[ik][j] for j in ref.group.indices]
            # Energy clustering is only a candidate selection. Every candidate
            # must independently pass metric, closure, group-law and character checks.
            clusters = np.split(np.arange(c.shape[1]), np.flatnonzero(np.diff(data.energies_ry) > energy_tol_ry) + 1)
            names = c3v_labels(ref).labels if ref.group.order == 6 else None
            for bands in clusters:
                result = analyze_little_group(c[:, bands], ops, ref, sk, atol=atol)
                if not result.match.valid:
                    raise ValueError("characters do not decompose into irreps")
                errors["ao_metric"] = max(errors["ao_metric"], *(r.ao_metric_error for r in result.representations))
                errors["subspace_closure"] = max(errors["subspace_closure"], *(r.closure_error for r in result.representations))
                errors["group_law"] = max(errors["group_law"], result.group_error)
                errors["character_reconstruction"] = max(errors["character_reconstruction"], result.match.residual)
                analyses.append(dict(spin=data.spin, k_index=ik+1, bands=(bands+1).tolist(),
                    energy_ry=data.energies_ry[bands].tolist(),
                    characters_real=result.characters.real.tolist(),
                    characters_imag=result.characters.imag.tolist(),
                    eigenvalues_real=[r.eigenvalues.real.tolist() for r in result.representations],
                    eigenvalues_imag=[r.eigenvalues.imag.tolist() for r in result.representations],
                    operation_rotations=ref.group.rotations.tolist(),
                    operation_translations=ref.group.translations.tolist(),
                    irreps={names[key] if names else key: n for key, n in result.match.multiplicities.items()}))
            for op in operations[ik]:
                delta = points.fractional-op.mapped_kpoint
                targets = np.flatnonzero(np.max(abs(delta-np.rint(delta)), axis=1) < 1e-8)
                if len(targets) != 1:
                    raise ValueError("mesh must contain exactly one representative of every star arm")
                jk = int(targets[0])
                b = op.matrix
                errors["hamiltonian_covariance_ry"] = max(errors["hamiltonian_covariance_ry"],
                    maxabs(b.conj().T @ hamiltonians[spin][jk] @ b-hk))
                if ik != jk:
                    # This calculation has a gap after band 16 in both spin blocks.
                    sew = project_sewing_matrix(c[:, :16], spin_wfc[jk].coefficients[:, :16],
                                                b, sk, overlaps[jk])
                    errors["cross_k_closure"] = max(errors["cross_k_closure"], sew.closure_error)
                    errors["cross_k_metric"] = max(errors["cross_k_metric"], sew.ao_metric_error)
                    cross_count += 1
    if max(errors.values()) > atol:
        raise ValueError(f"validation tolerance exceeded: {errors}")

    # In this origin all six NV operations have zero atom shifts. Exercise
    # the essential Bloch phase by choosing another periodic image of atom 0,
    # with the EXACT accompanying basis gauge change of the real C and S.
    ik = 0
    k, ref, sk = points.fractional[ik], refs[ik], overlaps[ik]
    if ref.group.order != 6:
        raise ValueError("first k point of this NV example must have C3v little group")
    image_shifts = np.zeros_like(cell.fractional_positions)
    image_shifts[0] = [1, 0, 0]
    ao_shifts = image_shifts[[label.atom_index for label in labels]]
    v = np.exp(2j * np.pi * (ao_shifts @ k))
    new_s = v[:, None] * sk * v.conj()[None, :]
    c = wfc[1][ik].coefficients[:, 14:16]  # beta bands 15/16, a complete E pair
    new_c = v[:, None] * c
    rewrapped_ops = [build_bloch_ao_operation(
        cell.fractional_positions+image_shifts, cell.lattice_angstrom, cell.species,
        labels, w, k, translation_fractional=t)
        for w, t in zip(ref.group.rotations, ref.group.translations)]
    rewrapped = analyze_little_group(new_c, rewrapped_ops, ref, new_s, atol=atol)
    baseline = analyze_little_group(c, [operations[ik][j] for j in ref.group.indices], ref, sk, atol=atol)
    wrong = [project_representation(new_c, op.ao_operation.matrix, new_s,
             operator_kind="coefficient") for op in rewrapped_ops]
    phase_check = dict(character_change=maxabs(rewrapped.characters-baseline.characters),
                      closure=max(r.closure_error for r in rewrapped.representations),
                      omitted_phase_metric_error=max(r.ao_metric_error for r in wrong),
                      omitted_phase_closure=max(r.closure_error for r in wrong))
    if phase_check["character_change"] > atol or phase_check["omitted_phase_closure"] < .1:
        raise ValueError("Bloch phase positive/negative control failed")
    files = [run/"INPUT", run/"STRU", run/"KPT", out/"INPUT.info", out/"KPT.info",
             out/"STRU_FINAL", out/"Orbital", out/"srg1_nao.csr",
             out/"hrs1g1_nao.csr", out/"hrs2g1_nao.csr", out/"running_scf.log"]
    files += [Path(d.source) for row in wfc for d in row]
    return dict(metadata=metadata, atoms=len(cell.species), nao=len(labels),
        kpoints=points.fractional.tolist(), little_group_orders=[r.group.order for r in refs],
        raw_weighted_electron_counts=[sum(float(d.occupations.sum()) for d in row) for row in wfc],
        tolerance=atol, candidate_energy_tolerance_ry=energy_tol_ry,
        csr_cell_rounding_error_angstrom=cell_error, errors=errors,
        accepted_subspaces=len(analyses), cross_k_maps_checked=cross_count,
        bloch_phase_control=phase_check, bands=analyses,
        sha256={str(p.relative_to(run)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--suffix", default="stateid")
    args = parser.parse_args()
    report = audit(args.run, suffix=args.suffix)
    rendered = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered)
    print(json.dumps({k: v for k, v in report.items() if k not in {"bands", "sha256"}}, indent=2))
