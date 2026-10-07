"""Synthetic scalar AO examples; no ABACUS wavefunctions are assumed.

python examples/bloch_little_group.py
python examples/bloch_little_group.py --irrep-labels
"""
import argparse
import numpy as np
from stateid.symmetry import (
    AOLabel, build_bloch_ao_operation, spgrep_irreps, analyze_little_group,
    c3v_labels, symmetry_from_structure, irrep_labels,
)


def main(use_irrep=False):
    # Gamma C3v: z is A1, the complete (-x,-y) subspace is E.
    r = np.array([[0,-1,0], [1,-1,0], [0,0,1]])
    f = np.array([[1,-1,0], [0,-1,0], [0,0,1]])
    rotations = np.array([np.eye(3, dtype=int), r, r@r, f, r@f, r@r@f])
    lattice = np.array([[1,0,0], [-.5,np.sqrt(3)/2,0], [0,0,2]])
    reference = spgrep_irreps(rotations, np.zeros((6,3)), [0,0,0])
    labels = [AOLabel(0, "C", 1, 1, m) for m in (0,1,-1)]
    operations = [build_bloch_ao_operation(
        [[0,0,0]], lattice, ["C"], labels, w, reference.group.kpoint)
        for w in reference.group.rotations]
    names = c3v_labels(reference)
    for columns in ([0], [1,2], [0,1,2]):
        result = analyze_little_group(np.eye(3)[:, columns], operations, reference)
        named = {names.labels[key]: count for key, count in result.match.multiplicities.items()}
        print("Gamma:", named, "characters:", np.round(result.characters.real, 8))

    # Geometric Pc (SG 7) in an a-glide setting. All AOs and C are synthetic.
    lattice = np.array([[1,0,0], [0,1.3,0], [.2,0,1.7]])
    positions = [[.13,.17,.23], [.63,-.17,.23], [.27,.32,.41], [.77,-.32,.41]]
    species, numbers = ["C","C","N","N"], [1,1,2,2]
    rotations, translations = symmetry_from_structure(lattice, positions, numbers)
    k = np.array([-.5,0,0])  # B in IrRep's reference setting for this cell
    reference = spgrep_irreps(rotations, translations, k)
    labels = [AOLabel(i,s,0,1,0) for i,s in enumerate(species)]
    operations = [build_bloch_ao_operation(
        positions, lattice, species, labels, w, k, translation_fractional=t)
        for w,t in zip(reference.group.rotations, reference.group.translations)]
    c = np.array([[1j],[1],[0],[0]]) / np.sqrt(2)
    result = analyze_little_group(c, operations, reference, np.eye(4))
    print("Glide boundary:", result.match.label,
          "characters:", np.round(result.characters, 8))
    print("Group residual:", result.group_error)
    if use_irrep:
        names = irrep_labels(reference, lattice, positions, numbers, "B")
        print("Standard label:", names.labels[result.match.irrep], "source:", names.source)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--irrep-labels", action="store_true")
    main(parser.parse_args().irrep_labels)
