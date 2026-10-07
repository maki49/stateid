"""Read the CellIndex::write_orb_info ABACUS Orbital table."""
from pathlib import Path

from stateid.symmetry.ao import AOLabel
from stateid.symmetry.harmonics import magnetic_order


def read_ao_labels(path):
    """Return labels in coefficient-row order, without trusting the #io legend.

    The first field is zero-based atom iat, NOT global AO index. The m field
    is a real-component index 0..2l, mapped here to signed magnetic_order(l).
    zeta is retained one-based. Completeness is checked by the AO builder.
    """
    labels = []
    for lineno, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        fields = line.split()
        if len(fields) != 6:
            raise ValueError(f"Orbital line {lineno}: expected six fields")
        try:
            atom, ell, component, zeta = (int(fields[i]) for i in (0, 2, 3, 4))
            order = magnetic_order(ell)
            if atom < 0 or zeta < 1 or not 0 <= component < len(order):
                raise ValueError("invalid index")
        except ValueError as exc:
            raise ValueError(f"Orbital line {lineno}: invalid AO indices") from exc
        labels.append(AOLabel(atom, fields[1], ell, zeta, order[component]))
    if not labels:
        raise ValueError("Orbital table contains no AO labels")
    return tuple(labels)
