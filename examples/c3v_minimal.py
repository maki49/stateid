"""Synthetic A1+E orbitals; this is not an NV- or ABACUS calculation."""

import numpy as np

from stateid.spin import determinant_s2
from stateid.symmetry import analyze_c3v


def operations():
    theta = 2 * np.pi / 3
    r = np.eye(3)
    r[1:, 1:] = [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
    f = np.diag([1.0, 1.0, -1.0])
    return {"E": np.eye(3), "C3": r, "C3^2": r @ r,
            "sv1": f, "sv2": r @ f, "sv3": r @ r @ f}


def main():
    for label, columns in [("A1", [0]), ("E", [1, 2])]:
        result = analyze_c3v(np.eye(3)[:, columns], operations(), operator_kind="coefficient")
        print(f"{label}: characters={np.round(result.characters.real, 8)}, identified={result.match.label}")
        assert result.match.irrep == label
    triplet = determinant_s2(np.eye(2), np.empty((2, 0)))
    print(f"two parallel spins: <S^2>={triplet:.1f}")


if __name__ == "__main__":
    main()
