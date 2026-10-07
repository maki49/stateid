# stateid

**English** | [简体中文](README.zh-CN.md)

**Evidence-based electronic-state identification from first-principles calculations, starting with ABACUS.**

`stateid` stands for electronic **state identification**. It combines symmetry,
spin, transition composition, and other physical evidence to explain what an
electronic state is and how strongly that assignment is supported. The name is
independent of a particular code, ground state, or defect. The package has not
yet been published to PyPI.

## Scope and current capabilities

The goal is to identify defect states, molecular states, and localized electronic
states in solids using energies, orbital and excited-state symmetry, spin,
transition composition, and localization. Tutorials explain the equations,
assumptions, and diagnostics as well as file handling. The initial focus is
**ABACUS**, with numerical physics separated from software-specific IO.

Version `0.2.0` supports ground-state orbital and Slater-determinant analysis,
legacy and modern ABACUS text CSR, Fourier transforms at arbitrary k points,
and complex multi-k wavefunction text. LR/TDA support includes a data model,
a standardized exchange format, explicit operator contractions, and matrix-free
C3v analysis of complete TDA root subspaces. **It does not yet automatically
identify excited states from a complete ABACUS LR calculation directory or
establish the NV⁻ ³E assignment.**

| Module | Implemented | Still needed |
|---|---|---|
| ABACUS wavefunction IO | Single-frame Γ real and multi-k complex LCAO text | Binary, appended frames, spinor semantics |
| ABACUS metadata | `OutputReader` / `AbacusReader` interfaces; `Orbital` AO labels | Logs, STRU, and calculation-version integration |
| Overlap IO | Legacy/modern text CSR S(R), single/multi-k S(k), NPY | Binary and other CSR dialects |
| LR eigenvector IO | `stateid_npz` v1: X/Y, energies, ph mapping | Native complete ABACUS LR files, MPI shard reconstruction |
| Symmetry | S-orthonormalization, D/χ, closure, C3v matching, Γ s/p/d AO operations, matrix-free TDA root projection | STRU integration and automatic symmetry discovery |
| Harmonic conventions | ABACUS m order, complex/real basis conversion, Γ s/p/d proper/improper operations and periodic atom mapping | Higher angular momentum and full Wigner D |
| Spin | Unrestricted-determinant ⟨S²⟩; TDA contraction with supplied S²_ph | Automatic S²_ph construction; full-LR spin response |

Unsupported formats and APIs explicitly raise `NotImplementedError` (or its IO
subclass `UnsupportedFormatError`). They do not guess matrix order or continue
with fabricated data.

Symmetry identification currently targets collinear spin without SOC, mainly
at Γ or in finite systems. IO and Fourier transforms support arbitrary k.
Non-Γ irrep analysis requires the little group of k and the corresponding
T(R,k), which are not constructed automatically. Magnetic groups, double groups,
antiunitary operations, and maps between k points remain unsupported.

## Installation and quick start

Requires Python ≥ 3.10, NumPy, and SciPy. Sparse matrices use SciPy; tests use
standard-library `unittest`. From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
python -m unittest discover -s tests -v
python examples/c3v_minimal.py
```

Expected example output:

```text
A1: characters=[1. 1. 1.], identified=A1
E: characters=[ 2. -1.  0.], identified=E
two parallel spins: <S^2>=2.0
```

You can also match numerical characters ordered by conjugacy class:

```python
from stateid.symmetry import match_characters

assert match_characters([1, 1, 1]).irrep == "A1"
assert match_characters([2, -1, 0]).irrep == "E"
assert not match_characters([2, -0.7, 0]).valid
```

These values are **representative characters or class averages**, rather than
sums over each class. Character matching alone checks compatibility; use
`analyze_c3v` to also validate subspace closure and group relations.

## Legacy/modern CSR and multiple k points

The two CSR formats use different header adapters and share the numerical
parser for values, column indices, and row pointers. `RealSpaceMatrix` is a
backend-neutral sparse X(R) container for S or H; it does not convert units.

```python
from stateid.io import AbacusReader, read_csr

sr = read_csr("path/to/sr_nao.csr")  # Also reads legacy data-SR-sparse_SPIN0.csr
kpoints = [[0, 0, 0], [0.25, 0, 0], [0.5, 0.25, 0]]  # Reciprocal fractional coordinates
Sk = sr.to_k(kpoints)               # (nk, nao, nao)
S_gamma = sr.to_k([0, 0, 0])        # (nao, nao)

reader = AbacusReader()
Sk_checked = reader.read_overlap("path/to/sr_nao.csr", format="abacus_csr",
                                 kpoints=kpoints)  # Hermitian/positive-definite checks per k
orbitals = reader.read_wavefunctions("path/to/wfs1k1_nao.txt", spin="alpha")
```

The cell gauge matches pyATB:

$$X(k)=\sum_R e^{+2\pi i k\cdot R}X(R).$$

There is no k weight or division by the number of R vectors, and individual
X(R) blocks need not be Hermitian. Files with multiple ionic frames require
an explicit `frame=0,1,...` (segment index, not the original step label).
Wavefunction-header `k_cartesian` values are Cartesian coordinates in
**2π/lat0 units**, not reciprocal fractional coordinates. Convert them with
`k_cartesian_to_fractional(k_cartesian, lattice_vectors)`.
See the [IO contract](docs/io-formats.md) and
[multi-k tutorial](docs/tutorials/04-csr-multik.md).

## Representations and characters

Orbitals are stored as columns:

$$|\psi_n\rangle=\sum_\mu|\phi_\mu\rangle C_{\mu n},\qquad
S_{\mu\nu}=\langle\phi_\mu|\phi_\nu\rangle.$$

For active AO coefficient transformations T, define

$$\hat R|\phi_\nu\rangle=\sum_\mu|\phi_\mu\rangle T_{\mu\nu}(R),\qquad
M_{\mu\nu}(R)=\langle\phi_\mu|\hat R|\phi_\nu\rangle=(ST)_{\mu\nu}.$$

If $C^\dagger SC=I$, the subspace representation and character are

$$D(R)=C^\dagger M(R)C=C^\dagger ST(R)C,\qquad
\chi(R)=\operatorname{Tr}D(R).$$

Choose `operator_kind="coefficient"` for T or `"matrix_element"` for M.
**This definition of M is not interchangeable with every object named M in
ABACUS source code.**

For an unnormalized subspace, form $G=C^\dagger SC$ and $Q=CG^{-1/2}$, then
use Q in place of C. A singular G requires a different subspace selection;
the library never silently drops orbitals. The AO coefficient projector is
$P=QQ^\dagger S$. Closure is measured by

$$\epsilon_R=\frac{\|T(R)Q-QD(R)\|_S}{\sqrt d},\qquad
\|A\|_S^2=\operatorname{Tr}(A^\dagger SA).$$

Only a complete invariant subspace with small residuals supports an exact
irrep assignment. Strain or Jahn–Teller distortions may remove C3v symmetry;
loosening tolerances does not restore it.

### Analyze a complete degenerate subspace

E is a two-dimensional irrep. A solver can return any unitary mixture of a
degenerate pair, $C'_E=C_EU$. Consequently,

$$D'(R)=U^\dagger D(R)U,\qquad \operatorname{Tr}D'(R)=\operatorname{Tr}D(R).$$

Individual orbital shapes, diagonal matrix elements, and x/y labels depend on
the selected basis. The complete subspace character does not. A single real E
component generally leaks into the other component under rotation.

| C3v | E (identity) | 2C3 | 3σv |
|---|---:|---:|---:|
| A1 | 1 | 1 | 1 |
| A2 | 1 | 1 | −1 |
| E | 2 | −1 | 0 |

Multiplicities are computed as

$$n_\alpha=\frac16\sum_c|c|\chi_\alpha(c)^*\chi(c).$$

The code checks nonnegative integer multiplicities and reconstruction error.
Reducible representations such as A1+E are reported without forcing a single
irrep label.

## Spin: interpreting ⟨S²⟩

For an integer-occupied, collinear unrestricted Slater determinant,

$$\frac{\langle\hat S^2\rangle}{\hbar^2}
=M_S^2+\frac{N_\alpha+N_\beta}{2}
-\sum_{ij}|(C_\alpha^\dagger SC_\beta)_{ij}|^2,\qquad
M_S=\frac{N_\alpha-N_\beta}{2}.$$

Pure singlets, doublets, and triplets have values 0, 3/4, and 2. The interface
requires **all occupied spatial orbitals** in both spin channels. A single KS
orbital or `nspin=2` cannot establish triplet character. For UKS, this is a
diagnostic of the auxiliary determinant, rather than a measurement of the exact
interacting state. An expectation value alone does not generally prove spin
purity.

`tda_s2` contracts a supplied complete operator in an orthonormal determinant
or CSF basis as $X^\dagger S^2_{ph}X/(X^\dagger X)$; `build_s2_ph` remains
unimplemented. Full-LR X/Y amplitudes have response-theory normalization and
cannot simply use an ordinary CI-vector formula. `lr_s2` explicitly refuses
that calculation. See the [spin tutorial](docs/tutorials/02-spin.md).

## NV⁻: an illustrative future workflow

The intended evidence chain combines an A1 orbital at ↓126, an E subspace at
↓127/128, a1→e LR transition components, reference-state symmetry, and spin
evidence to assess ³E. Orbital labels a1/e are returned as `A1`/`E`; the library
does not attach spin superscripts automatically.

```python
from stateid.io import AbacusReader
from stateid.symmetry import analyze_c3v

# Illustration only: supply real, matching files and verified AO operations.
reader = AbacusReader()
down = reader.read_wavefunctions("path/to/down_gamma.txt", spin="beta")
S = reader.read_overlap("path/to/sr_nao.csr", format="abacus_csr")
# operations: all six AO coefficient transformations T(R).
# Verify that 126/127/128 are one-based file band labels before selecting columns.
a1 = analyze_c3v(down.coefficients[:, [125]], operations, S,
                 operator_kind="coefficient")
e = analyze_c3v(down.coefficients[:, [126, 127]], operations, S,
                operator_kind="coefficient")
```

This is an illustration requiring supplied `operations`, not a runnable NV⁻
example. A runnable synthetic example is available in `examples/c3v_minimal.py`.
See the [NV⁻ tutorial](docs/tutorials/03-nv-minus.md) for the evidence boundaries.

## ABACUS TDA transition composition

`read_transition_analysis` reads thresholded `trans_analysis_*_tda.dat` files,
preserving original weights and spin/band labels. These are composition reports,
not complete X vectors or automatic irrep/spin assignments. Usage and real NV⁻
root-pair composition results are in the
[transition tutorial](docs/tutorials/05-transitions.md).

## Γ AO operations and TDA root subspaces

`read_ao_labels` reads ABACUS `Orbital`: its first column is a zero-based atom
index, and its m column is a real-harmonic component index converted to signed
m. Data-row order is the AO row order of C.

`build_gamma_ao_operation` uses explicit fractional positions, row-vector
lattice vectors, and active Cartesian operations to build Γ-point T for s/p/d
shells. It returns atom mappings, lattice return vectors, and geometric errors.
STRU parsing and automatic symmetry discovery remain unimplemented; validate
T†ST and orbital closure separately.

`analyze_tda_c3v` supports one or two complete spin-conserving ph product blocks
without materializing a ph² matrix. It derives the reference phase from
explicitly supplied complete alpha/beta occupied representations and returns
root orthonormality, reference phases, closure, group relations, and character
diagnostics. Thresholded composition tables cannot replace complete X. Full-LR
X/Y and automatic spin-multiplicity inference are unsupported. See the
[AO/TDA tutorial](docs/tutorials/06-ao-tda-symmetry.md).

## Roadmap and repository layout

1. Completed: numerical core, C3v, determinant spin, Γ/multi-k text, legacy/modern
   CSR, general Fourier transforms, array exchange, Γ s/p/d AO operations, and
   matrix-free TDA subspace analysis.
2. Next priority: collect C(k), S(R), structure, and AO metadata from the same
   calculation for a real joint $C(k)^\dagger S(k)C(k)\approx I$ regression.
3. Integrate STRU/orbital metadata and symmetry discovery; extend T(R,k).
4. Recover native ABACUS LR global ph indices, spin blocks, MPI shards, and X/Y
   conventions, starting with TDA.
5. Add real NV⁻ regressions, more point groups, and many-electron analysis;
   introduce other backends as the data model stabilizes.

```text
stateid/
├── README.md / README.zh-CN.md
├── LICENSE
├── pyproject.toml
├── src/stateid/
│   ├── io/          # Backend adapters and explicit data contracts
│   ├── realspace.py # Sparse X(R) and arbitrary-k Fourier transforms
│   ├── symmetry/    # Representations, C3v, harmonics, AO and TDA operations
│   └── spin/        # Determinant S², TDA contractions, LR placeholders
├── examples/c3v_minimal.py
├── tests/          # Mathematical invariants, error paths, and IO fixtures
└── docs/
    ├── abacus-conventions.md
    ├── io-formats.md
    ├── development/
    └── tutorials/  # Derivations, diagnostics, and physical evidence chains
```

Most detailed tutorials are currently in Chinese. Local ABACUS reference paths,
commits, and conventions are recorded in
[the convention notes](docs/abacus-conventions.md). Development has not modified
ABACUS source or copied its C++ implementation into this project. New backends
should implement `OutputReader` and keep file layouts out of the physics core.

## License

[MIT License](LICENSE), copyright © 2026 maki49.
