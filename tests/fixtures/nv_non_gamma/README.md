# NV⁻ non-Gamma regression data

Local ABACUS PBE SCF, 2026-10-07, 7 atoms, 91 AOs, scalar nspin=2.
Source run and SHA256 hashes are recorded in provenance.json.
See [tutorial 08](../../../docs/tutorials/08-real-abacus-nv.md) for physics, inputs,
units, residual definitions, and limitations.

All files are unmodified bytes from the successful calculation.
srg1_nao.csr.gz stores the complete 93-R-block S(R), not a fitted S(k).
wfk{1,2,3}s2g1_nao.txt.gz stores all 91 beta bands at three distinct k points.
gzip is only fixture packaging (deterministic mtime=0); tests decompress to
temporary files before invoking the ordinary text parsers. No ABACUS binary,
third-party source, pseudopotential, or orbital radial data is redistributed.

The complete 16 spin/k files and both H(R) matrices remain in the source run.
The integration example validates all of them; lightweight repository regression
uses the selected files above. provenance.json records the full-run summary,
the selected C3v results, and SHA256 of each uncompressed fixture.
