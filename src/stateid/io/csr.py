"""Read old/LTS and new ABACUS text CSR with one shared R-block parser."""

from pathlib import Path
import re

import numpy as np
from scipy.sparse import csr_matrix, vstack

from stateid.realspace import RealSpaceMatrix


class _Lines:
    """Streaming nonblank lines with one lookahead; keep comments for headers."""

    def __init__(self, stream):
        self.lines = (line.strip() for line in stream if line.strip())
        self.peek = next(self.lines, None)

    def pop(self):
        if self.peek is None:
            raise ValueError("truncated CSR file")
        line, self.peek = self.peek, next(self.lines, None)
        return line

    def skip_comments(self):
        while self.peek is not None and self.peek.startswith("#"):
            self.pop()

    def tokens(self, count, convert):
        result = []
        while len(result) < count:
            self.skip_comments()
            # Parenthesized complex values may contain whitespace.
            tokens = re.findall(r"\([^()]*\)|[^\s]+", self.pop().partition("#")[0])
            result.extend(convert(token) for token in tokens)
        if len(result) != count:
            raise ValueError("incorrect CSR block length")
        return np.asarray(result)


def _value(token):
    token = token.replace("D", "E").replace("d", "e")
    if token.startswith("(") and token.endswith(")"):
        real, imag = token[1:-1].split(",")
        return complex(float(real), float(imag))
    return complex(float(token))


def _match(pattern, line):
    match = re.fullmatch(pattern, line)
    if match is None:
        raise ValueError(f"unexpected CSR header: {line[:120]}")
    return match


def _header(lines):
    line = lines.pop()
    meta = {}
    if line.startswith("---"):
        meta["format"] = "abacus_new"
        meta["ionic_step"] = int(_match(r"---\s*Ionic Step\s+(\d+)\s*---", line)[1])
        _match(r"#.*", lines.pop())  # matrix description
        nspin = int(_match(r"(\d+)\s*#\s*number of spin directions", lines.pop())[1])
        meta["spin_index"] = int(_match(r"(\d+)\s*#\s*spin index", lines.pop())[1])
        if not 1 <= meta["spin_index"] <= nspin:
            raise ValueError("invalid CSR spin index")
        size = int(_match(r"(\d+)\s*#\s*number of localized basis", lines.pop())[1])
        nr = int(_match(r"(\d+)\s*#\s*number of Bravais lattice vector R", lines.pop())[1])
        cell = []
        while lines.peek is not None and "CSR Format" not in lines.peek:
            if lines.peek.startswith("---"):
                raise ValueError("missing CSR Format marker")
            cell.append(lines.pop())
        if lines.peek is None:
            raise ValueError("missing CSR Format marker")
        lines.pop()
        # The new writer embeds a structure block before its CSR banner.
        if len(cell) < 5:
            raise ValueError("missing lattice metadata in new CSR header")
        meta["lattice_constant_bohr"] = float(cell[1])
        lattice = np.array([[float(x) for x in row.split()] for row in cell[2:5]])
        if (lattice.shape != (3, 3) or not np.all(np.isfinite(lattice))
                or np.linalg.matrix_rank(lattice) != 3
                or not np.isfinite(meta["lattice_constant_bohr"]) or meta["lattice_constant_bohr"] <= 0):
            raise ValueError("invalid lattice metadata in new CSR header")
        meta["lattice_vectors"] = lattice
    else:
        meta["format"] = "abacus_legacy"
        if line.startswith("STEP:"):
            meta["ionic_step"] = int(_match(r"STEP:\s*(\d+)", line)[1])
            line = lines.pop()
        size = int(_match(r"Matrix Dimension of .+?:\s*(\d+)", line)[1])
        nr = int(_match(r"Matrix number of .+?:\s*(\d+)", lines.pop())[1])
    if size < 1:
        raise ValueError("CSR matrix dimension must be positive")
    return size, nr, meta


def _blocks(lines, size, nr):
    translations, blocks = [], []
    seen = set()
    for _ in range(nr):
        header = lines.tokens(4, int)
        r, nnz = tuple(header[:3]), int(header[3])
        if nnz < 0 or r in seen:
            raise ValueError("negative nnz or duplicate R block")
        seen.add(r)
        translations.append(r)
        if nnz:
            values = lines.tokens(nnz, _value)
            columns = lines.tokens(nnz, int)
            pointers = lines.tokens(size + 1, int)
        else:
            # Legacy empty blocks have no arrays. The new writer may include
            # labelled empty arrays plus size+1 zero row pointers.
            labelled_pointers = False
            while lines.peek is not None and lines.peek.startswith("#"):
                labelled_pointers |= "CSR row pointers" in lines.pop()
            values, columns = np.array([], complex), np.array([], int)
            pointers = lines.tokens(size + 1, int) if labelled_pointers else np.zeros(size + 1, int)
        if (not np.all(np.isfinite(values)) or pointers[0] != 0 or pointers[-1] != nnz
                or np.any(np.diff(pointers) < 0) or np.any(columns < 0) or np.any(columns >= size)):
            raise ValueError("invalid CSR values, column indices or row pointers")
        block = csr_matrix((values, columns, pointers), shape=(size, size))
        block.check_format(full_check=True)
        block.sum_duplicates()
        blocks.append(block.reshape((1, size * size)))
    packed = vstack(blocks, format="csr") if blocks else csr_matrix((0, size * size), dtype=complex)
    return np.asarray(translations, dtype=int).reshape(-1, 3), packed


def iter_csr(path):
    """Yield all text CSR sections in file order, preserving ionic step labels.

    Both real numbers and C++ (real,imag) complex values are supported.
    Binary output and unrelated CSR dialects are not supported.
    """
    with Path(path).open(encoding="utf-8") as stream:
        lines = _Lines(stream)
        lines.skip_comments()
        while lines.peek is not None:
            size, nr, meta = _header(lines)
            r, values = _blocks(lines, size, nr)
            yield RealSpaceMatrix(r, values, size, source=str(path), **meta)
            lines.skip_comments()


def read_csr(path, *, frame=None):
    """Read one CSR section. Multiple sections require an explicit frame index.

    frame is a nonnegative ZERO-BASED section index, not the raw ionic-step
    label (old/new ABACUS versions use different step numbering conventions).
    """
    if frame is not None and (isinstance(frame, bool) or not isinstance(frame, (int, np.integer)) or frame < 0):
        raise ValueError("frame must be a nonnegative integer index")
    result = None
    for index, data in enumerate(iter_csr(path)):
        if frame is None and index > 0:
            raise ValueError("multiple CSR sections: select frame=0, frame=1, ... explicitly")
        if index == (0 if frame is None else frame):
            result = data
    if result is None:
        raise ValueError("requested CSR frame does not exist")
    return result
