from .abacus import AbacusReader, UnsupportedFormatError, read_gamma_wavefunctions, read_wavefunctions
from .arrays import read_lr_npz, read_overlap_npy
from .csr import iter_csr, read_csr
from .models import LRData, OrbitalData, OutputReader

__all__ = ["AbacusReader", "UnsupportedFormatError", "read_gamma_wavefunctions", "read_wavefunctions", "iter_csr", "read_csr",
           "read_lr_npz", "read_overlap_npy", "LRData", "OrbitalData", "OutputReader"]

__all__ += ["read_transition_analysis", "Transition", "TransitionState"]

from .ao_labels import read_ao_labels
__all__ += ["read_ao_labels"]


def __getattr__(name):
    # Keep the CLI module unloaded for python -m stateid.io.abacus_transitions.
    if name in {"read_transition_analysis", "Transition", "TransitionState"}:
        from importlib import import_module
        return getattr(import_module(".abacus_transitions", __name__), name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
