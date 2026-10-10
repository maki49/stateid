"""Explicit errors for formats outside an adapter's documented scope."""


class UnsupportedFormatError(NotImplementedError):
    """The input requires an adapter that has not yet been implemented."""
