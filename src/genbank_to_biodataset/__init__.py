"""GenBank to BioDataset package metadata."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("genbank-to-biodataset")
except PackageNotFoundError:  # pragma: no cover - local dev
    __version__ = "0.0.0"

__all__ = ["__version__"]
