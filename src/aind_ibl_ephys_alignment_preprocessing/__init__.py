"""Library to prepare histology and ephys for the IBL ephys alignment GUI"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("aind-ibl-ephys-alignment-preprocessing")
except PackageNotFoundError:
    __version__ = "0.0.0.dev0"
