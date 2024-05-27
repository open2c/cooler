"""
Cooler
~~~~~~

A cool place to store your Hi-C.

:author: Nezar Abdennur
:license: BSD-3-Clause

"""
from . import fileops, parallel
from ._balance import balance_cooler
from ._logging import get_verbosity_level, set_verbosity_level
from ._reduce import coarsen_cooler, merge_coolers, zoomify_cooler
from ._version import __format_version__, __version__
from .api import Cooler, annotate
from .create import create_cooler, create_scool, rename_chroms
from .util import binnify, fetch_chromsizes, read_chromsizes

__all__ = [
    "Cooler",
    "annotate",
    "balance_cooler",
    "create_cooler",
    "create_scool",
    "rename_chroms",
    "coarsen_cooler",
    "merge_coolers",
    "zoomify_cooler",
    "fileops",
    "parallel",
    "binnify",
    "fetch_chromsizes",
    "read_chromsizes",
    "get_verbosity_level",
    "set_verbosity_level",
    "__version__",
    "__format_version__",
]
