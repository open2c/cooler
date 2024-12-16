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
    "__format_version__",
    "__version__",
    "annotate",
    "balance_cooler",
    "binnify",
    "coarsen_cooler",
    "create_cooler",
    "create_scool",
    "fetch_chromsizes",
    "fileops",
    "get_verbosity_level",
    "merge_coolers",
    "parallel",
    "read_chromsizes",
    "rename_chroms",
    "set_verbosity_level",
    "zoomify_cooler",
]
