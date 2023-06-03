"""
Cooler
~~~~~~

A cool place to store your Hi-C.

:copyright: (c) 2016 Massachusetts Institute of Technology
:author: Nezar Abdennur
:license: BSD

"""
from . import balance, create, fileops, parallel, tools
from ._logging import get_verbosity_level, set_verbosity_level
from ._version import __format_version__, __version__
from .api import Cooler, annotate
from .balance import balance_cooler
from .create import create_cooler, create_scool, rename_chroms
from .reduce import coarsen_cooler, merge_coolers, zoomify_cooler
from .util import binnify, fetch_chromsizes, read_chromsizes

ice = balance  # alias
