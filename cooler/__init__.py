"""
Cooler
~~~~~~

A cool place to store your Hi-C.

:copyright: (c) 2016 Massachusetts Institute of Technology
:author: Nezar Abdennur
:license: BSD

"""
from ._version import __version__, __format_version__
from ._logging import get_verbosity_level, set_verbosity_level
from .api import Cooler, annotate
from .create import create_cooler, rename_chroms, create_scool
from .reduce import merge_coolers, coarsen_cooler, zoomify_cooler
from .balance import balance_cooler
from .util import binnify, read_chromsizes, fetch_chromsizes
from . import parallel
from . import fileops

from . import create
from . import balance

ice = balance  # alias
