# -*- coding: utf-8 -*-
"""
Cooler
~~~~~~

A cool place to store your Hi-C.

:copyright: (c) 2016 Massachusetts Institute of Technology
:author: Nezar Abdennur
:license: BSD

"""
from ._version import __version__, __format_version__
from .api import Cooler, annotate
from .create import create_cooler, rename_chroms
from .reduce import merge_coolers, coarsen_cooler, zoomify_cooler
from .balance import balance_cooler
from .util import binnify, read_chromsizes, fetch_chromsizes
from . import tools
from . import fileops

from . import create
from . import balance
from . import io  # deprecated module
ice = balance  # alias
