# -*- coding: utf-8 -*-
"""
Cooler
~~~~~~

A cool place to store your Hi-C.

:copyright: (c) 2016 Massachusetts Institute of Technology
:author: Nezar Abdennur
:license: BSD

"""
__version__ = '0.4.0'
__format_version__ = 2

from .api import Cooler, get, info, chroms, bins, pixels, matrix, annotate
from .util import read_chromsizes, binnify
from .io import open_hdf5
from . import util
from . import ice
from . import io
