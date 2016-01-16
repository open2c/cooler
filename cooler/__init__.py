# -*- coding: utf-8 -*-
"""
Cooler
~~~~~~

A cool place to store your Hi-C.

:copyright: (c) 2016 Massachusetts Institute of Technology
:author: Nezar Abdennur
:license: MIT

"""
__version__ = '0.2dev'
__format_version__ = 0

from .api import Cooler, get, info, chromtable, bintable, pixeltable, matrix
from .util import read_chrominfo, make_bintable
from .io import open_hdf5
#from . import balancing
from . import util
from . import io


__all__ = [
    'Cooler',
    'get',
    'info',
    'chromtable',
    'bintable',
    'pixeltable',
    'matrix',
    'read_chrominfo',
    'make_bintable',
    'open_hdf5',
    'format',
    'balancing',
    'util',
]
