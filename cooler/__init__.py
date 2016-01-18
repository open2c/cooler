# -*- coding: utf-8 -*-
"""
Cooler
~~~~~~

A cool place to store your Hi-C.

:copyright: (c) 2016 Massachusetts Institute of Technology
:author: Nezar Abdennur
:license: BSD

"""
__version__ = '0.2'
__format_version__ = 0

from .api import Cooler, get, info, chromtable, bintable, pixeltable, matrix
from .util import read_chrominfo, make_bintable
from .io import open_hdf5
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
    'util',
    'io'
]
