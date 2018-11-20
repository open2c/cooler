# -*- coding: utf-8 -*-
"""
Cooler
~~~~~~

A cool place to store your Hi-C.

:copyright: (c) 2016 Massachusetts Institute of Technology
:author: Nezar Abdennur
:license: BSD

"""
import logging
__version__ = '0.8.0-dev'
__format_version__ = 2
_loggers = {}


def get_logger(name='cooler'):
    # Based on ipython traitlets
    global _loggers

    if name not in _loggers:
        _loggers[name] = logging.getLogger(name)
        # Add a NullHandler to silence warnings about not being
        # initialized, per best practice for libraries.
        _loggers[name].addHandler(logging.NullHandler())

    return _loggers[name]


from .api import Cooler, get, info, chroms, bins, pixels, matrix, annotate
from .util import read_chromsizes, binnify, open_hdf5
from . import util
from . import ice
from . import io
