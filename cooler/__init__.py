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
__version__ = '0.7.7'
__format_version__ = 2
_logger = None


def get_logger():
    # Based on ipython traitlets
    global _logger

    if _logger is None:
        _logger = logging.getLogger('cooler')
        # Add a NullHandler to silence warnings about not being
        # initialized, per best practice for libraries.
        _logger.addHandler(logging.NullHandler())

    return _logger


from .api import Cooler, get, info, chroms, bins, pixels, matrix, annotate
from .util import read_chromsizes, binnify, open_hdf5
from . import util
from . import ice
from . import io
