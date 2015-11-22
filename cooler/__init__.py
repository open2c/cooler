# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals
from .format import (from_dense, get_slice, get_matrix, 
                     get_scaffolds, get_bins, get_metadata)
from .util import open_hdf5
from . import genome

__version__ = '0.1dev'
__all__ = ['from_dense', 'get_slice', 'get_matrix', 'get_scaffolds', 
           'get_bins', 'get_metadata', 'open_hdf5', 'genome']
