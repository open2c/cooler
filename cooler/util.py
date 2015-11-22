# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals
from contextlib import contextmanager
import gzip
import six
import re

import numpy as np
import pandas
import h5py


@contextmanager
def open_hdf5(fp, mode='r', *args, **kwargs):
    """
    Context manager like ``h5py.File`` but accepts already open HDF5 file 
    handles which do not get closed on teardown.

    Parameters
    ----------
    fp : str or ``h5py.File`` object
        If an open file object is provided, it passes through unchanged, 
        provided that the requested mode is compatible.
        If a filepath is passed, the context manager will close the file on
        tear down.

    mode : str
        r        Readonly, file must exist
        r+       Read/write, file must exist
        a        Read/write if exists, create otherwise
        w        Truncate if exists, create otherwise
        w- or x  Fail if exists, create otherwise

    """
    if isinstance(fp, six.string_types):
        own_fh = True
        fh = h5py.File(fp, mode, *args, **kwargs)
    else:
        own_fh = False
        if mode =='r' and fp.mode == 'r+':
            raise ValueError("File object provided is not in readonly mode")
        elif mode in ('r+', 'a', 'w') and fp.mode == 'r':
            raise ValueError("File object provided is not writeable")
        else:
            raise ValueError("File exists")
        fh = fp
    try:
        yield fh
    finally:
        if own_fh:
            fh.close()


def atoi(s):
    return int(s.replace(',', ''))


def to_tsv(fp, table, index=False, **kwargs):
    kwargs['index'] = index
    if isinstance(fp, six.string_types) and fp.endswith('.gz'):
        with gzip.open(fp, 'wb') as f:
            table.to_csv(f, sep=b'\t', **kwargs)
    else:
        table.to_csv(fp, sep=b'\t', **kwargs)



def read_tsv(fp, **kwargs):
    if isinstance(fp, six.string_types) and fp.endswith('.gz'):
        kwargs['compression'] = 'gzip'
    return pandas.read_csv(fp, sep=b'\t', **kwargs)


_NS_REGEX = re.compile(r'(\d+)', re.U)


def natsort_key(s):
    return tuple([int(x) if x.isdigit() else x for x in _NS_REGEX.split(s) if x])


def natsorted(iterable):
    return sorted(iterable, key=natsort_key)


# def validate(h5):
#     pass


# class CoolerValidator(object):
#     FORMAT_URL = 'https://bitbucket.org/nvictus/cooler'
#     FORMAT_VERSIONS = set([(0, 1)])

#     ATTRS = {(0, 1): (
#         'id', 'generated-by', 'creation-date',
#         'format-url', 'format-version',
#         'bin-type', 'bin-size', 
#         'genome-assembly',
#         'shape', 'nnz')}
#     GROUPS = {(0, 1): ('contigs', 'bins', 'matrix', 'indexes')}
