# -*- coding: utf-8 -*-
from __future__ import division, print_function
from contextlib import contextmanager
import warnings
import sys
import six

import numpy as np
import h5py

from ._reader import (ContactReader, HDF5Aggregator, TabixAggregator, CoolerAggregator,
                      SparseLoader, DenseLoader)
from ._writer import write_chroms, write_bins, write_pixels, write_indexes, write_info
from ..util import get_binsize


def create(h5, chroms, lengths, bins, reader, metadata=None, assembly=None, h5opts=None):
    """
    Create a new Cooler file.

    Parameters
    ----------
    h5 : h5py.File or h5py.Group
        Open HDF5 group handle which will serve as the root of the cooler
        hierarchy.
    chroms : sequence of str
        Names of chromosomes or contigs ordered as they will appear in the
        contact matrix.
    lengths : sequence of int
        Lengths of chromosomes of contigs in base pairs.
    bins : pandas.DataFrame
        Segmentation of the chromosomes into genomic bins as a BED-like
        DataFrame with columns ``chrom``, ``start`` and ``end``.
    reader : cooler.io.ContactReader
        bla
    metadata : dict, optional
        Experiment metadata to store in the file. Must be JSON compatible.
    assembly : str, optional
        Name of genome assembly.
    h5opts : dict, optional
        HDF5 dataset filter options to use (compression, shuffling,
        checksumming, etc.). Default is to use autochunking and GZIP
        compression, level 6.

    Result
    ------
    Cooler hierarchy stored under ``h5``.

    """
    if h5opts is None:
        h5opts = dict(compression='gzip', compression_opts=6)
    n_chroms = len(chroms)
    n_bins = len(bins)

    print('chroms')
    grp = h5.create_group('chroms')
    write_chroms(grp, chroms, lengths, h5opts)

    print('bins')
    binsize = get_binsize(bins)
    grp = h5.create_group('bins')
    chrom_offset = write_bins(grp, chroms, bins, h5opts)

    print('pixels')
    grp = h5.create_group('pixels')
    bin1_offset, nnz = write_pixels(grp, n_bins, reader, h5opts)

    print('indexes')
    grp = h5.create_group('indexes')
    write_indexes(grp, chrom_offset, bin1_offset, h5opts)

    print('info')
    info = {}
    info['bin-type'] = 'fixed' if binsize is not None else 'variable'
    info['bin-size'] = binsize if binsize is not None else 'null'
    info['nchroms'] = n_chroms
    info['nbins'] = n_bins
    info['nnz'] = nnz
    if assembly is not None:
        info['genome-assembly'] = assembly
    if metadata is not None:
        info['metadata'] = metadata
    write_info(h5, info)


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
        * r        Readonly, file must exist
        * r+       Read/write, file must exist
        * a        Read/write if exists, create otherwise
        * w        Truncate if exists, create otherwise
        * w- or x  Fail if exists, create otherwise

    """
    if isinstance(fp, six.string_types):
        own_fh = True
        fh = h5py.File(fp, mode, *args, **kwargs)
    else:
        own_fh = False
        if mode == 'r' and fp.file.mode == 'r+':
            #warnings.warn("File object provided is writeable but intent is read-only")
            pass
        elif mode in ('r+', 'a') and fp.file.mode == 'r':
            raise ValueError("File object provided is not writeable")
        elif mode == 'w':
            raise ValueError("Cannot truncate open file")
        elif mode in ('w-', 'x'):
            raise ValueError("File exists")
        fh = fp
    try:
        yield fh
    finally:
        if own_fh:
            fh.close()
