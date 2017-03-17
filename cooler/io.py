# -*- coding: utf-8 -*-
from __future__ import division, print_function
import posixpath
import warnings
import six

import numpy as np
import h5py

from ._reader import (ContactReader, HDF5Aggregator, TabixAggregator,
                      PairixAggregator, CoolerAggregator, SparseLoader, 
                      DenseLoader)
from ._writer import (write_chroms, write_bins, write_pixels, write_indexes,
                      write_info, MAGIC, URL)
from .util import get_binsize
from . import get_logger


logger = get_logger()


def ls(filepath):
    """
    Traverse file hierarchy and list all cooler nodes.

    Parameters
    ----------
    filepath : str

    Returns
    -------
    list of Cooler group paths in the file
    
    """
    listing = []
    keys = ['chroms', 'bins', 'pixels', 'indexes']
    def _check_group(pth, grp):
        fmt = grp.attrs.get('format', None)
        url = grp.attrs.get('format-url', None)
        if fmt == MAGIC or url == URL:
            if not all(name in grp.keys() for name in keys):
                warnings.warn(
                    'Cooler path /{} appears to be corrupt'.format(pth))
            listing.append('/' + pth)

    with h5py.File(filepath, 'r') as f:
        _check_group('/', f)
        f.visititems(_check_group)

    return listing


def is_cooler(filepath, group=None):
    """
    Determine if a file contains a valid Cooler data hierarchy.

    Parameters
    ----------
    filepath : str
    group : str, optional
        Path to specific group to check. Otherwise returns True
        if any Cooler paths are detected.

    Returns
    -------
    bool
    
    """
    if not h5py.is_hdf5(filepath):
        return False
    if group is None:
        return len(ls(filepath)) > 0
    if not group.startswith('/'):
        group = '/' + group
    return group in ls(filepath)


def create(filepath, chromsizes, bins, iterator, metadata=None, assembly=None,
           h5opts=None, group='/', append=False, lock=None):
    """
    Create a new Cooler file.

    Parameters
    ----------
    filepath : str
        Name of output HDF5 file. If the file does not exist, it will be 
        created.
    chromsizes : OrderedDict, pandas.Series or sequence of pairs
        Names of chromosomes or contigs ordered as they will appear in the
        contact matrix mapped to their lengths in base pairs.
    bins : pandas.DataFrame
        Segmentation of the chromosomes into genomic bins as a BED-like
        DataFrame with columns ``chrom``, ``start`` and ``end``.
    iterator : cooler.io.ContactReader
        An iterator that yields all the processed pixel data as sequence of 
        chunks.
    metadata : dict, optional
        Experiment metadata to store in the file. Must be JSON compatible.
    assembly : str, optional
        Name of genome assembly.
    h5opts : dict, optional
        HDF5 dataset filter options to use (compression, shuffling,
        checksumming, etc.). Default is to use autochunking and GZIP
        compression, level 6.
    group : str, optional
        Target path of the Cooler group. Default is the root group "/". If the
        group already exists it will be overwritten.
    lock : multiprocessing.Lock, optional
        Optional lock to synchronize concurrent HDF5 file access.

    Result
    ------
    Cooler hierarchy stored in ``filepath`` under ``group``.

    """
    mode = 'a' if append else 'w'
    if h5opts is None:
        h5opts = dict(compression='gzip', compression_opts=6)

    try:
        chromsizes = six.iteritems(chromsizes)
    except AttributeError:
        pass
    chroms, lengths = zip(*chromsizes)
    binsize = get_binsize(bins)
    n_chroms = len(chroms)
    n_bins = len(bins)

    with h5py.File(filepath, mode) as f:
        logger.info('Creating cooler at "{}::{}"'.format(filepath, group))
        if group is '/':
            for name in ['chroms', 'bins', 'pixels', 'indexes']:
                if name in f:
                    del f[name]
        else:
            try:
                f.create_group(group)
            except ValueError:
                del h5[group]
                f.create_group(group)

    with h5py.File(filepath, 'r+') as f:
        h5 = f[group]

        logger.info('Writing chroms')
        grp = h5.create_group('chroms')
        write_chroms(grp, chroms, lengths, h5opts)

        logger.info('Writing bins')
        grp = h5.create_group('bins')
        chrom_offset = write_bins(grp, chroms, bins, h5opts)
        h5.create_group('pixels')

    # Multiprocess HDF5 reading is supported only while absolutely no HDF5 file
    # is open in write mode anywhere.
    # If provided with a lock shared with HDF5-reading processes, `write_pixels` 
    # will acquire it and open the file for writing for the duration of each
    # write step only. After it closes the file and releases the lock, the 
    # reading processes will have to acquire the lock and re-open the file to 
    # obtain the correct file state for reading.
    logger.info('Writing pixels')
    target = posixpath.join(group, 'pixels')
    bin1_offset, nnz, ncontacts = write_pixels(
        filepath, target, n_bins, iterator, h5opts, lock)

    with h5py.File(filepath, 'r+') as f:
        h5 = f[group]

        logger.info('Writing indexes')
        grp = h5.create_group('indexes')
        write_indexes(grp, chrom_offset, bin1_offset, h5opts)

        logger.info('Writing info')
        info = {}
        info['bin-type'] = 'fixed' if binsize is not None else 'variable'
        info['bin-size'] = binsize if binsize is not None else 'null'
        info['nchroms'] = n_chroms
        info['nbins'] = n_bins
        info['sum'] = ncontacts
        info['nnz'] = nnz
        if assembly is not None:
            info['genome-assembly'] = assembly
        if metadata is not None:
            info['metadata'] = metadata
        write_info(h5, info)

    logger.info('Done')
