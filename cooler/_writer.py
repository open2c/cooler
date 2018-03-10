# -*- coding: utf-8 -*-
"""
HDF5 writers
~~~~~~~~~~~~

This module defines and implements the Cooler schema.

"""
from __future__ import division, print_function
from datetime import datetime
import json
import six

import numpy as np
import pandas as pd
import h5py

from . import __version__, __format_version__, get_logger
from .core import put
from .util import rlencode


logger = get_logger()


MAGIC = u"HDF5::Cooler"
URL = u"https://github.com/mirnylab/cooler"
CHROM_DTYPE = np.dtype('S')
CHROMID_DTYPE = np.int32
CHROMSIZE_DTYPE = np.int32
COORD_DTYPE = np.int32
BIN_DTYPE = np.int64
COUNT_DTYPE = np.int32
CHROMOFFSET_DTYPE = np.int64
BIN1OFFSET_DTYPE = np.int64
PIXEL_FIELDS = ('bin1_id', 'bin2_id', 'count')
PIXEL_DTYPES = (('bin1_id', BIN_DTYPE), 
                ('bin2_id', BIN_DTYPE), 
                ('count', COUNT_DTYPE))


def write_chroms(grp, chroms, h5opts):
    """
    Write the chromosome table.

    Parameters
    ----------
    grp : h5py.Group
        Group handle of an open HDF5 file with write permissions.
    chroms : DataFrame
        Chromosome table containing at least 'chrom' and 'length' columns
    h5opts : dict
        HDF5 dataset filter options.

    """
    n_chroms = len(chroms)
    names = np.array(chroms['name'], dtype=CHROM_DTYPE)  # auto-adjusts char length
    grp.create_dataset('name',
                       shape=(n_chroms,),
                       dtype=names.dtype,
                       data=names,
                       **h5opts)
    grp.create_dataset('length',
                       shape=(n_chroms,),
                       dtype=CHROMSIZE_DTYPE,
                       data=chroms['length'],
                       **h5opts)
    
    # Extra columns
    columns = list(chroms.keys())
    for col in ['name', 'length']:
        columns.remove(col)
    if columns:
        put(grp, chroms[columns])


def write_bins(grp, bins, chromnames, h5opts, chrom_as_enum=True):
    """
    Write the genomic bin table.

    Parameters
    ----------
    grp : h5py.Group
        Group handle of an open HDF5 file with write permissions.
    bins : pandas.DataFrame
        BED-like data frame with at least three columns: ``chrom``, ``start``,
        ``end``, sorted by ``chrom`` then ``start``, and forming a complete
        genome segmentation. The ``chrom`` column must be sorted according to
        the ordering in ``chroms``.
    chromnames : sequence of str
        Contig names.
    h5opts : dict
        HDF5 dataset filter options.

    """
    n_chroms = len(chromnames)
    n_bins = len(bins)
    idmap = dict(zip(chromnames, range(n_chroms)))

    # Convert chrom names to enum
    chrom_ids = [idmap[chrom] for chrom in bins['chrom']]
    if chrom_as_enum:
        chrom_dtype = h5py.special_dtype(enum=(CHROMID_DTYPE, idmap))
    else:
        chrom_dtype = CHROMID_DTYPE

    # Store bins
    try:
        chrom_dset = grp.create_dataset('chrom',
                           shape=(n_bins,),
                           dtype=chrom_dtype,
                           data=chrom_ids,
                           **h5opts)
    except ValueError:
        # If too many scaffolds for HDF5 enum header, 
        # try storing chrom IDs as raw int instead
        if chrom_as_enum:
            chrom_as_enum = False
            chrom_dtype = CHROMID_DTYPE
            chrom_dset = grp.create_dataset('chrom',
                               shape=(n_bins,),
                               dtype=chrom_dtype,
                               data=chrom_ids,
                               **h5opts)
        else:
            raise
    if not chrom_as_enum:
        chrom_dset.attrs['enum_path'] = u'/chroms/name'

    grp.create_dataset('start',
                       shape=(n_bins,),
                       dtype=COORD_DTYPE,
                       data=bins['start'],
                       **h5opts)
    grp.create_dataset('end',
                       shape=(n_bins,),
                       dtype=COORD_DTYPE,
                       data=bins['end'],
                       **h5opts)
    
    # Extra columns
    columns = list(bins.keys())
    for col in ['chrom', 'start', 'end']:
        columns.remove(col)
    if columns:
        put(grp, bins[columns])


def prepare_pixels(grp, n_bins, columns, dtypes, h5opts):
    columns = list(columns)
    max_size = n_bins * (n_bins - 1) // 2 + n_bins
    init_size = min(5 * n_bins, max_size)
    grp.create_dataset('bin1_id',
                       dtype=dtypes.get('bin1_id', BIN_DTYPE),
                       shape=(init_size,),
                       maxshape=(max_size,),
                       **h5opts)
    grp.create_dataset('bin2_id',
                       dtype=dtypes.get('bin2_id', BIN_DTYPE),
                       shape=(init_size,),
                       maxshape=(max_size,),
                       **h5opts)
    grp.create_dataset('count',
                       dtype=dtypes.get('count', COUNT_DTYPE),
                       shape=(init_size,),
                       maxshape=(max_size,),
                       **h5opts)
    for col in ['bin1_id', 'bin2_id', 'count']:
        columns.remove(col)
    if columns:
        for col in columns:
            grp.create_dataset(col,
                               dtype=dtypes.get(col, float),
                               shape=(init_size,),
                               maxshape=(max_size,),
                               **h5opts)


def write_pixels(filepath, grouppath, columns, iterable, h5opts, lock):
    """
    Write the non-zero pixel table.

    Parameters
    ----------
    filepath : str
        Path to HDF5 output file.
    grouppath : str
        Qualified path to destination HDF5 group.
    columns : sequence
        Sequence of column names
    iterable : an iterable object
        An object that processes and yields binned contacts from some input 
        source as a stream of chunks. The chunks must be either pandas 
        DataFrames or mappings of column names to arrays.
    h5opts : dict
        HDF5 filter options.
    lock : multiprocessing.Lock, optional
        Optional lock to synchronize concurrent HDF5 file access.

    """
    nnz = 0
    total = 0
    for i, chunk in enumerate(iterable):
        
        if isinstance(chunk, pd.DataFrame):
            chunk = {k: v.values for k, v in six.iteritems(chunk)}
        
        try:
            if lock is not None:
                lock.acquire()

            logger.debug("writing chunk {}".format(i))
            
            with h5py.File(filepath, 'r+') as fw:
                grp = fw[grouppath]
                dsets = [grp[col] for col in columns]

                n = len(chunk[columns[0]])
                for col, dset in zip(columns, dsets):
                    dset.resize((nnz + n,))
                    dset[nnz:nnz+n] = chunk[col]
                nnz += n
                total += chunk['count'].sum()
                
                fw.flush()

        finally:
            if lock is not None:
                lock.release()

    return nnz, total


def index_pixels(grp, n_bins, nnz):
    bin1 = grp['bin1_id']
    bin1_offset = np.zeros(n_bins + 1, dtype=BIN1OFFSET_DTYPE)
    curr_val = 0
    for start, length, value in zip(*rlencode(bin1, 1000000)):
        bin1_offset[curr_val:value + 1] = start
        curr_val = value + 1
    bin1_offset[curr_val:] = nnz
    return bin1_offset


def index_bins(grp, n_chroms, n_bins):
    chrom_ids = grp['chrom']
    chrom_offset = np.zeros(n_chroms + 1, dtype=CHROMOFFSET_DTYPE)
    curr_val = 0
    for start, length, value in zip(*rlencode(chrom_ids)):
        chrom_offset[curr_val:value + 1] = start
        curr_val = value + 1
    chrom_offset[curr_val:] = n_bins
    return chrom_offset


def write_indexes(grp, chrom_offset, bin1_offset, h5opts):
    """
    Write the indexes.

    Parameters
    ----------
    grp : h5py.Group
        Group handle of an open HDF5 file with write permissions.
    chrom_offset : sequence
        Lookup table: chromosome ID -> first row in bin table (bin ID)
        corresponding to that chromosome.
    bin1_offset : sequence
        Lookup table: genomic bin ID -> first row in pixel table (pixel ID)
        having that bin on the first axis.

    """
    grp.create_dataset(
        "chrom_offset",
        shape=(len(chrom_offset),), 
        dtype=CHROMOFFSET_DTYPE,
        data=chrom_offset, 
        **h5opts)
    grp.create_dataset(
        "bin1_offset",
        shape=(len(bin1_offset),), 
        dtype=BIN1OFFSET_DTYPE,
        data=bin1_offset, 
        **h5opts)


def write_info(grp, info):
    """
    Write the file description and metadata attributes.

    Parameters
    ----------
    grp : h5py.Group
        Group handle of an open HDF5 file with write permissions.
    info : dict
        Dictionary, unnested with the possible exception of the ``metadata``
        key. ``metadata``, if present, must be JSON-serializable.

    Required keys
    -------------
    nbins : int
        number of genomic bins
    nnz : int
        number of non-zero pixels

    """
    assert 'nbins' in info
    assert 'nnz' in info
    info.setdefault('genome-assembly', 'unknown')
    info['metadata'] = json.dumps(info.get('metadata', {}))
    info['creation-date'] = datetime.now().isoformat()
    info['generated-by'] = 'cooler-' + __version__
    info['format'] = MAGIC
    info['format-version'] = __format_version__
    info['format-url'] = URL
    grp.attrs.update(info)
