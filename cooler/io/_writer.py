# -*- coding: utf-8 -*-
"""
HDF5 writers
~~~~~~~~~~~~

This module defines and implements the Cooler schema.

"""
from __future__ import division, print_function
from datetime import datetime
import json

import numpy as np
import pandas
import h5py

from .. import __version__, __format_version__
from ..util import rlencode


CHROM_DTYPE = np.dtype('S32')
CHROMID_DTYPE = np.int32
CHROMSIZE_DTYPE = np.int32
COORD_DTYPE = np.int32
BIN_DTYPE = np.int64
COUNT_DTYPE = np.int32
CHROMOFFSET_DTYPE = np.int64
BIN1OFFSET_DTYPE = np.int64


def write_chroms(grp, chroms, lengths, h5opts):
    """
    Write the chromosome table.

    Parameters
    ----------
    grp : h5py.Group
        Group handle of an open HDF5 file with write permissions.
    chroms : sequence of str
        Contig names.
    lengths : sequence of int
        Contig lengths in base pairs.
    h5opts : dict
        HDF5 dataset filter options.

    """
    n_chroms = len(chroms)
    names = np.array(chroms, dtype=CHROM_DTYPE)
    grp.create_dataset('name',
                       shape=(n_chroms,),
                       dtype=CHROM_DTYPE,
                       data=names,
                       **h5opts)
    grp.create_dataset('length',
                       shape=(n_chroms,),
                       dtype=CHROMSIZE_DTYPE,
                       data=lengths,
                       **h5opts)


def write_bins(grp, chroms, bins, h5opts):
    """
    Write the genomic bin table.

    Parameters
    ----------
    grp : h5py.Group
        Group handle of an open HDF5 file with write permissions.
    chroms : sequence of str
        Contig names.
    bins : pandas.DataFrame
        BED-like data frame with at least three columns: ``chrom``, ``start``,
        ``end``, sorted by ``chrom`` then ``start``, and forming a complete
        genome segmentation. The ``chrom`` column must be sorted according to
        the ordering in ``chroms``.
    h5opts : dict
        HDF5 dataset filter options.

    """
    n_chroms = len(chroms)
    n_bins = len(bins)
    idmap = dict(zip(chroms, range(n_chroms)))

    # Convert chrom names to enum
    chrom_ids = [idmap[chrom] for chrom in bins['chrom']]
    enum_dtype = h5py.special_dtype(enum=(CHROMID_DTYPE, idmap))

    # Store bins
    grp.create_dataset('chrom',
                       shape=(n_bins,),
                       dtype=enum_dtype,
                       data=chrom_ids,
                       **h5opts)
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

    # Index the chromosome offsets
    chrom_offset = np.zeros(n_chroms + 1, dtype=CHROMOFFSET_DTYPE)
    curr_val = 0
    for start, length, value in zip(*rlencode(chrom_ids)):
        chrom_offset[curr_val:value + 1] = start
        curr_val = value + 1
    chrom_offset[curr_val:] = n_bins

    return chrom_offset


def write_pixels(grp, n_bins, reader, h5opts):
    """
    Write the non-zero pixel table.

    Parameters
    ----------
    grp : h5py.Group
        Group handle of an open HDF5 file with write permissions.
    n_bins : int
        Number of genomic bins.
    reader : object
        Reader object that reads and/or aggregates contacts from
        the input file(s). A reader returns chunks of binned contacts (bin1_id,
        bin2_id, count) sorted by ``bin1_id`` then ``bin2_id``.
    h5opts : dict
        HDF5 filter options.

    """
    n_pairs = reader.size()
    init_size = 5 * n_bins
    max_size = min(n_pairs, n_bins * (n_bins - 1) // 2 + n_bins)

    # Preallocate
    bin1 = grp.create_dataset('bin1_id',
                              dtype=BIN_DTYPE,
                              shape=(init_size,),
                              maxshape=(max_size,),
                              **h5opts)
    bin2 = grp.create_dataset('bin2_id',
                              dtype=BIN_DTYPE,
                              shape=(init_size,),
                              maxshape=(max_size,),
                              **h5opts)
    count = grp.create_dataset('count',
                              dtype=COUNT_DTYPE,
                              shape=(init_size,),
                              maxshape=(max_size,),
                              **h5opts)

    # Store the pixels
    nnz = 0
    for chunk in reader:
        n = len(chunk['bin1_id'])
        for dset in [bin1, bin2, count]:
            dset.resize((nnz + n,))
        bin1[nnz:nnz+n] = chunk['bin1_id']
        bin2[nnz:nnz+n] = chunk['bin2_id']
        count[nnz:nnz+n] = chunk['count']
        nnz += n

    # Index the first axis (matrix row) offsets
    bin1_offset = np.zeros(n_bins + 1, dtype=BIN1OFFSET_DTYPE)
    curr_val = 0
    for start, length, value in zip(*rlencode(bin1, 1000000)):
        bin1_offset[curr_val:value + 1] = start
        curr_val = value + 1
    bin1_offset[curr_val:] = nnz

    return bin1_offset, nnz


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
    grp.create_dataset("chrom_offset",
                       shape=(len(chrom_offset),), dtype=np.int32,
                       data=chrom_offset, **h5opts)
    grp.create_dataset("bin1_offset",
                       shape=(len(bin1_offset),), dtype=np.int32,
                       data=bin1_offset, **h5opts)


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
    info.setdefault('id', "null")
    info.setdefault('genome-assembly', 'unknown')
    info['metadata'] = json.dumps(info.get('metadata', {}))
    info['creation-date'] = datetime.now().isoformat()
    info['library-version'] = __version__
    info['format-version'] = __format_version__
    info['format-url'] = "https://github.com/mirnylab/cooler"
    grp.attrs.update(info)
