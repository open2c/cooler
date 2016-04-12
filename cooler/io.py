# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals
from contextlib import contextmanager
from collections import OrderedDict
from datetime import datetime
import json
import six

from pandas.algos import is_lexsorted
import numpy as np
import pandas
import h5py

from . import __version__, __format_version__
from .util import lexbisect


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
        if mode == 'r' and fp.mode == 'r+':
            raise ValueError("File object provided is not in readonly mode")
        elif mode in ('r+', 'a', 'w') and fp.mode == 'r':
            raise ValueError("File object provided is not writeable")
        elif mode != 'r':
            raise ValueError("File exists")
        fh = fp
    try:
        yield fh
    finally:
        if own_fh:
            fh.close()


def write_chromtable(grp, chromtable, h5opts):
    n_chroms = len(chromtable)
    grp.create_dataset('name',
                       shape=(n_chroms,),
                       dtype='S32',
                       data=np.array(chromtable['name'], dtype='S32'),
                       **h5opts)
    grp.create_dataset('length',
                       shape=(n_chroms,),
                       dtype=np.int64,
                       data=chromtable['length'],
                       **h5opts)


def write_bintable(grp, chromtable, bintable, h5opts):
    n_chroms = len(chromtable)
    n_bins = len(bintable)
    if 'id' not in chromtable.columns:
        chromtable['id'] = np.arange(n_chroms)
    chrom_ids = chromtable['id'].loc[bintable['chrom']]
    grp.create_dataset('chrom_id',
                       shape=(n_bins,),
                       dtype=np.int32,
                       data=chrom_ids,
                       **h5opts)
    grp.create_dataset('start',
                       shape=(n_bins,),
                       dtype=np.int64,
                       data=bintable['start'],
                       **h5opts)
    grp.create_dataset('end',
                       shape=(n_bins,),
                       dtype=np.int64,
                       data=bintable['end'],
                       **h5opts)


def write_indexes(grp, chrom_offset, bin1_offset, h5opts):
    grp.create_dataset("chrom_offset",
                       shape=(len(chrom_offset),), dtype=np.int32,
                       data=chrom_offset, **h5opts)
    grp.create_dataset("bin1_offset",
                       shape=(len(bin1_offset),), dtype=np.int32,
                       data=bin1_offset, **h5opts)


def write_info(h5, info):
    info.setdefault('id', "null")
    info.setdefault('genome-assembly', 'unknown')
    info['metadata'] = json.dumps(info.get('metadata', {}))
    info['creation-date'] = datetime.now().isoformat()
    info['library-version'] = __version__
    info['format-version'] = __format_version__
    info['format-url'] = "https://github.com/mirnylab/cooler"
    h5.attrs.update(info)


def _aggregate(grp, chromtable, bintable, h5read, binsize, h5opts, chunksize,
               check_sorted):
    # TODO: do a first sweep to do the lexsort check, or else chunk selection
    # can potentially hang on unsorted data
    def _load_chunk(h5read, lo, hi):
        data = OrderedDict([
            ('chrom_id1', h5read['chrms1'][lo:hi]),
            ('cut1', h5read['cuts1'][lo:hi]),
            ('chrom_id2', h5read['chrms2'][lo:hi]),
            ('cut2', h5read['cuts2'][lo:hi]),
        ])
        if check_sorted and not is_lexsorted(
              [arr.astype(np.int64) for arr in data.values()]):
              raise ValueError("Paired read coordinates are not lexically sorted.")
        return pandas.DataFrame(data)

    n_reads = len(h5read['chrms1'])
    n_bins = len(bintable)
    chromtable = chromtable.set_index('name', drop=False)

    # convert genomic coords of bin starts to absolute
    accum_length = np.r_[0, np.cumsum(chromtable['length'])]
    bin_chrom_ids = chromtable['id'][bintable['chrom']].values
    abs_start_coords = accum_length[bin_chrom_ids] + bintable['start']

    # chrom offset index: chrom_id -> offset in bins
    chrom_nbins =  bintable.groupby(bin_chrom_ids, sort=False).size()
    chrom_offset = np.r_[0, np.cumsum(chrom_nbins)]

    # initialize the bin1 offset index
    bin1_offset = np.zeros(n_bins + 1, dtype=np.int64)

    # fill the pixel table
    init_size = 5 * n_bins
    max_size = min(n_reads, n_bins * (n_bins - 1) // 2 + n_bins)

    Bin1  = grp.create_dataset('bin1_id',
                               dtype=np.int64,
                               shape=(init_size,),
                               maxshape=(max_size,),
                               **h5opts)
    Bin2  = grp.create_dataset('bin2_id',
                               dtype=np.int64,
                               shape=(init_size,),
                               maxshape=(max_size,),
                               **h5opts)
    Count = grp.create_dataset('count',
                               dtype=np.int64,
                               shape=(init_size,),
                               maxshape=(max_size,),
                               **h5opts)

    hi = np.searchsorted(h5read['chrms1'], -1, 'right')
    bin1hi = 0
    nnz = 0 

    while hi < n_reads:
        # fetch next chunk, making sure our selection doesn't split a bin1
        lo, hi = hi, min(hi + chunksize, n_reads)
        cid, cut = h5read['chrms1'][hi-1], h5read['cuts1'][hi-1]
        abs_pos = accum_length[cid] + cut
        i = int(np.searchsorted(abs_start_coords, abs_pos, side='right')) - 1
        end_pos = bintable['end'][i]
        hi = lexbisect(
            (h5read["chrms1"], h5read["cuts1"]),
            (cid, end_pos),
            'left')
        print(lo, hi)
        
        # assign bins to reads
        table = _load_chunk(h5read, lo, hi)
        abs_pos1 = accum_length[h5read['chrms1'][lo:hi]] + h5read['cuts1'][lo:hi]
        abs_pos2 = accum_length[h5read['chrms2'][lo:hi]] + h5read['cuts2'][lo:hi]
        if check_sorted and np.any(abs_pos1 > abs_pos2):
            raise ValueError("Found a read pair in descending order.")

        if binsize is None:
            table['bin1'] = np.searchsorted(abs_start_coords, abs_pos1, side='right') - 1
            table['bin2'] = np.searchsorted(abs_start_coords, abs_pos2, side='right') - 1
        else:
            rel_bin1 = np.floor(table['cut1']/binsize).astype(int)
            rel_bin2 = np.floor(table['cut2']/binsize).astype(int)
            table['bin1'] = chrom_offset[table['chrom_id1'].values] + rel_bin1
            table['bin2'] = chrom_offset[table['chrom_id2'].values] + rel_bin2

        # reduce
        gby = table.groupby(['bin1', 'bin2'])
        agg = (gby['chrom_id1'].count()
                               .reset_index()
                               .rename(columns={'chrom_id1': 'count'}))
        n = len(agg)

        # insert new pixels
        for dset in [Bin1, Bin2, Count]:
            dset.resize((nnz + n,))
        Bin1[nnz:nnz+n] = agg['bin1'].values
        Bin2[nnz:nnz+n] = agg['bin2'].values
        Count[nnz:nnz+n] = agg['count'].values

        # add new pixeltable rows to the bin1 offset index
        bin1lo, bin1hi = bin1hi, i + 1
        bin1_range = np.arange(bin1lo, bin1hi)
        bin1_offset[bin1lo:bin1hi] = nnz + np.searchsorted(
            agg['bin1'].values, bin1_range, side='left')
 
        nnz += n

    bin1_offset[bin1hi:] = nnz

    return chrom_offset, bin1_offset, nnz


def from_readhdf5(h5, chromtable, bintable, h5read,
                  binsize=None, info=None, h5opts=None, chunksize=40000000,
                  check_sorted=True):
    h5opts = {'compression': 'lzf'} if h5opts is None else h5opts
    info = {} if info is None else info

    n_chroms = len(chromtable)
    n_bins = len(bintable)

    print('scaffolds')
    grp = h5.create_group('scaffolds')
    write_chromtable(grp, chromtable, h5opts)

    print('bins')
    grp = h5.create_group('bins')
    write_bintable(grp, chromtable, bintable, h5opts)

    print('matrix')
    grp = h5.create_group('matrix')
    chrom_offset, bin1_offset, nnz = _aggregate(
        grp, chromtable, bintable, h5read, binsize, h5opts, chunksize,
        check_sorted)

    print('indexes')
    grp = h5.create_group('indexes')
    write_indexes(grp, chrom_offset, bin1_offset, h5opts)

    print('info')
    info['bin-type'] = 'fixed' if binsize is not None else 'variable'
    info['bin-size'] = binsize if binsize is not None else 'null'
    info['nchroms'] = n_chroms
    info['nbins'] = n_bins
    info['nnz'] = nnz
    write_info(h5, info)


def from_dense(h5, chromtable, bintable, heatmap,
               binsize=None, h5opts=None, info=None):
    h5opts = {'compression': 'lzf'} if h5opts is None else h5opts
    info = {} if info is None else info

    n_chroms = len(chromtable)
    n_bins = len(bintable)
    if len(heatmap) != n_bins:
        raise ValueError(
            "length mismatch:" +
            " heatmap length is {0}, bin table length is {1}".format(
                len(heatmap), n_bins))

    print('scaffolds')
    grp = h5.create_group('scaffolds')
    write_chromtable(grp, chromtable, h5opts)

    print('bins')
    bintype = 'variable' if binsize is None else 'fixed'
    grp = h5.create_group('bins')
    write_bintable(grp, chromtable, bintable, h5opts)

    print('matrix')
    # TRIU sparsify the matrix
    i, j = np.nonzero(heatmap)
    mask = i <= j
    triu_i, triu_j = i[mask], j[mask]
    values = heatmap[triu_i, triu_j]
    nnz = len(values)
    grp = h5.create_group('matrix')
    grp.create_dataset('bin1_id',
                       shape=(len(values),),
                       dtype=np.int32,
                       data=triu_i, **h5opts)
    grp.create_dataset('bin2_id',
                       shape=(len(values),),
                       dtype=np.int32,
                       data=triu_j, **h5opts)
    grp.create_dataset('count',
                       shape=(len(values),),
                       dtype=np.int32,
                       data=values, **h5opts)

    print('indexes')
    grp = h5.create_group('indexes')
    chrom_offset = np.r_[0, np.cumsum(np.ceil(chromtable['length']/binsize))]
    bin1_offset = np.r_[np.searchsorted(triu_i, np.arange(n_bins), side='left'), nnz]
    write_indexes(grp, chrom_offset, bin1_offset, h5opts)

    print('info')
    info['bin-type'] = bintype
    info['bin-size'] = binsize
    info['nchroms'] = n_chroms
    info['nbins'] = n_bins
    info['nnz'] = nnz
    write_info(h5, info)
