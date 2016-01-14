# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals
from contextlib import contextmanager
from collections import OrderedDict
from datetime import datetime
import six

import numpy as np
import pandas
import h5py

from . import __version__
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
        if mode =='r' and fp.mode == 'r+':
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
    n_bins  = len(bintable)
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
    n_chroms = len(chrom_offset)
    n_bins  = len(bin1_offset)
    grp.create_dataset("chrom_offset",
                        shape=(len(chrom_offset),), dtype=np.int32, 
                        data=chrom_offset, **h5opts)
    grp.create_dataset("bin_offset",
                        shape=(len(bin1_offset),), dtype=np.int32,
                        data=bin1_offset, **h5opts)


def _from_fraghdf5(grp, chromtable, bintable, h5frag, binsize, h5opts, chunksize):
    n_records = len(h5frag["chrms1"])
    n_chroms = len(chromtable)
    n_bins  = len(bintable)

    def load_chunk(h5frag, lo, hi):
        return pandas.DataFrame(OrderedDict([
            ('chrom1', h5frag['chrms1'][lo:hi]),
            ('cut1', h5frag['cuts1'][lo:hi]),
            ('strand1', h5frag['strands1'][lo:hi]),
            ('chrom2', h5frag['chrms2'][lo:hi]),
            ('cut2', h5frag['cuts2'][lo:hi]),
            ('strand2', h5frag['strands2'][lo:hi]),
        ]))
    
    init_size = 5 * n_bins
    max_size  = min(n_records, n_bins * (n_bins - 1) // 2) + 10000
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
    
    chrom_lengths_bin = np.ceil(chromtable['length'].values/binsize).astype(int)
    chrom_offset = np.r_[0, np.cumsum(chrom_lengths_bin)]
    binedges = np.arange(0, int(chromtable['length'].max()) + binsize, binsize)
    bin_offset = np.zeros(n_bins+1, dtype=np.int64)
    bin_hi = 0
    hi = 0
    i = 0
    while True:
        lo, hi = hi, min(hi + chunksize, n_records)
        chrom = h5frag["chrms1"][hi - 1]
        pos = h5frag["cuts1"][hi - 1]
        pos_floor = int(np.ceil(pos/binsize)) * binsize
        hi = lexbisect((h5frag["chrms1"], h5frag["cuts1"]), (chrom, pos_floor), 'right')
        
        table = load_chunk(h5frag, lo, hi)
        coord_bin1 = np.floor(table['cut1']/binsize).astype(int)
        coord_bin2 = np.floor(table['cut2']/binsize).astype(int)
        table['bin1'] = chrom_offset[table['chrom1']] + coord_bin1
        table['bin2'] = chrom_offset[table['chrom2']] + coord_bin2
        print(lo, hi)
        
        gby = table.groupby(['bin1', 'bin2'])
        agg = gby['chrom1'].count().reset_index().rename(columns={'chrom1':'count'})
        n_unique = len(agg)

        # insert new matrix elements
        for dset in [Bin1, Bin2, Count]:
            dset.resize((i + n_unique,))
        Bin1[i:i+n_unique] = agg['bin1']
        Bin2[i:i+n_unique] = agg['bin2']
        Count[i:i+n_unique] = agg['count']
        
        # add new matrix rows to the offset index
        bin_lo, bin_hi = bin_hi, agg['bin1'].max()
        bin_offset[bin_lo:bin_hi] = i + np.searchsorted(agg['bin1'], np.arange(bin_lo, bin_hi), side='left')
    
        i += n_unique
        if hi == n_records:
            break

    nnz = len(Count)
    bin_offset[n_bins] = nnz

    return chrom_offset, bin_offset, nnz


def from_fraghdf5(h5, chromtable, bintable, h5frag, binsize=None, info=None, h5opts=None, chunksize=40000000):
    h5opts = {'compression': 'lzf'} if h5opts is None else h5opts
    info = {} if info is None else info

    n_records = len(h5frag["chrms1"])
    n_chroms = len(chromtable)
    n_bins  = len(bintable)
    
    print('sequence assemblies')
    grp = h5.create_group('scaffolds')
    write_chromtable(grp, chromtable, h5opts)
    
    print('bins')
    bintype = 'fixed'
    if binsize is None:
        raise ValueError("variable bin size not yet supported")
    grp = h5.create_group('bins')
    write_bintable(grp, chromtable, bintable, h5opts)
    
    print('matrix')
    grp = h5.create_group('matrix')
    chrom_offset, bin1_offset, nnz = _from_fraghdf5(grp, chromtable, bintable, h5frag, binsize, h5opts, chunksize)
    
    print('indexes')
    grp = h5.create_group('indexes') 
    write_indexes(grp, chrom_offset, bin1_offset, h5opts)
    
    print('info')
    h5.attrs['id'] = info.get('id', "No ID")
    h5.attrs['generated-by'] = info.get('generated-by', "cooler")
    h5.attrs['creation-date'] = datetime.now().isoformat()
    h5.attrs['format-version'] = __version__
    h5.attrs['format-url'] = "https://github.com/mirnylab/cooler"
    h5.attrs['genome-assembly'] = info.get('genome-assembly', 'unknown')
    h5.attrs['bin-type'] = bintype
    h5.attrs['bin-size'] = binsize
    h5.attrs['nchroms'] = n_chroms
    h5.attrs['nbins'] = n_bins
    h5.attrs['nnz'] = nnz


def from_dense(h5, chromtable, bintable, heatmap, binsize=None, h5opts=None, info=None):
    h5opts = {'compression': 'lzf'} if h5opts is None else h5opts
    info = {} if info is None else info

    n_chroms = len(chromtable)
    n_bins = len(bintable)
    if len(heatmap) != n_bins:
        raise ValueError(
            "length mismatch:" + 
            " heatmap length is {0}, bin table length is {1}".format(
                len(heatmap), n_bins))

    print('sequence assemblies')
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
    h5.attrs['id'] = info.get('id', "No ID")
    h5.attrs['generated-by'] = info.get('generated-by', "cooler")
    h5.attrs['creation-date'] = datetime.now().isoformat()
    h5.attrs['format-version'] = __version__
    h5.attrs['format-url'] = "https://github.com/mirnylab/cooler"
    h5.attrs['genome-assembly'] = info.get('genome-assembly', 'unknown')
    h5.attrs['bin-type'] = bintype
    h5.attrs['bin-size'] = binsize
    h5.attrs['nchroms'] = len(chromtable)
    h5.attrs['nbins'] = heatmap.shape[0]
    h5.attrs['nnz'] = nnz
