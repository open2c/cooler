# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals
import six

import numpy as np
import numexpr
import pandas
import h5py

from hiclib.hicShared import h5dictBinarySearch, binarySearch
from mirnylib.genome import Genome


def make_chrom_table(genome):
    names = ['chr{}'.format(i) for i in genome.chrmLabels]
    return pandas.DataFrame({
        'name': names,
        'id': np.arange(genome.chrmCount),
        'length': genome.chrmLens},
        columns=['name', 'id', 'length'],
        index=names)


def make_bin_table(genome, binsize):
    if binsize is not "auto":
        genome.setResolution(binsize)
    chrmIDs = genome.chrmIdxBinCont
    chrmNames = ['chr{}'.format(genome.idx2label[i]) for i in chrmIDs]
    starts = genome.posBinCont
    ends = np.minimum(genome.posBinCont + binsize, genome.chrmLens[chrmIDs])
    #GC = np.concatenate(genome.GCBin)
    #GC[GC == -1] = np.nan
    return pandas.DataFrame({
        "chrom": chrmNames, 
        "start": starts, 
        "end": ends},
        columns=['chrom', 'start', 'end'])


def _load(h5frag, lo, hi):
    return pandas.DataFrame(OrderedDict([
        ('chrom1', h5frag['chrms1'][lo:hi]),
        ('cut1', h5frag['cuts1'][lo:hi]),
        ('strand1', h5frag['strands1'][lo:hi]),
        ('chrom2', h5frag['chrms2'][lo:hi]),
        ('cut2', h5frag['cuts2'][lo:hi]),
        ('strand2', h5frag['strands2'][lo:hi]),
    ]))


def write_chromtable(grp, chrom_table, h5opts):
    n_chroms = len(chrom_table)
    grp.create_dataset('name', 
                        shape=(n_chroms,), 
                        dtype='S32',
                        data=np.array(chrom_table['name'], dtype='S32'),
                        **h5opts)
    grp.create_dataset('length', 
                        shape=(n_chroms,), 
                        dtype=np.int64,
                        data=chrom_table['length'],
                        **h5opts)

    
def write_bintable(grp, chrom_table, bin_table, h5opts):
    n_chroms = len(chrom_table)
    n_bins  = len(bin_table)
    if 'id' not in chrom_table.columns:
        chrom_table['id'] = np.arange(n_chroms)
    chrom_ids = chrom_table['id'].loc[bin_table['chrom']]
    grp.create_dataset('chrom_id',
                        shape=(n_bins,),
                        dtype=np.int32,
                        data=chrom_ids,
                        **h5opts)
    grp.create_dataset('start',
                        shape=(n_bins,),
                        dtype=np.int64,
                        data=bin_table['start'], 
                        **h5opts)
    grp.create_dataset('end',
                        shape=(n_bins,),
                        dtype=np.int64,
                        data=bin_table['end'], 
                        **h5opts)
    
def write_matrix(grp, chrom_table, bin_table, h5frag, h5opts, chunksize):
    n_records = len(h5frag["chrms1"])
    n_chroms = len(chrom_table)
    n_bins  = len(bin_table)
    
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
    
    chrom_lengths_bin = np.ceil(chrom_table['length'].values/binsize).astype(int)
    chrom_offset = np.r_[0, np.cumsum(chrom_lengths_bin)]
    binedges = np.arange(0, int(chrom_table['length'].max()) + binsize, binsize)
    binlabels = np.arange(len(binedges)-1)
    
    mat_offset = np.zeros(n_bins+1, dtype=np.int64)
    bin_hi = 0
    hi = 0
    i = 0
    while True:
        lo, hi = hi, min(hi + chunksize, n_records)
        chrom = h5frag["chrms1"][hi - 1]
        pos = h5frag["cuts1"][hi - 1]
        pos_floor = int(np.ceil(pos/binsize)) * binsize
        hi = h5dictBinarySearch(h5frag["chrms1"], h5frag["cuts1"], (chrom, pos_floor), 'right')
        
        table = _load(h5frag, lo, hi)
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
        
        # add to the bin_to_matrix index
        bin_lo, bin_hi = bin_hi, agg['bin1'].max()
        mat_offset[bin_lo:bin_hi] = i + np.searchsorted(agg['bin1'], np.arange(bin_lo, bin_hi), side='left')
    
        i += n_unique
        if hi == n_records:
            break
    nnz = len(Count)
    mat_offset[n_bins] = nnz
    return chrom_offset, mat_offset, nnz
    

def from_readtable(h5, chrom_table, bin_table, binsize, h5frag, h5opts, metadata, chunksize=40000000):
    n_records = len(h5frag["chrms1"])
    n_chroms = len(chrom_table)
    n_bins  = len(bin_table)
    
    print('sequence assemblies')
    grp = h5.create_group('scaffolds')
    write_chromtable(grp, chrom_table, h5opts)
    
    print('bins')
    grp = h5.create_group('bins')
    write_bintable(grp, chrom_table, bin_table, h5opts)
    
    print('matrix')
    grp = h5.create_group('matrix')
    chrom_offset, mat_offset, nnz = write_matrix(grp, chrom_table, bin_table, h5frag, h5opts, chunksize)
    
    print('indexes')
    grp = h5.create_group('indexes') 
    idx1 = grp.create_group('chrom_to_bin')
    idx1["bin_lo"] = chrom_offset[:-1]
    idx1["bin_hi"] = chrom_offset[1:]
    idx2 = grp.create_group('bin_to_matrix')
    idx2["mat_lo"] = mat_offset[:-1]
    idx2["mat_hi"] = mat_offset[1:]
    
    # Attributes
    print('metadata')
    h5.attrs['id'] = metadata.get('id', "No ID")
    h5.attrs['bin-type'] = 'fixed'
    h5.attrs['bin-size'] = binsize
    h5.attrs['genome-assembly'] = metadata.get('genome-assembly', 'unknown')
    h5.attrs['format-url'] = "https://bitbucket.org/nvictus/cooler"
    h5.attrs['format-version'] = (0, 1)
    h5.attrs['generated-by'] = metadata.get('generated-by', "cooler")
    h5.attrs['creation-date'] = datetime.now().isoformat()
    h5.attrs['shape'] = (n_bins, n_bins)
    h5.attrs['nnz'] = nnz

