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


class _PixelGrouper(object):
    def __init__(self, chromtable, bintable, binsize):
        # Create a reverse lookup table to assign reads to bins:
        # hash -> binID where hash = chrmID * mult + pos / resolution
        mult = int(np.ceil(chromtable['length'].max() / binsize)) # + 1
        chrom_ids = chromtable['id'].loc[bintable['chrom']].values
        chrom_ids = chrom_ids.astype(np.int64)
        offset = bintable['start'] // binsize
        hashes = chrom_ids * mult + offset
        n_bins = len(bintable)
        self.binsize = binsize
        self.mult = mult
        self.lookup = np.zeros(hashes.max() + 1, dtype=np.int32)
        self.lookup[hashes] = np.arange(n_bins)

    def count(self, chrom1, pos1, chrom2, pos2):
        """
        Count number of read pairs provided that fall inside the pixels of the 
        2-D grid.

        """
        # assign read pairs to pixels
        binsize = self.binsize
        mult = self.mult
        bin1_ids = self.lookup[
            numexpr.evaluate("chrom1 * mult + pos1 / binsize").astype(int)]
        bin2_ids = self.lookup[
            numexpr.evaluate("chrom2 * mult + pos2 / binsize").astype(int)]
        # make a two-element hash to aggregate and count pixel occurrences
        S = len(self.lookup) + 1000
        pair_hashes = numexpr.evaluate("bin1_ids * S + bin2_ids")
        pair_hashes_u, counts = np.unique(pair_hashes, return_counts=True)
        bin1_ids_u = numexpr.evaluate("pair_hashes_u / S").astype(int)
        bin2_ids_u = numexpr.evaluate("pair_hashes_u % S").astype(int)
        return bin1_ids_u, bin2_ids_u, counts


def from_read_table(h5, chrom_table, bin_table, binsize, h5frag, h5opts, 
                    chunksize=40000000):
    n_bins  = len(bin_table)
    n_chroms = len(chrom_table)
    FragChrms1, FragCuts1 = h5frag["chrms1"], h5frag["cuts1"]
    FragChrms2, FragCuts2 = h5frag["chrms2"], h5frag["cuts2"]
    n_records = len(FragChrms1)

    print('sequence assemblies')
    grp = h5.create_group('scaffolds')
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

    print('bins')
    grp = h5.create_group('bins')
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

    print('matrix')
    grp = h5.create_group('matrix')
    grouper = _PixelGrouper(chrom_table, bin_table, binsize)
    init_size = 5 * n_bins
    max_size  = min(n_records, n_bins * (n_bins - 1) // 2) + 10000
    HMBin1  = grp.create_dataset('bin1_id', 
                                 dtype=np.int64, 
                                 shape=(init_size,), 
                                 maxshape=(max_size,))
    HMBin2  = grp.create_dataset('bin2_id',
                                 dtype=np.int64, 
                                 shape=(init_size,),
                                 maxshape=(max_size,))
    HMCount = grp.create_dataset('count', 
                                 dtype=np.int64, 
                                 shape=(init_size,), 
                                 maxshape=(max_size,))
    mat_lo = np.zeros(n_bins, dtype=np.int64)
    bin_lo, bin_hi = 0, 0
    frag_lo, frag_hi = 0, 0
    i = 0
    while True:
        frag_hi = min(frag_hi + chunksize, n_records)
        chrom = FragChrms1[frag_hi - 1]
        pos = int(np.ceil(FragCuts1[frag_hi - 1]/binsize)) * binsize
        frag_hi = h5dictBinarySearch(
            FragChrms1, FragCuts1, (chrom, pos), 'right')
        # bin the fragment counts
        c1, p1 = FragChrms1[frag_lo:frag_hi], FragCuts1[frag_lo:frag_hi]
        c2, p2 = FragChrms2[frag_lo:frag_hi], FragCuts2[frag_lo:frag_hi]
        bin1_ids, bin2_ids, counts = grouper.count(c1, p1, c2, p2)
        n_unique = len(counts)
        for dset in [HMBin1, HMBin2, HMCount]:
            dset.resize((i + n_unique,))
        # insert new matrix elements
        HMBin1[i:i+n_unique] = bin1_ids
        HMBin2[i:i+n_unique] = bin2_ids
        HMCount[i:i+n_unique] = counts
        # add to the bin_to_matrix index
        bin_lo, bin_hi = bin_hi, bin1_ids.max()
        mat_lo[bin_lo:bin_hi] = i + np.searchsorted(
            bin1_ids, np.arange(bin_lo, bin_hi), side='left')
        i += n_unique
        if frag_hi == n_records:
            break
    
    print('indexes')
    grp = h5.create_group('indexes')
    nnz = len(HMCount)
    mat_hi = np.r_[mat_lo[1:], nnz].astype(np.int32)
    chrom_lengths_bin = np.ceil(chrom_table['length'].values/binsize)
    chrom_binedges = np.r_[0, np.cumsum(chrom_lengths_bin)].astype(np.int64)
    idx1 = grp.create_group('chrom_to_bin')
    idx1["bin_lo"] = chrom_binedges[:-1]
    idx1["bin_hi"] = chrom_binedges[1:]
    idx2 = grp.create_group('bin_to_matrix')
    idx2["mat_lo"] = mat_lo
    idx2["mat_hi"] = mat_hi


