# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals
from contextlib import contextmanager
from multiprocessing import Pool
from datetime import datetime
from copy import copy
import warnings
import six

from scipy.sparse import coo_matrix
import numpy as np
import numexpr
import pandas
import h5py

from hiclib.hicShared import h5dictBinarySearch, binarySearch
from mirnylib.genome import Genome


def make_chrom_table(genome):
    return pandas.DataFrame({
        'chrom_id': np.arange(genome.chrmCount),
        'name': genome.chrmLabels,
        'length': genome.chrmLens},
        columns=['chrom_id', 'name', 'length'])


def make_bin_table(genome, binsize):
    if binsize is not "auto":
        genome.setResolution(binsize)
    chromID = genome.chrmIdxBinCont
    chrmName = [genome.idx2label[i] for i in genome.chrmIdxBinCont]
    start = genome.posBinCont
    end = np.minimum(genome.posBinCont + binsize, genome.chrmLens[chromID])
    GC = np.concatenate(genome.GCBin)
    GC[GC == -1] = np.nan
    return pandas.DataFrame({
        "chrom_id": chromID, 
        "start": start, 
        "end": end, 
        "GC": GC},
        columns=['chrom_id', 'start', 'end', 'GC'])


class _FragmentGrouper(object):
    def __init__(self, chromtable, bintable, binsize):
        # Create a reverse lookup table to assign fragments to bins:
        # hash -> binID where hash = chrmID * mult + pos / resolution
        mult = int(np.ceil(chromtable['length'].max() / binsize)) # + 1
        chrom_ids = bintable['chrom_id'].astype(np.int64)
        offset = bintable['start'] // binsize
        hashes = chrom_ids * mult + offset
        n_bins = len(bintable)
        self.binsize = binsize
        self.mult = mult
        self.lookup = np.zeros(hashes.max() + 1, dtype=np.int32)
        self.lookup[hashes] = np.arange(n_bins)

    def count(self, chrom1, pos1, chrom2, pos2):
        # assign pairs to bins
        binsize = self.binsize
        mult = self.mult
        bin1_ids = self.lookup[numexpr.evaluate("chrom1 * mult + pos1 / binsize").astype(int)]
        bin2_ids = self.lookup[numexpr.evaluate("chrom2 * mult + pos2 / binsize").astype(int)]
        # make a 2-element hash to aggregate and count unique bin pairs
        S = len(self.lookup) + 1000
        pair_hashes = numexpr.evaluate("bin1_ids * S + bin2_ids")
        pair_hashes_u, counts = np.unique(pair_hashes, return_counts=True)
        bin1_ids_u = numexpr.evaluate("pair_hashes_u / S").astype(int)
        bin2_ids_u = numexpr.evaluate("pair_hashes_u % S").astype(int)
        return bin1_ids_u, bin2_ids_u, counts


# Still broken...
def store_from_fragment_map(h5frag, h5out, genome, binsize, chunksize=40000000):
    chrom_table = make_chrom_table(genome)
    bin_table = make_bin_table(genome, binsize)
    n_bins  = len(bin_table)

    grp = h5out.create_group('chromosomes')
    grp["name"] = np.array(chrom_table['name'], dtype='S')
    grp["length"] = chrom_table['length']

    grp = h5out.create_group('bins')
    grp["chrom_id"] = bin_table['chrom_id']
    grp["start"] = bin_table['start']
    grp["end"] = bin_table['end']

    grp = h5out.create_group('heatmap')
    FragChrms1, FragCuts1 = h5frag["chrms1"], h5frag["cuts1"]
    FragChrms2, FragCuts2 = h5frag["chrms2"], h5frag["cuts2"]
    n_records = len(FragChrms1)

    init_size = 5 * n_bins
    max_size  = min(n_records, n_bins * (n_bins - 1) // 2) + 10000
    HMBin1  = grp.create_dataset(
        'bin1', dtype=np.int32, shape=(init_size,), maxshape=(max_size,))
    HMBin2  = grp.create_dataset(
        'bin2', dtype=np.int32, shape=(init_size,), maxshape=(max_size,))
    HMCount = grp.create_dataset(
        'count', dtype=np.int32, shape=(init_size,), maxshape=(max_size,))
    grouper = _FragmentGrouper(chrom_table, bin_table, binsize)
    starts = np.zeros(n_bins, dtype=int)
    bin_lo, bin_hi = 0, 0
    frag_lo, frag_hi = 0, 0
    ptr = 0
    while True:
        frag_hi = min(frag_hi + chunksize, n_records)
        chrom = FragChrms1[frag_hi-1]
        pos = int(np.ceil(FragCuts1[frag_hi-1]/binsize)) * binsize
        frag_hi = h5dictBinarySearch(FragChrms1, FragCuts1, (chrom, pos), 'right')
        
        # bin the fragment counts
        c1, p1 = FragChrms1[frag_lo:frag_hi], FragCuts1[frag_lo:frag_hi]
        c2, p2 = FragChrms2[frag_lo:frag_hi], FragCuts2[frag_lo:frag_hi]
        bin1_ids, bin2_ids, counts = grouper.count(c1, p1, c2, p2)
        n_unique = len(counts)
        
        # insert new heatmap records
        HMBin1[ptr:ptr+n_unique] = bin1_ids
        HMBin2[ptr:ptr+n_unique] = bin2_ids
        HMCount[ptr:ptr+n_unique] = counts
        for dset in [HMBin1, HMBin2, HMCount]:
            dset.resize((ptr + n_unique,))

        # update the bin_to_heatmap index
        bin_lo, bin_hi = bin_hi, bin1_ids.max() + 1
        starts[bin_lo:bin_hi] = ptr + np.searchsorted(
            bin1_ids, np.arange(bin_lo, bin_hi), side='left')
        ptr += n_unique
        
        if frag_hi == n_records:
            break
        
    grp = h5out.create_group('indexes')
    chrom_binedges = np.r_[0, np.cumsum(g.chrmLensBin)]
    idx1 = grp.create_group('chrom2bin')
    idx1["bin_lo"] = chrom_binedges[:-1]
    idx1["bin_hi"] = chrom_binedges[1:]
    
    nnz = len(HMCount)
    ends = np.r_[starts[1:], nnz]
    idx2 = grp.create_group('bin2heatmap')
    idx2["heatmap_lo"] = starts
    idx2["heatmap_hi"] = ends

