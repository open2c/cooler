from __future__ import division, print_function, unicode_literals
from contextlib import contextmanager
from multiprocessing import Pool
from copy import copy
import warnings

from scipy.sparse import coo_matrix
import numpy as np
import numexpr
import pandas
import h5py
import six

from hiclib.hicShared import h5dictBinarySearch, binarySearch
from mirnylib.genome import Genome


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


def to_tsv(fp, table, index=False, **kwargs):
    kwargs['index'] = index
    table.to_csv(fp, sep=b'\t', **kwargs)


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


def _get_chunks(N, chunksize):
    borders = range(0, N, chunksize) + [N]
    if len(borders) > 2:
        if borders[-1] - borders[-2] < 0.3 * chunksize:
            borders.pop(-2)
    return zip(borders[:-1], borders[1:])


def _get_random_chunks(N, chunksize):
    a = np.random.randint(chunksize / 2, chunksize, 2 * N / chunksize)
    b = np.cumsum(a)
    c = np.searchsorted(b, N - chunksize / 3)
    borders = np.r_[0, b[:c], N]
    return zip(borders[:-1], borders[1:])


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


# This works...
def store_from_dense(heatmap, h5out, genome, binsize):
    chrom_table = make_chrom_table(genome)
    bin_table = make_bin_table(genome, binsize)
    n_bins = len(bin_table)
    if len(heatmap) != n_bins:
        raise ValueError(
            "length mismatch: heatmap length is {0}, table length is {1}".format(
                len(heatmap), len(bin_table)))

    i, j = np.nonzero(heatmap)
    mask = i >= j
    tril_i = i[mask]
    tril_j = j[mask]
    values = heatmap[tril_i, tril_j]
    nnz = len(values)
        
    grp = h5out.create_group('chromosomes')
    grp["name"] = np.array(chrom_table['name'], dtype='S')
    grp["length"] = chrom_table['length']

    grp = h5out.create_group('bins')
    grp["chrom_id"] = bin_table['chrom_id']
    grp["start"] = bin_table['start']
    grp["end"] = bin_table['end']

    grp = h5out.create_group('heatmap')
    grp["bin1"]  = np.array(tril_i, dtype=np.int32)
    grp["bin2"]  = np.array(tril_j, dtype=np.int32)
    grp["count"] = np.array(values, dtype=np.int32)

    grp = h5out.create_group('indexes')
    
    chrom_binedges = np.r_[0, np.cumsum(g.chrmLensBin)]
    idx1 = grp.create_group('chrom2bin')
    idx1["bin_lo"] = chrom_binedges[:-1]
    idx1["bin_hi"] = chrom_binedges[1:]
    
    starts = np.searchsorted(triu_i, np.arange(n_bins), side='left')
    ends = np.r_[starts[1:], nnz]
    idx2 = grp.create_group('bin2heatmap')
    idx2["heatmap_lo"] = starts
    idx2["heatmap_hi"] = ends


def parse_genomic_region(chromtable, region):
    chrom_lengths = chromtable['length']
    chrom_ids = chromtable['id']
    try:
        if isinstance(region, six.string_types):
            chrom, start = region, 0
            end = chrom_lengths[chrom]
        else:
            chrom, start, end = region
        chrom_id = chrom_ids[chrom]
        chrom_len = chrom_lengths[chrom]
    except KeyError:
        raise ValueError("Invalid chromosome name or integer ID")
    if start < 0 or end >= chrom_len:
        raise ValueError(
            "Genomic coordinates out of bounds: [0, {})".format(chrom_len))
    return chrom_id, start, end


def get_bin_offset(h5, chrom_id):
    return h5['indexes']['chrom2bin']['bin_lo'][chrom_id]


def get_bin_extent(h5, region, binsize=None):
    chrom_id, start, end = region
    if binsize is not None:
        bin_chrom_lo = h5['indexes']['chrom2bin']['bin_lo'][chrom_id]
        bin_lo = bin_chrom_lo + int(np.floor(start/binsize))
        bin_hi = bin_chrom_lo + int(np.ceil(end/binsize))
    else:
        bin_chrom_lo = h5['indexes']['chrom2bin']['bin_lo'][chrom_id]
        bin_chrom_hi = h5['indexes']['chrom2bin']['bin_hi'][chrom_id]
        bin_chrom = h5['bins']['start'][bin_chrom_lo:bin_chrom_hi]
        bin_lo = bin_chrom_lo + np.searchsorted(bin_chrom, start, 'left')
        bin_hi = bin_chrom_lo + np.searchsorted(bin_chrom, end, 'right')
    return bin_lo, bin_hi


def slice_heatmap(h5, bin1_lo, bin1_hi, bin2_lo, bin2_hi):
    bin1_ids = np.arange(bin1_lo, bin1_hi)    
    hm1_lo = h5['indexes']['bin2heatmap']['heatmap_lo'][bin1_lo:bin1_hi]
    hm1_hi = h5['indexes']['bin2heatmap']['heatmap_hi'][bin1_lo:bin1_hi]
    i, j, v = [], [], []
    for bin1_id, lo1, hi1 in zip(bin1_ids, hm1_lo, hm1_hi):
        hm_bin2 = h5['heatmap']['bin2'][lo1:hi1]
        lo = lo1 + np.searchsorted(hm_bin2, bin2_lo)
        hi = lo1 + np.searchsorted(hm_bin2, bin2_hi)
        i.append(np.zeros(hi-lo) + bin1_id)
        j.append(h5['heatmap']['bin2'][lo:hi])
        v.append(h5['heatmap']['count'][lo:hi])
    i = np.concatenate(i, axis=0).astype(int)
    j = np.concatenate(j, axis=0).astype(int)
    v = np.concatenate(v, axis=0)
    return i, j, v


class Cooler(object):
    def __init__(self, h5path, binsize=None):
        self.filepath = h5path
        with open_hdf5(self.filepath, 'r') as h5:
            self._binsize = h5.attrs.get('binsize', binsize)
            self._chrom_table = pandas.DataFrame({
                'id': np.arange(h5['chromosomes']['length'].shape[0]),
                'length': h5['chromosomes']['length'][:]},
                    columns=['id', 'length'],
                    index=h5['chromosomes']['name'][:])
            self._bin_table = pandas.DataFrame({
                'chrom': self._chrom_table.index[h5['bins']['chrom_id'][:]].values,
                'start': h5['bins']['start'],
                'end': h5['bins']['end'],
                }, columns=['chrom', 'start', 'end'])

    @property
    def binsize(self):
        return self._binsize

    @property
    def chromosomes(self):
        return self._chrom_table.copy()

    @property
    def bins(self):
        return self._bin_table.copy()

    def load_bins(self, region):
        pass

    def load_track(self, region):
        pass

    def load(self, region, region2=None, sparse=False, as_triples=False):
        region1 = parse_genomic_region(self._chrom_table, region)
        if region2 is not None:
            region2 = parse_genomic_region(self._chrom_table, region2)
        else:
            region2 = region1
        if region1 != region2 and overlaps(region1, region2):
            i, j, v = [], [], []
            m, n = 0, 0
            for r1, r2 in split(region1, region2):
                (ii, jj, vv), (mm, nn) = self._load(r1, r2)
                i.append(ii)
                j.append(jj)
                v.append(vv)
                m += mm
                n += nn
            i = np.concatenate(i, axis=0)
            j = np.concatenate(j, axis=0)
            v = np.concatenate(v, axis=0)
        else:
            (i, j, v), (m, n) = self._load(region, region2)
        if as_triples:
            return i, j, v
        A = coo_matrix((v, (i, j)), (m, n))
        if not sparse:
            A = A.toarray()
        return A

    def _load(self, region, region2):
        binsize = self._binsize
        symmetric = False
        transpose = False
        if region1 == region2:
            symmetric = True
        elif region1 < region2:
            region1, region2 = region2, region1
            transpose = True
                    
        with open_hdf5(self.filepath) as h5:
            bin1_lo, bin1_hi = get_bin_extent(h5, region1, self._binsize)
            if symmetric:
                bin2_lo, bin2_hi = bin1_lo, bin1_hi
            else:
                bin2_lo, bin2_hi = get_bin_extent(h5, region2, self._binsize)
            i, j, v = slice_heatmap(h5, bin1_lo, bin1_hi, bin2_lo, bin2_hi)
        m = bin1_hi - bin1_lo
        n = bin2_hi - bin2_lo
        i -= bin1_lo
        j -= bin2_lo

        if symmetric:
            i, j, v = np.r_[i, j], np.r_[j, i], np.r_[v, v]
        if transpose:
            return (j, i, v), (n, m)
        return (i, j, v), (m, n)

