# -*- coding: utf-8 -*-
"""
Contact Readers
~~~~~~~~~~~~~~~

Reader classes convert input data of various flavors into a chunked stream of
binned contacts.

"""
from __future__ import division, print_function
from collections import OrderedDict, Counter
from contextlib import contextmanager
from bisect import bisect_left
from multiprocess import Pool
import subprocess
import itertools
import warnings
import json
import sys
import six

from pandas.algos import is_lexsorted
import numpy as np
import pandas
import h5py

from . import get_logger
from .util import rlencode, get_binsize


logger = get_logger()


class ContactReader(object):
    """
    Interface of a contact reader.

    """
    def size(self):
        """ Total number of contacts """
        raise NotImplementedError

    def __iter__(self):
        """ Iterator over chunks of binned contacts

        Chunks are expected to have the following format:

        * dict of 1D arrays
        * keys `bin1_id`, `bin2_id`, `count`
        * arrays lexically sorted by `bin_id` then `bin2_id`

        """
        raise NotImplementedError


def check_bins(bins, chromsizes):
    if bins['chrom'].dtype.name != 'category':
        bins['chrom'] = pandas.Categorical(
            bins.chrom, 
            categories=list(chromsizes.index), 
            ordered=True)
    else:
        assert (bins['chrom'].cat.categories == chromsizes.index).all()
    return bins


def parse_genomic_segmentation(chromsizes, bins):
    idmap = pandas.Series(index=chromsizes.keys(), data=range(len(chromsizes)))
    bins_chrom_id = bins['chrom'].cat.codes
    chrom_posoffset = np.r_[0, np.cumsum(chromsizes.values)]
    start_abspos = chrom_posoffset[bins_chrom_id] + bins['start'].values
    chrom_binoffset = np.r_[0, np.cumsum(
        bins.groupby('chrom', sort=False).size().values)]
    return idmap, chrom_binoffset, chrom_posoffset, start_abspos


class HDF5Aggregator(ContactReader):
    """
    Aggregate contacts from a hiclib-style HDF5 contacts file.

    """
    def __init__(self, h5pairs, chromsizes, bins, chunksize, **kwargs):
        self.h5 = h5pairs
        self.C1 = kwargs.pop('C1', 'chrms1')
        self.P1 = kwargs.pop('P1', 'cuts1')
        self.C2 = kwargs.pop('C2', 'chrms2')
        self.P2 = kwargs.pop('P2', 'cuts2')
        self.bins = bins
        self.binsize = get_binsize(bins)
        self.n_contacts = len(self.h5[self.C1])
        self.n_bins = len(bins)
        self.chunksize = chunksize
        # convert genomic coords of bin starts to absolute
        self.idmap = pandas.Series(index=chromsizes.keys(), 
                                   data=range(len(chromsizes)))
        bin_chrom_ids = self.idmap[bins['chrom']].values
        self.cumul_length = np.r_[0, np.cumsum(chromsizes)]
        self.abs_start_coords = self.cumul_length[bin_chrom_ids] + bins['start']
        # chrom offset index: chrom_id -> offset in bins
        chrom_nbins =  bins.groupby(bin_chrom_ids, sort=False).size()
        self.chrom_offset = np.r_[0, np.cumsum(chrom_nbins)]
        # index extents of chromosomes on first axis of contact list
        self.partition = self._index_chroms()

    def _index_chroms(self):
        starts, lengths, values = rlencode(self.h5[self.C1], self.chunksize)
        if len(set(values)) != len(values):
            raise ValueError(
                "Read pair coordinates are not sorted on the first axis")
        return dict(zip(values, zip(starts, starts + lengths)))

    def _load_chunk(self, lo, hi):
        data = OrderedDict([
            ('chrom_id1', self.h5[self.C1][lo:hi]),
            ('cut1', self.h5[self.P1][lo:hi]),
            ('chrom_id2', self.h5[self.C2][lo:hi]),
            ('cut2', self.h5[self.P2][lo:hi]),
        ])
        return pandas.DataFrame(data)

    def aggregate(self, chrom):
        h5pairs = self.h5
        C1, P1, C2, P2 = self.C1, self.P1, self.C2, self.P2
        bins = self.bins
        binsize = self.binsize
        chunksize = self.chunksize
        chrom_offset = self.chrom_offset
        cumul_length = self.cumul_length
        abs_start_coords = self.abs_start_coords

        cid = self.idmap[chrom]
        chrom_lo, chrom_hi = self.partition.get(cid, (-1, -1))
        lo = chrom_lo
        hi = lo
        while hi < chrom_hi:
            # fetch next chunk, making sure our selection doesn't split a bin1
            lo, hi = hi, min(hi + chunksize, chrom_hi)
            abs_pos = cumul_length[cid] + h5pairs[P1][hi-1]
            i = int(np.searchsorted(abs_start_coords, abs_pos, 
                    side='right')) - 1
            bin_end = bins['end'][i]
            hi = bisect_left(h5pairs[P1], bin_end, lo, chrom_hi)
            if lo == hi:
                hi = chrom_hi

            logger.info('{} {}'.format(lo, hi))

            # assign bins to reads
            table = self._load_chunk(lo, hi)
            abs_pos1 = (cumul_length[h5pairs[C1][lo:hi]] +
                        h5pairs[P1][lo:hi])
            abs_pos2 = (cumul_length[h5pairs[C2][lo:hi]] +
                        h5pairs[P2][lo:hi])
            if np.any(abs_pos1 > abs_pos2):
                raise ValueError(
                    "Found a read pair that maps to the lower triangle of the "
                    "contact map (side1 > side2). Check that the provided "
                    "chromosome ordering and read pair file are consistent "
                    "such that all pairs map to the upper triangle with "
                    "respect to the given chromosome ordering.")

            if binsize is None:
                table['bin1_id'] = np.searchsorted(
                    abs_start_coords, abs_pos1, side='right') - 1
                table['bin2_id'] = np.searchsorted(
                    abs_start_coords, abs_pos2, side='right') - 1
            else:
                rel_bin1 = np.floor(table['cut1']/binsize).astype(int)
                rel_bin2 = np.floor(table['cut2']/binsize).astype(int)
                table['bin1_id'] = (
                    chrom_offset[table['chrom_id1'].values] + rel_bin1)
                table['bin2_id'] = (
                    chrom_offset[table['chrom_id2'].values] + rel_bin2)

            # reduce
            gby = table.groupby(['bin1_id', 'bin2_id'])
            agg = (gby['chrom_id1'].count()
                                   .reset_index()
                                   .rename(columns={'chrom_id1': 'count'}))
            yield agg

    def size(self):
        return len(self.h5['chrms1'])

    def __iter__(self):
        for chrom in self.idmap.keys():
            for df in self.aggregate(chrom):
                yield {k: v.values for k, v in six.iteritems(df)}


class TabixAggregator(ContactReader):
    """
    Aggregate contacts from a sorted, BGZIP-compressed and tabix-indexed
    tab-delimited text file.

    """
    def __init__(self, filepath, chromsizes, bins, map=map, **kwargs):
        try:
            import pysam
        except ImportError:
            raise ImportError("pysam is required to read tabix files")
        
        self.C2 = kwargs.pop('C2', 3)
        self.P2 = kwargs.pop('P2', 4)
        self._map = map
        self.chromsizes = chromsizes
        self.bins = check_bins(bins, chromsizes)
        self.binsize = get_binsize(bins)
        self.bins_grouped = bins.groupby('chrom', sort=False)
        self.idmap, \
        self.chrom_binoffset, \
        self.chrom_abspos, \
        self.start_abspos = parse_genomic_segmentation(chromsizes, bins)
        
        # read pair records
        self.filepath = filepath
        self.n_records = None
        self.contigs = list(chromsizes.keys())
        with pysam.TabixFile(filepath, 'r', encoding='ascii') as f:
            try:
                file_contigs = [c.decode('ascii') for c in f.contigs]
            except AttributeError:
                file_contigs = f.contigs
        
        for chrom in self.contigs:
            if chrom not in file_contigs:
                warnings.warn(
                    "Did not find contig '{}' in contact list file.".format(chrom))
        self.contigs = [ctg for ctg in self.contigs if ctg in file_contigs]

    def __getstate__(self):
        d = self.__dict__.copy()
        d.pop('_map', None)
        return d

    def _size(self, chrom):
        import pysam
        with pysam.TabixFile(self.filepath, 'r', encoding='ascii') as f:
            return sum(1 for line in f.fetch(chrom))
    
    def size(self):
        if self.n_records is None:
            self.n_records = sum(self._map(self._size, self.contigs))
        return self.n_records
    
    def aggregate(self, chrom):
        import pysam
        filepath = self.filepath
        binsize = self.binsize
        idmap = self.idmap
        chromsizes = self.chromsizes
        chrom_binoffset = self.chrom_binoffset
        chrom_abspos = self.chrom_abspos
        start_abspos = self.start_abspos
        C2, P2 = self.C2, self.P2
        these_bins = self.bins_grouped.get_group(chrom)
        
        rows = []
        with pysam.TabixFile(filepath, 'r', encoding='ascii') as f:
            parser = pysam.asTuple()
            accumulator = Counter()
            
            for bin1_id, bin1 in these_bins.iterrows():
                for line in f.fetch(chrom, bin1.start, bin1.end, parser=parser):
                    chrom2, pos2 = line[C2], int(line[P2])
                    cid2 = idmap[chrom2]
                    if binsize is None:
                        lo, hi = chrom_binoffset[cid2], chrom_binoffset[cid2+1]
                        bin2_id = lo + np.searchsorted(
                            start_abspos[lo:hi], 
                            chrom_abspos[cid2] + pos2,
                            side='right') - 1
                    else:
                        bin2_id = chrom_binoffset[cid2] + (pos2 // binsize)
                    accumulator[bin2_id] += 1

                if not accumulator:
                    continue

                rows.append(
                    pandas.DataFrame({
                        'bin1_id': bin1_id,
                        'bin2_id': list(accumulator.keys()),
                        'count':   list(accumulator.values())},
                        columns=['bin1_id', 'bin2_id', 'count'])
                          .sort_values('bin2_id')
                )
                accumulator.clear()
        
        logger.info(chrom)
        return pandas.concat(rows, axis=0) if len(rows) else None
    
    def __iter__(self):
        for df in self._map(self.aggregate, list(self.contigs)):
            if df is not None:
                yield {k: v.values for k, v in six.iteritems(df)}


class PairixAggregator(ContactReader):
    def __init__(self, filepath, chromsizes, bins, map=map, **kwargs):
        try:
            import pypairix
        except ImportError:
            raise ImportError("pypairix is required to read pairix-indexed files")
        
        self.P2 = kwargs.pop('P2', 3)
        self._map = map
        self.chromsizes = chromsizes
        self.bins = check_bins(bins, chromsizes)
        self.binsize = get_binsize(bins)
        self.bins_grouped = bins.groupby('chrom', sort=False)
        self.idmap, \
        self.chrom_binoffset, \
        self.chrom_abspos, \
        self.start_abspos = parse_genomic_segmentation(chromsizes, bins)

        # read pair records
        self.filepath = filepath
        self.n_records = None
        self.contigs = list(chromsizes.keys())

        f = pypairix.open(filepath, 'r')
        file_contigs = set(
            itertools.chain.from_iterable(
                [b.split('|') for b in f.get_blocknames()]))

        for chrom in self.contigs:
            if chrom not in file_contigs:
                warnings.warn(
                    "Did not find contig '{}' in contact list file.".format(chrom))
        self.contigs = [ctg for ctg in self.contigs if ctg in file_contigs]

    def __getstate__(self):
        d = self.__dict__.copy()
        d.pop('_map', None)
        return d

    def _size(self, block):
        import pypairix
        f = pypairix.open(self.filepath, 'r')
        chrom1, chrom2 = block
        return sum(1 for line in f.query2D(
            chrom1, 0, self.chromsizes[chrom1],
            chrom2, 0, self.chromsizes[chrom2]))
    
    def size(self):
        if self.n_records is None:
            blocks = itertools.combinations_with_replacement(self.contigs, 2)
            self.n_records = sum(self._map(self._size, blocks))
        return self.n_records
    
    def aggregate(self, chrom1):
        import pypairix
        filepath = self.filepath
        chromsizes = self.chromsizes
        binsize = self.binsize
        idmap = self.idmap
        chrom_binoffset = self.chrom_binoffset
        chrom_abspos = self.chrom_abspos
        start_abspos = self.start_abspos
        P2 = self.P2
        these_bins = self.bins_grouped.get_group(chrom1)

        rows = []

        f = pypairix.open(filepath, 'r')
        accumulator = Counter()
        i = self.contigs.index(chrom1)
        remaining = self.contigs[i:]

        for bin1_id, bin1 in these_bins.iterrows():
            for cid2, chrom2 in enumerate(remaining, i):
                chrom2_size = chromsizes[chrom2]
                for line in f.query2D(
                        chrom1, bin1.start, bin1.end,
                        chrom2, 0, chrom2_size):
                    pos2 = int(line[P2])
                    if binsize is None:
                        lo, hi = chrom_binoffset[cid2], chrom_binoffset[cid2+1]
                        bin2_id = lo + np.searchsorted(
                            start_abspos[lo:hi], 
                            chrom_abspos[cid2] + pos2,
                            side='right') - 1
                    else:
                        bin2_id = chrom_binoffset[cid2] + (pos2 // binsize)
                    accumulator[bin2_id] += 1
            if not accumulator:
                continue
            rows.append(
                pandas.DataFrame({
                    'bin1_id': bin1_id,
                    'bin2_id': list(accumulator.keys()),
                    'count':   list(accumulator.values())},
                    columns=['bin1_id', 'bin2_id', 'count'])
                      .sort_values('bin2_id')
            )
            accumulator.clear()
        
        logger.info(chrom1)

        return pandas.concat(rows, axis=0) if len(rows) else None
    
    def __iter__(self):
        for df in self._map(self.aggregate, list(self.contigs)):
            if df is not None:
                yield {k: v.values for k, v in six.iteritems(df)}


class CoolerAggregator(ContactReader):
    """
    Aggregate contacts from an existing Cooler file.

    """
    def __init__(self, cooler_path, bins, chunksize, cooler_root="/", map=map):
        self._map = map
        
        self.cooler_path = cooler_path
        self.cooler_root = cooler_root
        with h5py.File(self.cooler_path, 'r') as h5:
            grp = h5[cooler_root]
            self._size = grp.attrs['nnz']
            chroms = grp['chroms/name'][:].astype('U')
            lengths = grp['chroms/length'][:]
            self.old_chrom_offset = grp['indexes/chrom_offset'][:]
            self.old_bin1_offset = grp['indexes/bin1_offset'][:]
            self.old_binsize = grp.attrs['bin-size']

        self.new_binsize = get_binsize(bins)
        assert self.new_binsize % self.old_binsize == 0
        self.factor = self.new_binsize // self.old_binsize
        self.chunksize = chunksize
        #self.new_bins = bins
        
        self.chroms = chroms
        self.idmap = pandas.Series(index=chroms, data=range(len(chroms)))
        bin_chrom_ids = self.idmap[bins['chrom']].values
        self.cumul_length = np.r_[0, np.cumsum(lengths)]
        self.abs_start_coords = self.cumul_length[bin_chrom_ids] + bins['start']

        # chrom offset index: chrom_id -> offset in bins
        self.chrom_offset = np.r_[0, 
            np.cumsum(bins.groupby('chrom', sort=False).size())]

    def __getstate__(self):
        d = self.__dict__.copy()
        d.pop('_map', None)
        return d

    def size(self):
        return self._size
    
    def _aggregate(self, span):
        from ..api import Cooler
        lo, hi = span
        looger.info('{} {}'.format(lo, hi))

        try:
            lock.acquire()
            with h5py.File(self.cooler_path, 'r') as h5:
                c = Cooler(h5[self.cooler_root])
                table = c.pixels(join=True, convert_enum=False)
                chunk = table[lo:hi]
                #chunk['chrom1'] = pandas.Categorical(chunk['chrom1'], categories=self.chroms)
                #chunk['chrom2'] = pandas.Categorical(chunk['chrom2'], categories=self.chroms)
        finally:
            lock.release()

        # use the "start" point as anchor for re-binning
        # XXX - alternatives: midpoint anchor, proportional re-binning
        binsize = self.new_binsize
        chrom_offset = self.chrom_offset
        cumul_length = self.cumul_length
        abs_start_coords = self.abs_start_coords

        chrom_id1 = chunk['chrom1'].values  #.cat.codes.values
        chrom_id2 = chunk['chrom2'].values  #.cat.codes.values
        start1 = chunk['start1'].values
        start2 = chunk['start2'].values
        if binsize is None:
            abs_start1 = cumul_length[chrom_id1] + start1
            abs_start2 = cumul_length[chrom_id2] + start2
            chunk['bin1_id'] = np.searchsorted(abs_start_coords, abs_start1, side='right') - 1
            chunk['bin2_id'] = np.searchsorted(abs_start_coords, abs_start2, side='right') - 1
        else:
            rel_bin1 = np.floor(start1/binsize).astype(int)
            rel_bin2 = np.floor(start2/binsize).astype(int)
            chunk['bin1_id'] = chrom_offset[chrom_id1] + rel_bin1
            chunk['bin2_id'] = chrom_offset[chrom_id2] + rel_bin2

        grouped = chunk.groupby(['bin1_id', 'bin2_id'], sort=False)
        return grouped['count'].sum().reset_index()

    def aggregate(self, span):
        try:
            chunk = self._aggregate(span)
        except MemoryError as e:
            raise RuntimeError(str(e))
        return chunk

    def __iter__(self):
        from itertools import chain
        old_chrom_offset = self.old_chrom_offset
        old_bin1_offset = self.old_bin1_offset
        chunksize = self.chunksize
        factor = self.factor
        
        spans = []
        for chrom, i in self.idmap.items():
            # it's important to extract some multiple of `factor` rows at a time
            c0, c1 = old_chrom_offset[i], old_chrom_offset[i+1]
            step = (chunksize // factor) * factor
            edges = np.arange(old_bin1_offset[c0], old_bin1_offset[c1]+step, step)
            edges[-1] = old_bin1_offset[c1]
            spans.append(zip(edges[:-1], edges[1:]))
        spans = list(chain.from_iterable(spans))
        
        for df in self._map(self.aggregate, spans):
            yield {k: v.values for k, v in six.iteritems(df)}


class SparseLoader(ContactReader):
    """
    Load binned contacts from a single 3-column sparse matrix text file.

    """
    def __init__(self, filepath, chunksize):
        # number of lines in file
        p1 = subprocess.Popen(['unpigz',  '-p', '8',  '-c', filepath], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(['wc', '-l'], stdin=p1.stdout, stdout=subprocess.PIPE)
        self.n_records = int(p2.communicate()[0])

        # file iterator
        self.reader = pandas.read_csv(filepath, sep='\t', iterator=True,
            names=['bin1_id', 'bin2_id', 'count'])
        self.chunksize = chunksize

    def size(self):
        return self.n_records

    def __iter__(self):
        while True:
            try:
                data = self.reader.read(self.chunksize)
                yield {k: v.values for k,v in six.iteritems(data)}
            except StopIteration:
                break


class SparseTileLoader(ContactReader):
    """
    Load binned contacts from a collection of 3-column sparse matrix files
    representing contig-contig contact matrix tiles.

    """
    pass


class DenseLoader(ContactReader):
    """
    Load a dense genome-wide numpy array contact matrix.

    """
    def __init__(self, heatmap):
        # TRIU sparsify the matrix
        i, j = np.nonzero(heatmap)
        mask = i <= j
        triu_i, triu_j = i[mask], j[mask]
        self.data = {
            'bin1_id': triu_i,
            'bin2_id': triu_j,
            'count': heatmap[triu_i, triu_j],
        }
        self.nnz = len(triu_i)

    def size(self):
        return self.nnz

    def __iter__(self):
        yield self.data


class DenseTileLoader(ContactReader):
    """
    Load a contact matrix from a collection of dense numpy array contig-contig
    tiles.

    """
    pass

