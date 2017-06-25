# -*- coding: utf-8 -*-
"""
Contact Binners
~~~~~~~~~~~~~~~

Binners are iterators that convert input data of various flavors into a 
properly sorted, chunked stream of binned contacts.

"""
from __future__ import division, print_function
from collections import OrderedDict, Counter
from bisect import bisect_left
from multiprocess import Pool
import itertools
import warnings
import sys
import six

#from pandas.algos import is_lexsorted
import numpy as np
import pandas
import h5py

from . import get_logger
from .util import rlencode, get_binsize, parse_region
from .io import parse_cooler_uri
from .tools import lock, partition


logger = get_logger()


class ContactBinner(object):
    """
    Base class for iterable contact binners.

    """
    def __getstate__(self):
        d = self.__dict__.copy()
        d.pop('_map', None)
        return d

    def size(self):
        """ Returns the total number of contacts.  **DEPRECATED** """
        raise NotImplementedError

    def __iter__(self):
        """ Iterator over chunks of binned contacts (i.e., nonzero pixels)

        Chunks are expected to have the following format:

        * dict of 1D arrays
        * keys `bin1_id`, `bin2_id`, `count`
        * arrays lexically sorted by `bin_id` then `bin2_id`

        """
        raise NotImplementedError


def check_bins(bins, chromsizes):
    is_cat = pandas.api.types.is_categorical(bins['chrom'])
    bins = bins.copy()
    if not is_cat:
        bins['chrom'] = pandas.Categorical(
            bins.chrom, 
            categories=list(chromsizes.index), 
            ordered=True)
    else:
        assert (bins['chrom'].cat.categories == chromsizes.index).all()

    return bins


def balanced_partition(gs, n_chunk_max, file_contigs, loadings=None):
    n_bins = len(gs.bins)
    grouped = gs._bins_grouped
    
    chrom_nbins = grouped.size()
    if loadings is None:
        loadings = chrom_nbins
    chrmax = loadings.argmax()
    loadings = loadings / loadings.loc[chrmax]
    const = chrom_nbins.loc[chrmax] / n_chunk_max

    granges = []
    for chrom, group in grouped:
        if chrom not in file_contigs:
            continue
        clen = gs.chromsizes[chrom]
        step = int(np.ceil(const / loadings.loc[chrom]))        
        anchors = group.start.values[::step]
        if anchors[-1] != clen:
            anchors = np.r_[anchors, clen]
        granges.extend( (chrom, start, end) 
            for start, end in zip(anchors[:-1], anchors[1:]))
    return granges


class GenomeSegmentation(object):
    def __init__(self, chromsizes, bins):
        bins = check_bins(bins, chromsizes)
        self._bins_grouped = bins.groupby('chrom', sort=False)
        nbins_per_chrom = self._bins_grouped.size().values

        self.chromsizes = chromsizes
        self.binsize = get_binsize(bins)
        self.contigs = list(chromsizes.keys())
        self.bins = bins
        self.idmap = pandas.Series(
            index=chromsizes.keys(), 
            data=range(len(chromsizes)))
        self.chrom_binoffset = np.r_[0, np.cumsum(nbins_per_chrom)]
        self.chrom_abspos = np.r_[0, np.cumsum(chromsizes.values)]
        self.start_abspos = (self.chrom_abspos[bins['chrom'].cat.codes] + 
                             bins['start'].values)
    
    def fetch(self, region):
        chrom, start, end = parse_region(region, self.chromsizes)
        result = self._bins_grouped.get_group(chrom)
        if start > 0 or end < self.chromsizes[chrom]:
            lo = result['end'].values.searchsorted(start, side='right')
            hi = lo + result['start'].values[lo:].searchsorted(end, side='left')
            result = result.iloc[lo:hi]
        return result


class HDF5Aggregator(ContactBinner):
    """
    Aggregate contacts from a hiclib-style HDF5 contacts file.

    """
    def __init__(self, h5pairs, chromsizes, bins, chunksize, **kwargs):
        self.h5 = h5pairs
        self.C1 = kwargs.pop('C1', 'chrms1')
        self.P1 = kwargs.pop('P1', 'cuts1')
        self.C2 = kwargs.pop('C2', 'chrms2')
        self.P2 = kwargs.pop('P2', 'cuts2')
        self.gs = GenomeSegmentation(chromsizes, bins)
        self.chunksize = chunksize
        self.partition = self._index_chroms()

    def _index_chroms(self):
        # index extents of chromosomes on first axis of contact list
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
        chunksize = self.chunksize
        bins = self.gs.bins
        binsize = self.gs.binsize
        chrom_binoffset = self.gs.chrom_binoffset
        chrom_abspos = self.gs.chrom_abspos
        start_abspos = self.gs.start_abspos
        cid = self.gs.idmap[chrom]

        chrom_lo, chrom_hi = self.partition.get(cid, (-1, -1))
        lo = chrom_lo
        hi = lo
        while hi < chrom_hi:
            # update `hi` to make sure our selection doesn't split a bin1
            lo, hi = hi, min(hi + chunksize, chrom_hi)
            abspos = chrom_abspos[cid] + h5pairs[P1][hi - 1]
            bin_id = int(np.searchsorted(
                start_abspos, abspos, side='right')) - 1
            bin_end = bins['end'][bin_id]
            hi = bisect_left(h5pairs[P1], bin_end, lo, chrom_hi)
            if lo == hi:
                hi = chrom_hi

            logger.info('{} {}'.format(lo, hi))

            # load chunk and assign bin IDs to each read side
            table = self._load_chunk(lo, hi)
            abspos1 = chrom_abspos[h5pairs[C1][lo:hi]] + h5pairs[P1][lo:hi]
            abspos2 = chrom_abspos[h5pairs[C2][lo:hi]] + h5pairs[P2][lo:hi]
            if np.any(abspos1 > abspos2):
                raise ValueError(
                    "Found a read pair that maps to the lower triangle of the "
                    "contact map (side1 > side2). Check that the provided "
                    "chromosome ordering and read pair file are consistent "
                    "such that all pairs map to the upper triangle with "
                    "respect to the given chromosome ordering.")

            if binsize is None:
                table['bin1_id'] = np.searchsorted(
                    start_abspos, abspos1, side='right') - 1
                table['bin2_id'] = np.searchsorted(
                    start_abspos, abspos2, side='right') - 1
            else:
                rel_bin1 = np.floor(table['cut1']/binsize).astype(int)
                rel_bin2 = np.floor(table['cut2']/binsize).astype(int)
                table['bin1_id'] = (
                    chrom_binoffset[table['chrom_id1'].values] + rel_bin1)
                table['bin2_id'] = (
                    chrom_binoffset[table['chrom_id2'].values] + rel_bin2)

            # reduce
            gby = table.groupby(['bin1_id', 'bin2_id'])
            agg = (gby['chrom_id1'].count()
                                   .reset_index()
                                   .rename(columns={'chrom_id1': 'count'}))
            yield agg

    def size(self):
        return len(self.h5['chrms1'])

    def __iter__(self):
        for chrom in self.gs.contigs:
            for df in self.aggregate(chrom):
                yield {k: v.values for k, v in six.iteritems(df)}


class TabixAggregator(ContactBinner):
    """
    Aggregate contacts from a sorted, BGZIP-compressed and tabix-indexed
    tab-delimited text file.

    """
    def __init__(self, filepath, chromsizes, bins, map=map, n_chunks=1, **kwargs):
        try:
            import pysam
        except ImportError:
            raise ImportError("pysam is required to read tabix files")

        import dill
        import pickle
        dill.settings['protocol'] = pickle.HIGHEST_PROTOCOL
        
        self._map = map
        self.n_chunks = n_chunks
        self.C2 = kwargs.pop('C2', 3)
        self.P2 = kwargs.pop('P2', 4)
        
        # all requested contigs will be placed in the output matrix
        self.gs = GenomeSegmentation(chromsizes, bins)
        
        # find available contigs in the contact list
        self.filepath = filepath
        self.n_records = None
        with pysam.TabixFile(filepath, 'r', encoding='ascii') as f:
            try:
                self.file_contigs = [c.decode('ascii') for c in f.contigs]
            except AttributeError:
                self.file_contigs = f.contigs
        
        # warn about requested contigs not seen in the contact list
        for chrom in self.gs.contigs:
            if chrom not in self.file_contigs:
                warnings.warn(
                    "Did not find contig " +
                    " '{}' in contact list file.".format(chrom))
        
        warnings.warn(
            "NOTE: When using the Tabix aggregator, make sure the order of "
            "chromosomes in the provided chromsizes agrees with the chromosome "
            "ordering of read ends in the contact list file.")

    def aggregate(self, grange):
        chrom1, start, end = grange
        import pysam
        filepath = self.filepath
        binsize = self.gs.binsize
        idmap = self.gs.idmap
        chromsizes = self.gs.chromsizes
        chrom_binoffset = self.gs.chrom_binoffset
        chrom_abspos = self.gs.chrom_abspos
        start_abspos = self.gs.start_abspos
        C2 = self.C2
        P2 = self.P2

        logger.info('Binning {}:{}-{}|*'.format(chrom1, start, end))
        
        these_bins = self.gs.fetch((chrom1, start, end))
        rows = []
        with pysam.TabixFile(filepath, 'r', encoding='ascii') as f:
            parser = pysam.asTuple()
            accumulator = Counter()
            
            for bin1_id, bin1 in these_bins.iterrows():
                for line in f.fetch(chrom1, bin1.start, bin1.end,
                                    parser=parser):
                    chrom2 = line[C2]
                    pos2 = int(line[P2])
                    
                    try:
                        cid2 = idmap[chrom2]
                    except KeyError:
                        # this chrom2 is not requested
                        continue
                    
                    if binsize is None:
                        lo = chrom_binoffset[cid2]
                        hi = chrom_binoffset[cid2 + 1]
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
        
        logger.info('Finished {}:{}-{}|*'.format(chrom1, start, end))

        return pandas.concat(rows, axis=0) if len(rows) else None
    
    def __iter__(self):
        granges = balanced_partition(self.gs, self.n_chunks, self.file_contigs)
        for df in self._map(self.aggregate, granges):
            if df is not None:
                yield {k: v.values for k, v in six.iteritems(df)}


class PairixAggregator(ContactBinner):
    """
    Aggregate contacts from a sorted, BGZIP-compressed and pairix-indexed
    tab-delimited text file.

    """
    def __init__(self, filepath, chromsizes, bins, map=map, n_chunks=1, **kwargs):
        try:
            import pypairix
        except ImportError:
            raise ImportError(
                "pypairix is required to read pairix-indexed files")

        import dill
        import pickle
        dill.settings['protocol'] = pickle.HIGHEST_PROTOCOL
        
        self._map = map
        self.n_chunks = n_chunks
        f = pypairix.open(filepath, 'r')
        self.C1 = f.get_chr1_col()
        self.C2 = f.get_chr2_col()
        self.P1 = f.get_startpos1_col()
        self.P2 = f.get_startpos2_col()
        self.file_contigs = set(
            itertools.chain.from_iterable(
                [b.split('|') for b in f.get_blocknames()]))
        
        # all requested contigs will be placed in the output matrix
        self.gs = GenomeSegmentation(chromsizes, bins)
        
        # find available contigs in the contact list
        self.filepath = filepath
        self.n_records = None
       
        # warn about requested contigs not seen in the contact list
        for chrom in self.gs.contigs:
            if chrom not in self.file_contigs:
                warnings.warn(
                    "Did not find contig " +
                    " '{}' in contact list file.".format(chrom))

    def aggregate(self, grange):
        chrom1, start, end = grange
        import pypairix
        filepath = self.filepath
        binsize = self.gs.binsize
        chromsizes = self.gs.chromsizes
        chrom_binoffset = self.gs.chrom_binoffset
        chrom_abspos = self.gs.chrom_abspos
        start_abspos = self.gs.start_abspos
        C1 = self.C1
        C2 = self.C2
        P1 = self.P1
        P2 = self.P2
        
        logger.info('Binning {}:{}-{}|*'.format(chrom1, start, end))

        f = pypairix.open(filepath, 'r')
        these_bins = self.gs.fetch((chrom1, start, end))
        remaining_chroms = self.gs.idmap[chrom1:]
        cid1 = self.gs.idmap[chrom1]

        accumulator = Counter()
        rows = []
        for bin1_id, bin1 in these_bins.iterrows():
            
            for chrom2, cid2 in six.iteritems(remaining_chroms):
                
                chrom2_size = chromsizes[chrom2]

                if chrom1 != chrom2 and f.exists2(chrom2, chrom1):  # flipped
                    iterator = f.query2D(chrom2, 0, chrom2_size, 
                                         chrom1, bin1.start, bin1.end)
                    pos2_col = P1
                else:
                    iterator = f.query2D(chrom1, bin1.start, bin1.end, 
                                         chrom2, 0, chrom2_size)
                    pos2_col = P2

                for line in iterator:
                    
                    pos2 = int(line[pos2_col])

                    if binsize is None:
                        lo = chrom_binoffset[cid2]
                        hi = chrom_binoffset[cid2 + 1]
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
        
        logger.info('Finished {}:{}-{}|*'.format(chrom1, start, end))

        return pandas.concat(rows, axis=0) if len(rows) else None
    
    def __iter__(self):
        granges = balanced_partition(self.gs, self.n_chunks, self.file_contigs)
        for df in self._map(self.aggregate, granges):
            if df is not None:
                yield {k: v.values for k, v in six.iteritems(df)}


class CoolerAggregator(ContactBinner):
    """
    Aggregate contacts from an existing Cooler file.

    """
    def __init__(self, source_uri, bins, chunksize, batchsize, map=map):
        from cooler.api import Cooler
        self._map = map
        self.source_uri = source_uri
        self.chunksize = chunksize
        self.batchsize = batchsize

        clr = Cooler(source_uri)
        self._size = clr.info['nnz']
        self.old_binsize = clr.binsize
        self.old_chrom_offset = clr._load_dset('indexes/chrom_offset')
        self.old_bin1_offset = clr._load_dset('indexes/bin1_offset')
        self.gs = GenomeSegmentation(clr.chromsizes, bins)
        self.new_binsize = get_binsize(bins)
        assert self.new_binsize % self.old_binsize == 0
        self.factor = self.new_binsize // self.old_binsize
    
    def _aggregate(self, span):
        from cooler.api import Cooler
        lo, hi = span

        clr = Cooler(self.source_uri)
        # convert_enum=False returns chroms as raw ints
        table = clr.pixels(join=True, convert_enum=False)
        chunk = table[lo:hi]
        logger.info('{} {}'.format(lo, hi))

        # use the "start" point as anchor for re-binning
        # XXX - alternatives: midpoint anchor, proportional re-binning
        binsize = self.gs.binsize
        chrom_binoffset = self.gs.chrom_binoffset
        chrom_abspos = self.gs.chrom_abspos
        start_abspos = self.gs.start_abspos

        chrom_id1 = chunk['chrom1'].values
        chrom_id2 = chunk['chrom2'].values
        start1 = chunk['start1'].values
        start2 = chunk['start2'].values
        if binsize is None:
            abs_start1 = chrom_abspos[chrom_id1] + start1
            abs_start2 = chrom_abspos[chrom_id2] + start2
            chunk['bin1_id'] = np.searchsorted(
                start_abspos, 
                abs_start1, 
                side='right') - 1
            chunk['bin2_id'] = np.searchsorted(
                start_abspos, 
                abs_start2, 
                side='right') - 1
        else:
            rel_bin1 = np.floor(start1/binsize).astype(int)
            rel_bin2 = np.floor(start2/binsize).astype(int)
            chunk['bin1_id'] = chrom_binoffset[chrom_id1] + rel_bin1
            chunk['bin2_id'] = chrom_binoffset[chrom_id2] + rel_bin2

        grouped = chunk.groupby(['bin1_id', 'bin2_id'], sort=False)
        return grouped['count'].sum().reset_index()

    def aggregate(self, span):
        try:
            chunk = self._aggregate(span)
        except MemoryError as e:
            raise RuntimeError(str(e))
        return chunk

    def __iter__(self):
        old_chrom_offset = self.old_chrom_offset
        old_bin1_offset = self.old_bin1_offset
        chunksize = self.chunksize
        batchsize = self.batchsize
        factor = self.factor
        
        # Partition pixels into chunks, respecting chrom1 boundaries
        spans = []
        for chrom, i in six.iteritems(self.gs.idmap):
            # it's important to extract some multiple of `factor` rows at a time
            c0 = old_chrom_offset[i]
            c1 = old_chrom_offset[i + 1]
            step = (chunksize // factor) * factor
            edges = np.arange(
                old_bin1_offset[c0], 
                old_bin1_offset[c1] + step, 
                step)
            edges[-1] = old_bin1_offset[c1]
            spans.append(zip(edges[:-1], edges[1:]))
        spans = list(itertools.chain.from_iterable(spans))
        
        # Process batches of k chunks at a time, then yield the results
        for i in range(0, len(spans), batchsize):
            try:
                lock.acquire()
                results = self._map(self.aggregate, spans[i:i+batchsize])
            finally:
                lock.release()
            for df in results:
                yield {k: v.values for k, v in six.iteritems(df)}


class CoolerMerger(ContactBinner):
    """
    Merge (i.e. sum) multiple cooler matrices with identical axes.

    """
    def __init__(self, coolers, chunksize, **kwargs):
        self.coolers = list(coolers)
        self.chunksize = chunksize

        binsize = coolers[0].binsize
        if binsize is not None:
            if len(set(c.binsize for c in coolers)) > 1:
                raise ValueError("Coolers must have the same resolution")
            chromsizes = coolers[0].chromsizes
            for i in range(1, len(coolers)):
                if not np.all(coolers[i].chromsizes == chromsizes):
                    raise ValueError("Coolers must have the same chromosomes")
        else:
            bins = coolers[0].bins()[['chrom', 'start', 'end']][:]
            for i in range(1, len(coolers)):
                if not np.all(
                    coolers[i].bins()[['chrom', 'start', 'end']][:] == bins):
                    raise ValueError("Coolers must have same bin structure")

    def __iter__(self):
        chunksize = self.chunksize
        indexes = [c._load_dset('indexes/bin1_offset') for c in self.coolers]
        nnzs = [len(c.pixels()) for c in self.coolers]
        logger.info('nnzs: {}'.format(nnzs))

        lo = 0
        starts = [0] * len(self.coolers)
        while True:
            hi = max(bisect_left(o[:-1], min(start + chunksize, nnz), lo=lo) 
                                 for start, nnz, o in zip(starts, nnzs, indexes))
            if hi == lo:
                break
            stops = [o[hi] for o in indexes]
            logger.info('current: {}'.format(stops))
            
            combined = pandas.concat(
                [c.pixels()[start:stop] 
                    for c, start, stop in zip(self.coolers, starts, stops)],
                axis=0,
                ignore_index=True)

            df = (combined.groupby(['bin1_id', 'bin2_id'], sort=True)
                          .aggregate({'count': np.sum})
                          .reset_index())
            yield {k: v.values for k, v in six.iteritems(df)}

            lo = hi
            starts = stops


class BedGraph2DLoader(ContactBinner):
    """
    Contact iterator for a sparse tsv Hi-C matrix with fields:
        "chrom1, start1, end1, chrom2, start2, end2, count"

    """
    FIELD_NUMBERS = OrderedDict([
        ('chrom1', 0),
        ('start1', 1),
        ('end1', 2),
        ('chrom2', 3),
        ('start2', 4),
        ('end2', 5),
        ('count', 6),
    ])

    FIELD_DTYPES = {
        'chrom1': object,
        'start1': int,
        'end1': int,
        'chrom2': object,
        'start2': int,
        'end2': int,
        'count': int,
    }

    def __init__(self, filepath, chromsizes, bins, field_numbers=None, 
                 field_dtypes=None):
        try:
            import pypairix
        except ImportError:
            raise ImportError(
                "pypairix is required to read pairix-indexed files")
        
        self._map = map
        self.filepath = filepath
        f = pypairix.open(self.filepath, 'r')
        self.file_contigs = set(
            itertools.chain.from_iterable(
                [b.split('|') for b in f.get_blocknames()]))
        
        # all requested contigs will be placed in the output matrix
        self.gs = GenomeSegmentation(chromsizes, bins)
        self.bins = self.gs.bins.copy()
        self.bins['chrom'] = self.bins['chrom'].astype(object)
        self.bins['bin'] = self.bins.index
       
        # warn about requested contigs not seen in the contact list
        for chrom in self.gs.contigs:
            if chrom not in self.file_contigs:
                warnings.warn(
                    "Did not find contig " +
                    " '{}' in bg2 file.".format(chrom))

        # Assign the column numbers
        self.field_numbers = self.FIELD_NUMBERS.copy()
        if field_numbers is not None:
            self.field_numbers.update(field_numbers)
        self.columns = list(self.field_numbers.keys())
        self.usecols = list(self.field_numbers.values())
        
        # Assign the column dtypes. Assume additional value fields are float.
        self.out_columns = ['bin1_id', 'bin2_id', 'count']
        self.dtypes = self.FIELD_DTYPES.copy()
        for col in self.columns:
            if col not in self.dtypes:
                self.out_columns.append(col)
                self.dtypes[col] = float

        # Override defaults
        if field_dtypes is not None:
            self.dtypes.update(field_dtypes)
    
    def aggregate(self, chrom1):
        import pypairix

        f = pypairix.open(self.filepath, 'r')
        cid1 = self.gs.idmap[chrom1]
        chromsizes = self.gs.chromsizes
        these_bins = self.gs.fetch(chrom1)
        remaining_chroms = self.gs.idmap[chrom1:]
        c1 = self.field_numbers['chrom1']
        c2 = self.field_numbers['chrom2']
        s1 = self.field_numbers['start1']
        s2 = self.field_numbers['start2']
        e1 = self.field_numbers['end1']
        e2 = self.field_numbers['end2']

        # read contact matrix one row at a time
        logger.info(chrom1)
        lines = []
        for bin1_id, bin1 in these_bins.iterrows():
            for chrom2, cid2 in six.iteritems(remaining_chroms):
                chrom2_size = chromsizes[chrom2]
                if chrom1 != chrom2 and f.exists2(chrom2, chrom1):  # flipped block
                    q = []
                    for line in f.query2D(
                                    chrom2, 0, chrom2_size, 
                                    chrom1, bin1.start, bin1.end):
                        line[c1], line[c2] = line[c2], line[c1]
                        line[s1], line[s2] = line[s2], line[s1]
                        line[e1], line[e2] = line[e2], line[e1]
                        q.append(line)
                else:
                    q = list(line for line in 
                             f.query2D(chrom1, bin1.start, bin1.end,
                                       chrom2, 0, chrom2_size))
                lines.extend(q)

        df = pandas.DataFrame(lines)
        df = df[self.usecols]
        df.columns = self.columns
        for col, dtype in self.dtypes.items():
            df[col] = df[col].astype(dtype)

        # assign bin IDs from bin table
        df = (df.merge(self.bins, 
                       left_on=['chrom1', 'start1', 'end1'], 
                       right_on=['chrom', 'start', 'end'])
                .merge(self.bins, 
                       left_on=['chrom2', 'start2', 'end2'], 
                       right_on=['chrom', 'start', 'end'], 
                       suffixes=('1', '2'))
                .rename(columns={'bin1': 'bin1_id', 
                                 'bin2': 'bin2_id'}))
        df = (df[self.out_columns]
                .sort_values(['bin1_id', 'bin2_id']))
        return df
    
    def __iter__(self):
        chroms = [ctg for ctg in self.gs.contigs 
                       if ctg in self.file_contigs]

        for df in self._map(self.aggregate, chroms):
            if df is not None:
                yield {k: v.values for k, v in six.iteritems(df)}


class SparseLoader(ContactBinner):
    """
    Load binned contacts from a single 3-column sparse matrix text file.

    """
    FIELD_NUMBERS = OrderedDict([
        ('bin1_id', 0),
        ('bin2_id', 1),
        ('count', 2),
    ])

    FIELD_DTYPES = OrderedDict([
        ('bin1_id', int),
        ('bin2_id', int),
        ('count', int),
    ])

    def __init__(self, filepath, bins, chunksize, field_numbers=None, 
                 field_dtypes=None):
        """
        Parameters
        ----------
        filepath : str
            Path to tsv file
        chunksize : number of rows of the matrix file to read at a time

        """
        self._map = map
        self.filepath = filepath
        self.chunksize = chunksize
        self.n_bins = len(bins)

        # Assign the column numbers
        self.field_numbers = self.FIELD_NUMBERS.copy()
        if field_numbers is not None:
            self.field_numbers.update(field_numbers)
        self.columns = list(self.field_numbers.keys())
        self.usecols = list(self.field_numbers.values())
        
        # Assign the column dtypes. Assume additional value fields are float.
        self.out_columns = ['bin1_id', 'bin2_id', 'count']
        self.dtypes = self.FIELD_DTYPES.copy()
        for col in self.columns:
            if col not in self.dtypes:
                self.out_columns.append(col)
                self.dtypes[col] = float

        # Override defaults
        if field_dtypes is not None:
            self.dtypes.update(field_dtypes)

    def __iter__(self):
        n_bins = self.n_bins
        
        iterator = pandas.read_csv(
            self.filepath, 
            sep='\t', 
            iterator=True,
            comment='#',
            chunksize=self.chunksize,
            usecols=self.usecols,
            names=self.columns,
            dtype=self.dtypes)

        for chunk in iterator:
            if np.any(chunk['bin1_id'] > chunk['bin2_id']):
                raise ValueError("Found bin1_id greater than bin2_id")
            if (np.any(chunk['bin1_id'] >= n_bins) or 
                np.any(chunk['bin2_id'] >= n_bins)):
                raise ValueError(
                    "Found a bin ID that exceeds the declared number of bins. " 
                    "Check whether your bin table is correct.")
            yield {k: v.values for k,v in six.iteritems(chunk)}


class SparseBlockLoader(ContactBinner):
    """

    """
    def __init__(self, chromsizes, bins, mapping, chunksize):
        bins = check_bins(bins, chromsizes)
        self.bins = bins
        self.chromosomes = list(chromsizes.index)
        self.chunksize = chunksize
        n_chroms = len(chromsizes)
        n_bins = len(bins)
        
        chrom_ids = bins['chrom'].cat.codes
        self.offsets = np.zeros(n_chroms + 1, dtype=int)
        curr_val = 0
        for start, length, value in zip(*rlencode(chrom_ids)):
            self.offsets[curr_val:value + 1] = start
            curr_val = value + 1
        self.offsets[curr_val:] = n_bins
        
        self.mapping = mapping
        
    def select_block(self, chrom1, chrom2):
        try:
            block = self.mapping[chrom1, chrom2]
        except KeyError:
            try:
                block = self.mapping[chrom2, chrom1].T
            except KeyError:
                warning.warn(
                    "Block for {{{}, {}}} not found".format(chrom1, chrom2))
                raise
        return block

    def __iter__(self):
        n_bins = len(self.bins)
        chromosomes = self.chromosomes
        
        for cid1, chrom1 in enumerate(chromosomes):
            offset = self.offsets[cid1]
            chrom1_nbins = self.offsets[cid1 + 1] - offset
            spans = partition(0, chrom1_nbins, self.chunksize)
            
            for lo, hi in spans:
                chunks = []
                for chrom2 in chromosomes[cid1:]:
                    try:
                        block = self.select_block(chrom1, chrom2)
                    except KeyError:
                        continue
                    chunks.append(block.tocsr()[lo:hi, :])
                X = scipy.sparse.hstack(chunks).tocsr().tocoo()
                
                i, j, v = X.row, X.col, X.data
                mask = (offset + i) <= (offset + j)
                triu_i, triu_j, triu_v = i[mask], j[mask], v[mask]

                yield {
                    'bin1_id': offset + triu_i,
                    'bin2_id': offset + triu_j,
                    'count': triu_v,
                }


class ArrayLoader(ContactBinner):
    """
    Load a dense genome-wide numpy array contact matrix.
    Works with array-likes such as h5py.Dataset and memmapped arrays

    See Also
    --------
    numpy.save, numpy.load (mmap_mode)

    """
    def __init__(self, bins, array, chunksize):
        if len(bins) != array.shape[0]:
            raise ValueError("Number of bins must equal the dimenion of the matrix")
        self.array = array
        self.chunksize = chunksize
    
    def __iter__(self):
        n_bins = self.array.shape[0]
        spans = partition(0, n_bins, self.chunksize)
        
        # TRIU sparsify the matrix
        for lo, hi in spans:
            X = self.array[lo:hi, :]
            i, j = np.nonzero(X)
            
            mask = (lo + i) <= j
            triu_i, triu_j = i[mask], j[mask]
            
            yield {
                'bin1_id': lo + triu_i,
                'bin2_id': triu_j,
                'count': X[triu_i, triu_j],
            }


class ArrayBlockLoader(ContactBinner):
    """
    
    """
    def __init__(self, chromsizes, bins, mapping, chunksize):
        bins = check_bins(bins, chromsizes)
        self.bins = bins
        self.chromosomes = list(chromsizes.index)
        self.chunksize = chunksize
        n_chroms = len(chromsizes)
        n_bins = len(bins)
        
        chrom_ids = bins['chrom'].cat.codes
        self.offsets = np.zeros(n_chroms + 1, dtype=int)
        curr_val = 0
        for start, length, value in zip(*rlencode(chrom_ids)):
            self.offsets[curr_val:value + 1] = start
            curr_val = value + 1
        self.offsets[curr_val:] = n_bins
        
        self.mapping = mapping
        
    def select_block(self, chrom1, chrom2):
        try:
            block = self.mapping[chrom1, chrom2]
        except KeyError:
            try:
                block = self.mapping[chrom2, chrom1].T
            except KeyError:
                warning.warn(
                    "Block for {{{}, {}}} not found".format(chrom1, chrom2))
                raise
        return block

    def __iter__(self):
        n_bins = len(self.bins)
        chromosomes = self.chromosomes
        
        for cid1, chrom1 in enumerate(chromosomes):
            offset = self.offsets[cid1]
            chrom1_nbins = self.offsets[cid1 + 1] - offset
            spans = partition(0, chrom1_nbins, self.chunksize)
            
            for lo, hi in spans:
                chunks = []
                for chrom2 in chromosomes[cid1:]:
                    try:
                        block = self.select_block(chrom1, chrom2)
                    except KeyError:
                        continue
                    chunks.append(block[lo:hi, :])
                X = np.concatenate(chunks, axis=1)
                
                i, j = np.nonzero(X)
                mask = (offset + i) <= (offset + j)
                triu_i, triu_j = i[mask], j[mask]

                yield {
                    'bin1_id': offset + triu_i,
                    'bin2_id': offset + triu_j,
                    'count': X[triu_i, triu_j],
                }
