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
import subprocess
import json
import sys
import six

from pandas.algos import is_lexsorted
import numpy as np
import pandas
import h5py

from ..util import rlencode, get_binsize


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


class HDF5Aggregator(ContactReader):
    """
    Aggregate contacts from a hiclib-style HDF5 contacts file.

    """
    def __init__(self, h5pairs, chromsizes, bins, chunksize):
        self.h5 = h5pairs
        self.bins = bins
        self.binsize = get_binsize(bins)
        self.n_contacts = len(self.h5['chrms1'])
        self.n_bins = len(bins)
        self.chunksize = chunksize
        # convert genomic coords of bin starts to absolute
        self.idmap = pandas.Series(index=chromsizes.keys(), data=range(len(chromsizes)))
        bin_chrom_ids = self.idmap[bins['chrom']].values
        self.cumul_length = np.r_[0, np.cumsum(chromsizes)]
        self.abs_start_coords = self.cumul_length[bin_chrom_ids] + bins['start']
        # chrom offset index: chrom_id -> offset in bins
        chrom_nbins =  bins.groupby(bin_chrom_ids, sort=False).size()
        self.chrom_offset = np.r_[0, np.cumsum(chrom_nbins)]
        # index extents of chromosomes on first axis of contact list
        self.partition = self._index_chroms()

    def _index_chroms(self):
        starts, lengths, values = rlencode(self.h5['chrms1'], self.chunksize)
        if len(set(values)) != len(values):
            raise ValueError("Read pair coordinates are not sorted on the first axis")
        return dict(zip(values, zip(starts, starts + lengths)))

    def _load_chunk(self, lo, hi):
        data = OrderedDict([
            ('chrom_id1', self.h5['chrms1'][lo:hi]),
            ('cut1', self.h5['cuts1'][lo:hi]),
            ('chrom_id2', self.h5['chrms2'][lo:hi]),
            ('cut2', self.h5['cuts2'][lo:hi]),
        ])
        return pandas.DataFrame(data)

    def _iterchunks(self, chrom):
        h5pairs = self.h5
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
            abs_pos = cumul_length[cid] + h5pairs['cuts1'][hi-1]
            i = int(np.searchsorted(abs_start_coords, abs_pos, side='right')) - 1
            bin_end = bins['end'][i]
            hi = bisect_left(h5pairs['cuts1'], bin_end, lo, chrom_hi)
            print(lo, hi)  # flush=True

            # assign bins to reads
            table = self._load_chunk(lo, hi)
            abs_pos1 = cumul_length[h5pairs['chrms1'][lo:hi]] + h5pairs['cuts1'][lo:hi]
            abs_pos2 = cumul_length[h5pairs['chrms2'][lo:hi]] + h5pairs['cuts2'][lo:hi]
            if np.any(abs_pos1 > abs_pos2):
                raise ValueError(
                    "Found a read pair that maps to the lower triangle of the contact map (side1 > side2). "
                    "Check that the provided chromsome ordering and read pair file are consistent "
                    "such that all pairs map to the upper triangle with respect to the given "
                    "chromosome ordering.")

            if binsize is None:
                table['bin1_id'] = np.searchsorted(abs_start_coords, abs_pos1, side='right') - 1
                table['bin2_id'] = np.searchsorted(abs_start_coords, abs_pos2, side='right') - 1
            else:
                rel_bin1 = np.floor(table['cut1']/binsize).astype(int)
                rel_bin2 = np.floor(table['cut2']/binsize).astype(int)
                table['bin1_id'] = chrom_offset[table['chrom_id1'].values] + rel_bin1
                table['bin2_id'] = chrom_offset[table['chrom_id2'].values] + rel_bin2

            # reduce
            gby = table.groupby(['bin1_id', 'bin2_id'])
            agg = (gby['chrom_id1'].count()
                                   .reset_index()
                                   .rename(columns={'chrom_id1': 'count'}))
            yield {k:v.values for k,v in six.iteritems(agg)}

    def size(self):
        return len(self.h5['chrms1'])

    def __iter__(self):
        for chrom in self.idmap.keys():
            for chunk in self._iterchunks(chrom):
                yield chunk


class TabixAggregator(ContactReader):
    """
    Aggregate contacts from a sorted, BGZIP-compressed and tabix-indexed
    tab-delimited text file.

    """
    def __init__(self, filepath, chromsizes, bins):
        try:
            import pysam
        except ImportError:
            raise ImportError("pysam is required to read tabix files")
        n_bins = len(bins)
        self.idmap = pandas.Series(index=chromsizes.keys(), data=range(len(chromsizes)))
        self.bins = bins
        self.binsize = get_binsize(bins)
        self.pairsfile = pysam.TabixFile(filepath, 'r')
        self.parser = pysam.asTuple()
        # number of lines in file
        p1 = subprocess.Popen(['unpigz',  '-p', '8',  '-c', filepath], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(['wc', '-l'], stdin=p1.stdout, stdout=subprocess.PIPE)
        self.n_records = int(p2.communicate()[0])
        # convert genomic coords of bin starts to absolute
        self.idmap = pandas.Series(index=chromsizes.keys(), data=range(len(chromsizes)))
        bin_chrom_ids = self.idmap[bins['chrom']].values
        self.cumul_length = np.r_[0, np.cumsum(chromsizes)]
        self.abs_start_coords = self.cumul_length[bin_chrom_ids] + bins['start']
        # chrom offset index: chrom_id -> offset in bins
        chrom_nbins =  bins.groupby(bin_chrom_ids, sort=False).size()
        self.chrom_offset = np.r_[0, np.cumsum(chrom_nbins)]

    def _iterchunks(self, chrom):
        pairsfile = self.pairsfile
        chrom_offset = self.chrom_offset
        idmap = self.idmap
        parser = self.parser
        binsize = self.binsize
        bins = self.bins[self.bins.chrom==chrom]

        print(chrom)  # flush=True

        for bin1_id, bin in bins.iterrows():
            accumulator = Counter()
            for record in pairsfile.fetch(chrom, bin.start, bin.end, parser=parser):
                chrom2 = record[3]
                pos2 = int(record[4])
                bin2_id = chrom_offset[idmap[chrom2]] + pos2 // binsize
                accumulator[bin2_id] += 1
            if len(accumulator) == 0:
                continue

            agg = pandas.DataFrame({
                'bin1_id': bin1_id,
                'bin2_id': list(accumulator.keys()),
                'count': list(accumulator.values()),
            }, columns=['bin1_id', 'bin2_id', 'count'])
            agg = agg.sort_values('bin2_id')

            yield {k: v.values for k,v in six.iteritems(agg)}

    def size(self):
        return self.n_records

    def __iter__(self):
        for chrom in self.idmap.keys():
            for chunk in self._iterchunks(chrom):
                yield chunk


class CoolerAggregator(ContactReader):
    """
    Aggregate contacts from an existing Cooler file.

    """
    def __init__(self, cool, bins, chunksize):
        self.cool = cool
        self.bins = bins
        self.binsize = get_binsize(bins)
        self.chunksize = chunksize
        
        # convert genomic coords of bin starts to absolute
        chroms = self.cool.chroms()['name'][:].values
        lengths = self.cool.chroms()['length'][:].values
        self.idmap = pandas.Series(index=chroms, data=range(len(chroms)))
        bin_chrom_ids = self.idmap[bins['chrom']].values
        self.cumul_length = np.r_[0, np.cumsum(lengths)]
        self.abs_start_coords = self.cumul_length[bin_chrom_ids] + bins['start']
        
        # chrom offset index: chrom_id -> offset in bins
        self.chrom_offset = np.r_[0, np.cumsum(bins.groupby('chrom', sort=False).size())]

    def _aggregate(self, chunk):
        # use the "start" point as anchor for re-binning
        # XXX - alternatives: midpoint anchor, proportional re-binning
        bins = self.bins
        binsize = self.binsize
        chrom_offset = self.chrom_offset
        cumul_length = self.cumul_length
        abs_start_coords = self.abs_start_coords

        chrom_id1 = self.idmap.loc[chunk['chrom1'].values]
        chrom_id2 = self.idmap.loc[chunk['chrom2'].values]
        if binsize is None:
            abs_start1 = cumul_length[chrom_id1] + chunk['start1'].values
            abs_start2 = cumul_length[chrom_id2] + chunk['start2'].values
            chunk['bin1_id'] = np.searchsorted(abs_start_coords, abs_start1, side='right') - 1
            chunk['bin2_id'] = np.searchsorted(abs_start_coords, abs_start2, side='right') - 1
        else:
            rel_bin1 = np.floor(chunk['start1']/binsize).astype(int)
            rel_bin2 = np.floor(chunk['start2']/binsize).astype(int)
            chunk['bin1_id'] = chrom_offset[chrom_id1] + rel_bin1
            chunk['bin2_id'] = chrom_offset[chrom_id2] + rel_bin2
        
        gby = chunk.groupby(['bin1_id', 'bin2_id'], sort=False)
        agg = gby['count'].sum().reset_index()
        return {k: v.values for k,v in six.iteritems(agg)}

    def size(self):
        return self.cool.info['nnz']

    def __iter__(self):
        table = self.cool.pixels(join=True)
        chunksize = self.chunksize
        edges = np.arange(0, self.size() + chunksize, chunksize)
        for lo, hi in zip(edges[:-1], edges[1:]):
            chunk = table[lo:hi]
            print(lo, hi)
            yield self._aggregate(chunk)


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

