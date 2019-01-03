# -*- coding: utf-8 -*-
"""
Contact Binners
~~~~~~~~~~~~~~~

Binners are iterators that convert input data of various flavors into a
properly sorted, chunked stream of binned contacts.

"""
from __future__ import division, print_function
from collections import OrderedDict, Counter
from bisect import bisect_left, bisect_right
from multiprocess import Pool
from functools import partial
import itertools
import warnings
import sys
import six

from pandas.api.types import is_integer_dtype
import numpy as np
import pandas as pd
import h5py

from .._logging import get_logger
from ..util import (
    parse_region, rlencode, get_binsize, get_chromsizes, GenomeSegmentation,
    balanced_partition
)
from ..tools import lock, partition


logger = get_logger('cooler.create')


class BadInputError(ValueError):
    pass


SANITIZE_PRESETS = {
    'bg2': dict(decode_chroms=True, is_one_based=False,  tril_action='reflect',
                chrom_field='chrom', anchor_field='start',
                sided_fields=('chrom', 'start', 'end'), suffixes=('1', '2'),
                sort=True, validate=True),

    'pairs': dict(decode_chroms=True, is_one_based=False, tril_action='reflect',
                  chrom_field='chrom', anchor_field='pos',
                  sided_fields=('chrom', 'pos'), suffixes=('1', '2'),
                  sort=False, validate=True)
}



def _sanitize_records(chunk, gs, decode_chroms, is_one_based, tril_action,
                      chrom_field, anchor_field, sided_fields, suffixes,
                      sort, validate):
    # Get integer contig IDs
    if decode_chroms:
        # Unspecified chroms get assigned category = NaN and integer code = -1
        chrom1_ids = np.array(pd.Categorical(
            chunk['chrom1'], gs.contigs, ordered=True).codes)
        chrom2_ids = np.array(pd.Categorical(
            chunk['chrom2'], gs.contigs, ordered=True).codes)
    else:
        chrom1_ids = chunk['chrom1'].values
        chrom2_ids = chunk['chrom2'].values
        if validate:
            for col, dt in [('chrom1', chrom1_ids.dtype),
                            ('chrom2', chrom2_ids.dtype)]:
                if not is_integer_dtype(dt):
                    raise BadInputError(
                        "`{}` column is non-integer. ".format(col) +
                        "If string, use `decode_chroms=True` to convert to enum")

    # Drop records from non-requested chromosomes
    to_drop = (chrom1_ids < 0) | (chrom2_ids < 0)
    if np.any(to_drop):
        mask = ~to_drop
        chrom1_ids = chrom1_ids[mask]
        chrom2_ids = chrom2_ids[mask]
        chunk = chunk[mask].copy()

    # Handle empty case
    if not len(chunk):
        chunk['bin1_id'] = []
        chunk['bin2_id'] = []
        return chunk

    # Find positional anchor columns, convert to zero-based if needed
    anchor1 = np.array(chunk[anchor_field + suffixes[0]])
    anchor2 = np.array(chunk[anchor_field + suffixes[1]])
    if is_one_based:
        anchor1 -= 1
        anchor2 -= 1

    # Check types and bounds
    if validate:
        for dt in [anchor1.dtype, anchor2.dtype]:
            if not is_integer_dtype(dt):
                raise BadInputError("Found a non-integer anchor column")

        is_neg = (anchor1 < 0) | (anchor2 < 0)
        if np.any(is_neg):
            err = chunk[is_neg]
            raise BadInputError(
                "Found an anchor position with negative value:\n{}".format(
                err.head().to_csv(sep='\t')))

        chromsizes1 = gs.chromsizes[chrom1_ids].values
        chromsizes2 = gs.chromsizes[chrom2_ids].values
        is_excess = (anchor1 > chromsizes1) | (anchor2 > chromsizes2)
        if np.any(is_excess):
            err = chunk[is_excess]
            raise BadInputError(
                "Found an anchor position exceeding chromosome length:\n{}".format(
                err.head().to_csv(sep='\t')))

    # Handle lower triangle records
    if tril_action is not None:
        is_tril = (
            (chrom1_ids > chrom2_ids) |
            ((chrom1_ids == chrom2_ids) & (anchor1 > anchor2))
        )
        if np.any(is_tril):
            if tril_action == 'reflect':
                chrom1_ids[is_tril], chrom2_ids[is_tril] = \
                    chrom2_ids[is_tril], chrom1_ids[is_tril]
                anchor1[is_tril], anchor2[is_tril] = \
                    anchor2[is_tril], anchor1[is_tril]
                for field in sided_fields:
                    chunk.loc[is_tril, field + suffixes[0]], \
                    chunk.loc[is_tril, field + suffixes[1]] = \
                        chunk.loc[is_tril, field + suffixes[1]], \
                        chunk.loc[is_tril, field + suffixes[0]]
            elif tril_action == 'drop':
                mask = ~is_tril
                chrom1_ids = chrom1_ids[mask]
                chrom2_ids = chrom2_ids[mask]
                anchor1 = anchor1[mask]
                anchor2 = anchor2[mask]
                chunk = chunk[mask].copy()
            elif tril_action == 'raise':
                err = chunk[is_tril]
                raise BadInputError("Found lower triangle pairs:\n{}".format(
                    err.head().to_csv(sep='\t')))
            else:
                raise ValueError("Unknown tril_action value: '{}'".format(
                    tril_action))

    # Assign bin IDs from bin table
    chrom_binoffset = gs.chrom_binoffset
    binsize = gs.binsize
    if binsize is None:
        chrom_abspos = gs.chrom_abspos
        start_abspos = gs.start_abspos
        bin1_ids = []
        bin2_ids = []
        for cid1, pos1, cid2, pos2 in zip(chrom1_ids, anchor1,
                                          chrom2_ids, anchor2):
            lo = chrom_binoffset[cid1]
            hi = chrom_binoffset[cid1 + 1]
            bin1_ids.append(lo + np.searchsorted(
                                start_abspos[lo:hi],
                                chrom_abspos[cid1] + pos1,
                                side='right') - 1)
            lo = chrom_binoffset[cid2]
            hi = chrom_binoffset[cid2 + 1]
            bin2_ids.append(lo + np.searchsorted(
                                start_abspos[lo:hi],
                                chrom_abspos[cid2] + pos2,
                                side='right') - 1)
        chunk['bin1_id'] = bin1_ids
        chunk['bin2_id'] = bin2_ids
    else:
        chunk['bin1_id'] = chrom_binoffset[chrom1_ids] + anchor1 // binsize
        chunk['bin2_id'] = chrom_binoffset[chrom2_ids] + anchor2 // binsize

    # Sort by bin IDs
    if sort:
        chunk = chunk.sort_values(['bin1_id', 'bin2_id'])

    # TODO: check for duplicate records and warn

    return chunk


def sanitize_records(bins, schema=None, **kwargs):
    """
    Builds a funtion to sanitize and assign bin IDs to a data frame of
    paired genomic positions based on a provided genomic bin segmentation.

    Parameters
    ----------
    bins : DataFrame
        Bin table to compare records against.
    schema : str, optional
        Use pre-defined parameters for a particular format. Any options can be
        overriden via kwargs. If not provided, values for all the options below
        must be given.

    decode_chroms : bool
        Convert string chromosome names to integer IDs based on the order given
        in the bin table. Set to False if the chromosomes are already given as
        an enumeration, starting at 0. Records with either chrom ID < 0 are
        dropped.
    is_one_based : bool
        Whether the input anchor coordinates are one-based, rather than
        zero-based. They will be converted to zero-based.
    tril_action : 'reflect', 'drop', 'raise' or None
        How to handle lower triangle ("tril") records.
        If set to 'reflect', tril records will be flipped or "reflected"
        to their mirror image: "sided" column pairs will have their values
        swapped.
        If set to 'drop', tril records will be discarded. This is useful if
        your input data is symmetric, i.e. contains mirror duplicates of every
        record.
        If set to 'raise', an exception will be raised if any tril record is
        encountered.
    chrom_field : str
        Base name of the two chromosome/scaffold/contig columns.
    anchor_field : str
        Base name of the positional anchor columns.
    sided_fields : sequence of str
        Base names of column pairs to swap values between when
        mirror-reflecting records.
    suffixes : pair of str
        Suffixes used to identify pairs of sided columns. e.g.: ('1', '2'),
        ('_x', '_y'), etc.
    sort : bool
        Whether to sort the output dataframe by bin_id and bin2_id.
    validate : bool
        Whether to do type- and bounds-checking on the anchor position
        columns. Raises BadInputError.

    Returns
    -------
    callable :
        Function of one argument that takes a raw dataframe and returns a sanitized dataframe with bin IDs assigned.

    """
    if schema is not None:
        try:
            options = SANITIZE_PRESETS[schema].copy()
        except KeyError:
            raise ValueError("Unknown schema: '{}'".format(schema))
    else:
        options = {}
    options.update(**kwargs)
    chromsizes = get_chromsizes(bins)
    options['gs'] = GenomeSegmentation(chromsizes, bins)
    return partial(_sanitize_records, **options)


def _sanitize_pixels(chunk, gs, is_one_based=False, tril_action='reflect',
                    bin1_field='bin1_id', bin2_field='bin2_id', sided_fields=(),
                    suffixes=('1', '2'), validate=True, sort=True):
    if is_one_based:
        chunk['bin1_id'] -= 1
        chunk['bin2_id'] -= 1

    if tril_action is not None:
        is_tril = chunk['bin1_id'] > chunk['bin2_id']
        if np.any(is_tril):
            if tril_action == 'reflect':
                chunk.loc[is_tril, 'bin1_id'], \
                    chunk.loc[is_tril, 'bin2_id'] = chunk.loc[is_tril, 'bin2_id'], \
                                                       chunk.loc[is_tril, 'bin1_id']
                for field in sided_fields:
                    chunk.loc[is_tril, field + suffixes[0]], \
                    chunk.loc[is_tril, field + suffixes[1]] = \
                        chunk.loc[is_tril, field + suffixes[1]], \
                        chunk.loc[is_tril, field + suffixes[0]]
            elif tril_action == 'drop':
                chunk = chunk[~is_tril]
            elif tril_action == 'raise':
                raise BadInputError("Found bin1_id greater than bin2_id")
            else:
                raise ValueError("Unknown tril_action value: '{}'".format(
                    tril_action))

    return chunk.sort_values(['bin1_id', 'bin2_id']) if sort else chunk


def _validate_pixels(chunk, n_bins, boundscheck, triucheck, dupcheck, ensure_sorted):
    if boundscheck:
        is_neg = (chunk['bin1_id'] < 0) | (chunk['bin2_id'] < 0)
        if np.any(is_neg):
            raise BadInputError("Found bin ID < 0")
        is_excess = (chunk['bin1_id'] >= n_bins) | (chunk['bin2_id'] >= n_bins)
        if np.any(is_excess):
            raise BadInputError(
                "Found a bin ID that exceeds the declared number of bins. "
                "Check whether your bin table is correct.")

    if triucheck:
        is_tril = chunk['bin1_id'] > chunk['bin2_id']
        if np.any(is_tril):
            raise BadInputError("Found bin1_id greater than bin2_id")

    if not isinstance(chunk, pd.DataFrame):
        chunk = pd.DataFrame(chunk)

    if dupcheck:
        is_dup = chunk.duplicated(['bin1_id', 'bin2_id'])
        if is_dup.any():
            err = chunk[is_dup]
            raise BadInputError("Found duplicate pixels:\n{}".format(err.head().to_csv(sep='\t')))

    if ensure_sorted:
        chunk = chunk.sort_values(['bin1_id', 'bin2_id'])

    return chunk


def validate_pixels(n_bins, boundscheck, triucheck, dupcheck, ensure_sorted):
    return partial(
        _validate_pixels,
        n_bins=n_bins,
        boundscheck=boundscheck,
        triucheck=triucheck,
        dupcheck=dupcheck,
        ensure_sorted=ensure_sorted)


def sanitize_pixels(bins, **kwargs):
    """
    Builds a function to sanitize an already-binned genomic data with
    genomic bin assignments.

    Parameters
    ----------
    bins : DataFrame
        Bin table to compare pixel records against.
    is_one_based : bool, optional
        Whether the input bin IDs are one-based, rather than zero-based.
        They will be converted to zero-based.
    tril_action : 'reflect', 'drop', 'raise' or None
        How to handle lower triangle ("tril") pixels.
        If set to 'reflect' [default], tril pixels will be flipped or
        "reflected" to their mirror image: "sided" column pairs will have their
        values swapped.
        If set to 'drop', tril pixels will be discarded. This is useful if
        your input data is duplexed, i.e. contains mirror duplicates of every
        record.
        If set to 'raise', an exception will be raised if any tril record is
        encountered.
    bin1_field : str
        Name of the column representing ith (row) axis of the matrix.
        Default is 'bin1_id'.
    bin2_field : str
        Name of the column representing jth (col) axis of the matrix.
        Default is 'bin2_id'.
    sided_fields : sequence of str
        Base names of column pairs to swap values between when mirror-reflecting
        pixels.
    suffixes : pair of str
        Suffixes used to identify pairs of sided columns. e.g.: ('1', '2'),
        ('_x', '_y'), etc.
    sort : bool
        Whether to sort the output dataframe by bin_id and bin2_id.
    validate : bool
        Whether to do type- and bounds-checking on the bin IDs.
        Raises BadInputError.

    Returns
    -------
    callable :
        Function of one argument that takes a raw dataframe and returns a
        sanitized dataframe.

    """
    chromsizes = get_chromsizes(bins)
    kwargs['gs'] = GenomeSegmentation(chromsizes, bins)
    return partial(_sanitize_pixels, **kwargs)


def _aggregate_records(chunk, sort, agg, rename):
    return (chunk.groupby(['bin1_id', 'bin2_id'], sort=sort)
                 .aggregate(agg)
                 .rename(columns=rename)
                 .reset_index())


def aggregate_records(sort=True, count=True, agg=None, rename=None):
    """
    Generates a function that aggregates bin-assigned records by pixel.

    Parameters
    ----------
    sort : bool, optional
        Sort group keys. Get better performance by turning this off.
        Note that this does not influence the order of observations within each
        group.
    count : bool, optional
        Output the number of records per pixel. Default is True.
    agg : dict, optional
        Dict of column names -> functions or names.
    rename : dict, optional
        Dict to rename columns after aggregating.

    Returns
    -------
    Function that takes a dataframe of records with bin IDs assigned, groups
    them by pixel, counts them, and optionally aggregates other value columns.

    Notes
    -----
    The GroupBy 'count' method ignores NaNs within groups, as opposed to 'size'.

    """
    if agg is None:
        agg = {}

    if rename is None:
        rename = {}

    # We use one of the grouper columns to count the number of pairs per pixel.
    # We always do count, even if 'count' isn't requested as output.
    if count and 'count' not in agg:
        agg['bin1_id'] = 'size'
        rename['bin1_id'] = 'count'

    def _aggregate_records(chunk):
        return (chunk.groupby(['bin1_id', 'bin2_id'], sort=sort)
                     .aggregate(agg)
                     .rename(columns=rename)
                     .reset_index())
    return _aggregate_records

    # return partial(_aggregate_records, sort=sort, agg=agg, rename=rename)


class ContactBinner(object):
    """
    Base class for iterable contact binners.

    """
    def __getstate__(self):
        d = self.__dict__.copy()
        d.pop('_map', None)
        return d

    def __iter__(self):
        """ Iterator over chunks of binned contacts (i.e., nonzero pixels)

        Chunks are expected to have the following format:

        * dict of 1D arrays
        * keys `bin1_id`, `bin2_id`, `count`
        * arrays lexically sorted by `bin_id` then `bin2_id`

        """
        raise NotImplementedError


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
        return pd.DataFrame(data)

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
    def __init__(self, filepath, chromsizes, bins, map=map, n_chunks=1, is_one_based=False, **kwargs):
        try:
            import pysam
        except ImportError:
            raise ImportError("pysam is required to read tabix files")

        import dill
        import pickle
        dill.settings['protocol'] = pickle.HIGHEST_PROTOCOL

        self._map = map
        self.n_chunks = n_chunks
        self.is_one_based = bool(is_one_based)
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
            if not len(self.file_contigs):
                raise RuntimeError("No reference sequences found.")

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
        decr = int(self.is_one_based)
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
                    pos2 = int(line[P2]) - decr

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
                    pd.DataFrame({
                        'bin1_id': bin1_id,
                        'bin2_id': list(accumulator.keys()),
                        'count':   list(accumulator.values())},
                        columns=['bin1_id', 'bin2_id', 'count'])
                          .sort_values('bin2_id')
                )

                accumulator.clear()

        logger.info('Finished {}:{}-{}|*'.format(chrom1, start, end))

        return pd.concat(rows, axis=0) if len(rows) else None

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
    def __init__(self, filepath, chromsizes, bins, map=map, n_chunks=1, is_one_based=False, **kwargs):
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
        self.is_one_based = bool(is_one_based)
        f = pypairix.open(filepath, 'r')
        self.C1 = f.get_chr1_col()
        self.C2 = f.get_chr2_col()
        self.P1 = f.get_startpos1_col()
        self.P2 = f.get_startpos2_col()
        self.file_contigs = set(
            itertools.chain.from_iterable(
                [b.split('|') for b in f.get_blocknames()]))

        if not len(self.file_contigs):
            raise RuntimeError("No reference sequences found.")
        for c1, c2 in itertools.combinations(self.file_contigs, 2):
            if f.exists2(c1, c2) and f.exists2(c2, c1):
                raise RuntimeError(
                    "Pairs are not triangular: found blocks " +
                    "'{0}|{1}'' and '{1}|{0}'".format(c1, c2))

        # dumb heuristic to prevent excessively large chunks on one worker
        if hasattr(f, 'get_linecount'):
            n_lines = f.get_linecount()
            if n_lines < 0:
                # correct int32 overflow bug
                MAXINT32 = 2147483647
                n_lines = MAXINT32 + MAXINT32 + n_lines
            max_chunk = int(100e6)
            n_chunks = n_lines // 2 // max_chunk
            old_n = self.n_chunks
            self.n_chunks = max(self.n_chunks, n_chunks)
            if self.n_chunks > old_n:
                logger.info(
                    "Pairs file has {} lines. Increasing max-split to {}.".format(
                    n_lines, self.n_chunks))

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
        decr = int(self.is_one_based)
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

                    pos2 = int(line[pos2_col]) - decr

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
                pd.DataFrame({
                    'bin1_id': bin1_id,
                    'bin2_id': list(accumulator.keys()),
                    'count':   list(accumulator.values())},
                    columns=['bin1_id', 'bin2_id', 'count'])
                      .sort_values('bin2_id')
            )

            accumulator.clear()

        logger.info('Finished {}:{}-{}|*'.format(chrom1, start, end))

        return pd.concat(rows, axis=0) if len(rows) else None

    def __iter__(self):
        granges = balanced_partition(self.gs, self.n_chunks, self.file_contigs)
        for df in self._map(self.aggregate, granges):
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
                 field_dtypes=None, one_based=False):
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
        self.one_based_bin_ids = one_based
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

        iterator = pd.read_csv(
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
            if self.one_based_bin_ids:
                # convert to zero-based
                if np.any(chunk['bin1_id'] <= 0) or np.any(chunk['bin2_id'] <= 0):
                    raise ValueError(
                        "Found bin ID <= 0. Are you sure bin IDs are one-based?")
                chunk['bin1_id'] -= 1
                chunk['bin2_id'] -= 1
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
