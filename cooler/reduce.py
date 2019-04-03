# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
from collections import OrderedDict, defaultdict
from bisect import bisect_right
from functools import partial
from six.moves import map
import multiprocess as mp
import os.path as op
import itertools
import shlex
import math
import sys
import six

import pandas as pd
import numpy as np
import h5py

from ._logging import get_logger
from .create import ContactBinner, create
from .util import parse_cooler_uri, binnify, get_binsize, GenomeSegmentation
from .tools import lock


__all__ = ['merge_coolers', 'coarsen_cooler', 'zoomify_cooler']


logger = get_logger(__name__)


HIGLASS_TILE_DIM = 256

ZOOMS_4DN = [
    1000,
    2000,
    5000,
    10000,
    25000,
    50000,
    100000,
    250000,
    500000,
    1000000,
    2500000,
    5000000,
    10000000
]


def merge_breakpoints(indexes, maxbuf):
    """
    Partition k offset arrays for performing a k-way external merge, such that
    no single merge pass loads more than ``maxbuf`` records  into memory, with
    one exception (see Notes).

    Parameters
    ----------
    indexes : sequence of 1D arrays of equal length
        These offset-array indexes map non-negative integers to their offset
        locations in a corresponding data table
    maxbuf : int
        Maximum cumulative number of records loaded into memory for a single
        merge pass

    Returns
    -------
    breakpoints : 1D array
        breakpoint locations to segment all the offset arrays
    cum_offset : 1D array
        cumulative number of records that will be processed at each breakpoint

    Notes
    -----
    The one exception to the post-condition is if any single increment of the
    indexes maps to more than ``maxbuf`` records, these will produce
    oversized chunks.

    """
    k = len(indexes)

    # the virtual cumulative index if no pixels were merged
    cumindex = np.vstack(indexes).sum(axis=0)
    cum_start = 0
    cum_nnz = cumindex[-1]
    n = len(cumindex)

    breakpoints = [0]
    cum_offsets = [0]
    lo = 0
    while True:
        # find the next mark
        hi = bisect_right(cumindex, min(cum_start + maxbuf, cum_nnz), lo=lo) - 1
        if hi == lo:
            # number of records to nearest mark exceeds `maxbuf`
            # check for oversized chunks afterwards
            hi += 1

        breakpoints.append(hi)
        cum_offsets.append(cumindex[hi])

        if cumindex[hi] == cum_nnz:
            break

        lo = hi
        cum_start = cumindex[hi]

    breakpoints = np.array(breakpoints)
    cum_offsets = np.array(cum_offsets)
    return breakpoints, cum_offsets


class CoolerMerger(ContactBinner):
    """
    Implementation of cooler merging.

    """
    def __init__(self, coolers, maxbuf, columns=None, agg=None):
        self.coolers = list(coolers)
        self.maxbuf = maxbuf
        self.columns = ['count'] if columns is None else columns
        self.agg = {col: 'sum' for col in self.columns}
        if agg is not None:
            self.agg.update(agg)

        # check compatibility between input coolers
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
        indexes = [c._load_dset('indexes/bin1_offset') for c in self.coolers]
        breakpoints, cum_offsets = merge_breakpoints(indexes, self.maxbuf)
        chunksizes = np.diff(cum_offsets)
        if chunksizes.max() > self.maxbuf:
            warnings.warn(
                'Some merge passes will use more than {} pixels'.format(
                    self.maxbuf))
        nnzs = [len(c.pixels()) for c in self.coolers]
        logger.info('nnzs: {}'.format(nnzs))

        starts = [0] * len(self.coolers)
        for bp in breakpoints[1:]:
            stops = [index[bp] for index in indexes]
            logger.info('current: {}'.format(stops))

            # extract, concat
            combined = pd.concat(
                [c.pixels()[start:stop]
                    for c, start, stop in zip(self.coolers, starts, stops)
                        if (stop - start) > 0],
                axis=0,
                ignore_index=True)

            # sort and aggregate
            df = (combined.groupby(['bin1_id', 'bin2_id'], sort=True)
                          .aggregate(self.agg)
                          .reset_index())

            yield {k: v.values for k, v in six.iteritems(df)}

            starts = stops


def merge_coolers(output_uri, input_uris, mergebuf, columns=None, dtypes=None,
                  agg=None, **kwargs):
    """
    Merge multiple coolers with identical axes.

    The merged cooler is stored at ``output_uri``.

    .. versionadded:: 0.8.0

    Parameters
    ----------
    output_uri : str
        Output cooler file path or URI.
    input_uris : list of str
        List of input file path or URIs of coolers to combine.
    mergebuf : int
        Maximum number of pixels processed at a time.
    columns : list of str, optional
        Specify which pixel value columns to include in the aggregation.
        Default is to use all available value columns.
    dtypes : dict, optional
        Specific dtypes to use for value columns. Default is to propagate
        the current dtypes of the value columns.
    agg : dict, optional
        Functions to use for aggregating each value column. Pass the same kind
        of dict accepted by ``pandas.DataFrame.groupby.agg``. Default is to
        apply 'sum' to every value column.
    kwargs
        Passed to ``cooler.create``.

    See also
    --------
    cooler.coarsen_cooler
    cooler.zoomify_cooler

    """
    #TODO: combine metadata from inputs
    from .api import Cooler
    logger.info("Merging:\n{}".format('\n'.join(input_uris)))

    clrs = [Cooler(path) for path in input_uris]

    is_symm = [clr.storage_mode == u'symmetric-upper' for clr in clrs]
    if all(is_symm):
        symmetric_upper = True
    elif not any(is_symm):
        symmetric_upper = False
    else:
        ValueError("Cannot merge symmetric and non-symmetric coolers.")

    if columns is None:
        columns = ['count']

    dtype_map = defaultdict(list)
    for clr in clrs:
        pixel_dtypes = clr.pixels().dtypes
        for col in columns:
            if col not in pixel_dtypes:
                raise ValueError(
                    "Pixel value column '{}' not found in "
                    "input '{}'.".format(col, clr.filename))
            else:
                dtype_map[col].append(pixel_dtypes[col])

    if dtypes is None:
        dtypes = {}
    for col in columns:
        if col not in dtypes:
            dtypes[col] = np.find_common_type(dtype_map[col], [])

    bins = clrs[0].bins()[['chrom', 'start', 'end']][:]
    assembly = clrs[0].info.get('genome-assembly', None)
    iterator = CoolerMerger(clrs, maxbuf=mergebuf, columns=columns, agg=agg)

    create(
        output_uri,
        bins,
        iterator,
        columns=columns,
        dtypes=dtypes,
        assembly=assembly,
        symmetric_upper=symmetric_upper,
        **kwargs
    )


def _optimal_prune_partition(edges, maxlen):
    """Given an integer interval partition ``edges``, find the coarsened
    partition with the longest subintervals such that no new subinterval
    created by removing edges exceeds ``maxlen``.
    """
    n = len(edges)
    if n < 2:
        raise ValueError("Partition must have 2 or more edges.")
    opt  = np.zeros(n, dtype=int)
    pred = np.zeros(n, dtype=int)

    opt[0] = 0
    for i in range(1, n):
        # default to immediate predecessor edge
        opt[i] = opt[i-1] + min(maxlen, edges[i] - edges[i-1])
        pred[i] = i - 1
        # try earlier predecessors until we exceed maxlen
        for k in range(i-2, -1, -1):
            length = edges[i] - edges[k]
            if length > maxlen:
                break
            s = opt[k] + length
            if s >= opt[i]:
                opt[i] = s
                pred[i] = k

    # backtrack to trace optimal path
    path = np.zeros(n, dtype=int)
    i = path[0] = n - 1
    j = 1
    while i > 0:
        i = path[j] = pred[i]
        j += 1
    path = path[:j][::-1]

    return edges[path]


def _greedy_prune_partition(edges, maxlen):
    """Given an integer interval partition ``edges`` from 0..nnz, prune the
    edges to make the new subintervals roughly ``maxlen`` in length.
    """
    edges = np.asarray(edges)
    assert len(edges) >= 2 and edges[0] == 0
    cumlen = np.r_[0, np.cumsum(np.diff(edges))]
    cuts = [maxlen * i  for i in range(0, int(np.ceil(cumlen[-1] / maxlen)))]
    cuts.append(cumlen[-1])
    idx = np.unique(np.searchsorted(cumlen, cuts))
    return edges[idx]


def get_quadtree_depth(chromsizes, base_binsize, bins_per_tile):
    """
    Number of zoom levels for a quad-tree tiling of a genomic heatmap.

    At the base resolution, we need N tiles, where N is the smallest power of
    2 such that the tiles fully cover the 1D data extent. From that starting
    point, determine the number of zoom levels required to "coarsen" the map
    up to 1 tile.

    """
    # length of a base-level tile in bp
    tile_length_bp = bins_per_tile * base_binsize

    # number of tiles required along 1 dimension at this base resolution
    total_bp = sum(chromsizes)
    n_tiles = math.ceil(total_bp / tile_length_bp)

    # number of aggregation levels required to reach a single tile
    n_zoom_levels = int(math.ceil(np.log2(n_tiles)))

    return n_zoom_levels


def get_multiplier_sequence(resolutions, bases=None):
    """
    From a set of target resolutions and one or more base resolutions
    deduce the most efficient sequence of integer multiple aggregations
    to satisfy all targets starting from the base resolution(s).

    Parameters
    ----------
    resolutions: sequence of int
        The target resolutions
    bases: sequence of int, optional
        The base resolutions for which data already exists.
        If not provided, the smallest resolution is assumed to be the base.

    Returns
    -------
    resn: 1D array
        Resolutions, sorted in ascending order.
    pred: 1D array
        Index of the predecessor resolution in `resn`. A value of -1 implies
        that the resolution is a base resolution.
    mult: 1D array
        Multiplier to go from predecessor to target resolution.

    """
    if bases is None:
        # assume the base resolution is the smallest one
        bases = {min(resolutions)}
    else:
        bases = set(bases)

    resn = np.array(sorted(bases.union(resolutions)))
    pred = -np.ones(len(resn), dtype=int)
    mult = -np.ones(len(resn), dtype=int)

    for i, target in list(enumerate(resn))[::-1]:
        p = i - 1
        while p >= 0:
            if target % resn[p] == 0:
                pred[i] = p
                mult[i] = target // resn[p]
                break
            else:
                p -= 1

    for i, p in enumerate(pred):
        if p == -1 and resn[i] not in bases:
            raise ValueError(
                "Resolution {} cannot be derived from "
                "the base resolutions: {}.".format(resn[i], bases))

    return resn, pred, mult


class CoolerCoarsener(ContactBinner):
    """
    Implementation of cooler coarsening.

    """
    def __init__(self, source_uri, factor, chunksize, columns, agg,
                 batchsize, map=map):
        from .api import Cooler
        self._map = map

        self.source_uri = source_uri
        self.batchsize = batchsize

        assert isinstance(factor, int) and factor > 1
        self.factor = factor
        self.chunksize = int(chunksize)

        self.index_columns = ['bin1_id', 'bin2_id']
        self.value_columns = list(columns)
        self.columns = self.index_columns + self.value_columns
        self.agg = {col: 'sum' for col in self.value_columns}
        if agg is not None:
            self.agg.update(agg)

        clr = Cooler(source_uri)
        chromsizes = clr.chromsizes

        # Info for the old bin segmentation
        self.old_binsize = clr.binsize
        self.old_chrom_offset = clr._load_dset('indexes/chrom_offset')
        self.old_bin1_offset  = clr._load_dset('indexes/bin1_offset')

        # Calculate the new bin segmentation
        if self.old_binsize is None:
            self.new_binsize = None
        else:
            self.new_binsize = self.old_binsize * factor
        old_bins = clr.bins()[['chrom', 'start', 'end']][:]
        self.new_bins = self.coarsen_bins(old_bins, chromsizes, factor)
        self.gs = GenomeSegmentation(chromsizes, self.new_bins)

        # Pre-compute the partition of bin1 offsets that groups the pixels into
        # coarsened bins along the i-axis. Then remove some of the internal
        # edges of this partition to make bigger groups of pixels. This way
        # we ensure that none of the original groups gets split.
        edges = []
        for chrom, i in six.iteritems(self.gs.idmap):
            # Respect chrom1 boundaries
            c0 = self.old_chrom_offset[i]
            c1 = self.old_chrom_offset[i + 1]
            edges.extend(self.old_bin1_offset[c0:c1:factor])
        edges.append(self.old_bin1_offset[-1])
        self.edges = _greedy_prune_partition(edges, self.chunksize)

    @staticmethod
    def coarsen_bins(old_bins, chromsizes, factor):
        def _each(group):
            out = group[['chrom', 'start']].copy().iloc[::factor]
            end = group['end'].iloc[factor-1::factor].values
            if len(end) < len(out):
                end = np.r_[end, chromsizes[group.name]]
            out['end'] = end
            return out
        return old_bins.groupby('chrom').apply(_each).reset_index(drop=True)

    def _aggregate(self, span):
        from .api import Cooler
        lo, hi = span

        clr = Cooler(self.source_uri)
        # convert_enum=False returns chroms as raw ints
        table = clr.pixels(join=True, convert_enum=False)[self.columns]
        chunk = table[lo:hi]
        logger.info('{} {}'.format(lo, hi))

        # use the "start" point as anchor for re-binning
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

        return (
            chunk.groupby(self.index_columns, sort=True)
                 .aggregate(self.agg)
                 .reset_index()
        )

    def aggregate(self, span):
        try:
            chunk = self._aggregate(span)
        except MemoryError as e:
            raise RuntimeError(str(e))
        return chunk

    def __iter__(self):
        # Distribute batches of `batchsize` pixel spans at once.
        batchsize = self.batchsize
        spans = list(zip(self.edges[:-1], self.edges[1:]))
        for i in range(0, len(spans), batchsize):
            try:
                if batchsize > 1:
                    lock.acquire()
                results = self._map(self.aggregate, spans[i:i+batchsize])
            finally:
                if batchsize > 1:
                    lock.release()
            for df in results:
                yield {k: v.values for k, v in six.iteritems(df)}


def coarsen_cooler(base_uri, output_uri, factor, chunksize, nproc=1,
                   columns=None, dtypes=None, agg=None, **kwargs):
    """
    Coarsen a cooler to a lower resolution by an integer factor *k*.

    This is done by pooling *k*-by-*k* neighborhoods of pixels and aggregating.
    Each chromosomal block is coarsened individually. Result is a coarsened
    cooler stored at ``output_uri``.

    .. versionadded:: 0.8.0

    Parameters
    ----------
    base_uri : str
        Input cooler file path or URI.
    output_uri : str
        Input cooler file path or URI.
    factor : int
        Coarsening factor.
    chunksize : int
        Number of pixels processed at a time per worker.
    nproc : int, optional
        Number of workers for batch processing of pixels. Default is 1,
        i.e. no process pool.
    columns : list of str, optional
        Specify which pixel value columns to include in the aggregation.
        Default is to use all available value columns.
    dtypes : dict, optional
        Specific dtypes to use for value columns. Default is to propagate
        the current dtypes of the value columns.
    agg : dict, optional
        Functions to use for aggregating each value column. Pass the same kind
        of dict accepted by ``pandas.DataFrame.groupby.agg``. Default is to
        apply 'sum' to every value column.
    kwargs
        Passed to ``cooler.create``.

    See also
    --------
    cooler.zoomify_cooler
    cooler.merge_coolers

    """
    #TODO: decide whether to default to 'count' or whatever is there besides bin1_id, bin2_id
    # dtypes = dict(clr.pixels().dtypes.drop(['bin1_id', 'bin2_id']))

    from .api import Cooler
    clr = Cooler(base_uri)

    factor = int(factor)

    if columns is None:
        columns = ['count']

    if dtypes is None:
        dtypes = {}

    input_dtypes = clr.pixels().dtypes
    for col in columns:
        if col not in input_dtypes:
            raise ValueError(
                "Pixel value column '{}' not found in "
                "input '{}'.".format(col, clr.filename))
        else:
            dtypes[col] = input_dtypes[col]

    try:
        # Note: fork before opening to prevent inconsistent global HDF5 state
        if nproc > 1:
            pool = mp.Pool(nproc)
            kwargs.setdefault('lock', lock)

        iterator = CoolerCoarsener(
            base_uri,
            factor,
            chunksize,
            columns=columns,
            agg=agg,
            batchsize=nproc,
            map=pool.map if nproc > 1 else map)

        new_bins = iterator.new_bins

        kwargs.setdefault('append', True)

        create(
            output_uri,
            new_bins,
            iterator,
            dtypes=dtypes,
            symmetric_upper=clr.storage_mode == u'symmetric-upper',
            **kwargs)

    finally:
        if nproc > 1:
            pool.close()


def zoomify_cooler(base_uris, outfile, resolutions, chunksize, nproc=1,
                   columns=None, dtypes=None, agg=None, **kwargs):
    """
    Generate multiple cooler resolutions by recursive coarsening.

    Result is a "zoomified" or "multires" cool file stored at ``outfile``
    using the MCOOL v2 layout, where coolers are stored under a hierarchy of
    the form ``resolutions/<r>`` for each resolution ``r``.

    .. versionadded:: 0.8.0

    Parameters
    ----------
    base_uris : str or sequence of str
        One or more cooler URIs to use as "base resolutions" for aggregation.
    outfile : str
        Output multires cooler (mcool) file path.
    resolutions : list of int
        A list of target resolutions to generate.
    chunksize : int
        Number of pixels processed at a time per worker.
    nproc : int, optional
        Number of workers for batch processing of pixels. Default is 1,
        i.e. no process pool.
    columns : list of str, optional
        Specify which pixel value columns to include in the aggregation.
        Default is to use only the column named 'count' if it exists.
    dtypes : dict, optional
        Specific dtypes to use for value columns. Default is to propagate
        the current dtypes of the value columns.
    agg : dict, optional
        Functions to use for aggregating each value column. Pass the same kind
        of dict accepted by ``pandas.DataFrame.groupby.agg``. Default is to
        apply 'sum' to every value column.
    kwargs
        Passed to ``cooler.create``.

    See also
    --------
    cooler.coarsen_cooler
    cooler.merge_coolers

    """
    # TODO: provide presets? {'pow2', '4dn'}
    from .api import Cooler

    if isinstance(base_uris, six.string_types):
        base_uris = [base_uris]

    parsed_uris = {}
    n_bins_longest_chrom = {}
    base_resolutions = set()
    for input_uri in base_uris:
        infile, ingroup = parse_cooler_uri(input_uri)
        clr = Cooler(infile, ingroup)
        base_binsize = clr.binsize
        parsed_uris[base_binsize] = (infile, ingroup)
        n_bins_longest_chrom[base_binsize] = clr.bins()[:].groupby('chrom').size().max()
        base_resolutions.add(base_binsize)

    # Determine the sequence of reductions.
    resn, pred, mult = get_multiplier_sequence(resolutions, base_resolutions)

    n_zooms = len(resn)

    logger.info(
        "Copying base matrices and producing {} new zoom levels.".format(n_zooms)
    )

    if columns is None:
        columns = ['count']

    # Copy base matrix
    for base_binsize in base_resolutions:
        logger.info("Bin size: " + str(base_binsize))
        infile, ingroup = parsed_uris[base_binsize]
        with h5py.File(infile, 'r') as src, \
             h5py.File(outfile, 'w') as dest:
            prefix = '/resolutions/{}'.format(base_binsize)

            src.copy(ingroup + '/chroms',
                     dest, prefix + '/chroms')
            src.copy(ingroup + '/bins',
                     dest, prefix + '/bins')
            for col in ['bin1_id', 'bin2_id'] + list(columns):
                src.copy(ingroup + '/pixels/{}'.format(col),
                         dest, prefix + '/pixels/{}'.format(col))
            src.copy(ingroup + '/indexes',
                     dest, prefix + '/indexes')
            dest[prefix].attrs.update(src[ingroup].attrs)

    # Aggregate
    # Use lock to sync read/write ops on same file
    for i in range(n_zooms):
        if pred[i] == -1:
            continue
        prev_binsize = resn[pred[i]]
        binsize = prev_binsize * mult[i]
        logger.info(
            "Aggregating from {} to {}.".format(prev_binsize, binsize))
        coarsen_cooler(
            outfile + '::resolutions/{}'.format(prev_binsize),
            outfile + '::resolutions/{}'.format(binsize),
            mult[i],
            chunksize,
            nproc=nproc,
            columns=columns,
            dtypes=dtypes,
            agg=agg,
            **kwargs
        )

    with h5py.File(outfile, 'r+') as fw:
        fw.attrs.update({
            'format': u'HDF5::MCOOL',
            'format-version': 2,
        })


def legacy_zoomify(input_uri, outfile, nproc, chunksize, lock=None):
    """
    Quad-tree tiling using legacy MCOOL layout (::0, ::1, ::2, etc.).

    """
    from .api import Cooler
    infile, ingroup = parse_cooler_uri(input_uri)

    clr = Cooler(infile, ingroup)
    n_zooms = get_quadtree_depth(clr.chromsizes, clr.binsize, HIGLASS_TILE_DIM)
    factor = 2

    logger.info("total_length (bp): {}".format(np.sum(clr.chromsizes)))
    logger.info("binsize: {}".format(clr.binsize))
    logger.info("n_zooms: {}".format(n_zooms))
    logger.info("quad tile cover: {}".format(2**n_zooms))
    logger.info(
        "Copying base matrix to level " +
        "{0} and producing {0} new zoom levels ".format(n_zooms) +
        "counting down to 0..."
    )

    zoom_levels = OrderedDict()
    zoomLevel = str(n_zooms)
    binsize = clr.binsize
    logger.info(
        "Zoom level: "
        + str(zoomLevel)
        + " bin size: "
        + str(binsize))

    # Copy base matrix
    with h5py.File(infile, 'r') as src, \
         h5py.File(outfile, 'w') as dest:

        src.copy(ingroup, dest, str(zoomLevel))
        zoom_levels[zoomLevel] = binsize

    # Aggregate
    # Use lock to sync read/write ops on same file
    for i in range(n_zooms - 1, -1, -1):
        prev_binsize = binsize
        binsize *= factor
        prevLevel = str(i+1)
        zoomLevel = str(i)
        logger.info(
            "Aggregating at zoom level: "
            + str(zoomLevel)
            + " bin size: "
            + str(binsize))

        coarsen_cooler(
            outfile + '::' + str(prevLevel),
            outfile + '::' + str(zoomLevel),
            factor,
            chunksize=chunksize,
            nproc=nproc,
            lock=lock
        )
        zoom_levels[zoomLevel] = binsize

    with h5py.File(outfile, 'r+') as fw:
        fw.attrs.update({'max-zoom': n_zooms})
        #grp = fw.require_group('.zooms')
        fw.attrs['max-zooms'] = n_zooms
        fw.attrs.update(zoom_levels)

    return n_zooms, zoom_levels
