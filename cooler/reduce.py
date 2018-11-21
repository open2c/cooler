# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
from collections import OrderedDict
from bisect import bisect_right
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

from . import get_logger
from .io import parse_cooler_uri, ContactBinner, create
from .util import binnify, get_binsize, GenomeSegmentation
from .tools import lock


__all__ = ['merge', 'coarsen', 'zoomify']


logger = get_logger(__name__)


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
    Merge (i.e. sum) multiple cooler matrices having identical axes.

    """
    def __init__(self, coolers, maxbuf, **kwargs):
        self.coolers = list(coolers)
        self.maxbuf = maxbuf

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
                          .aggregate({'count': np.sum})
                          .reset_index())

            yield {k: v.values for k, v in six.iteritems(df)}

            starts = stops


class CoolerCoarsener(ContactBinner):
    """
    Aggregate contacts from an existing Cooler file.

    """
    def __init__(self, source_uri, bins, chunksize, batchsize, map=map):
        from .api import Cooler
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
        from .api import Cooler
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


QUAD_TILE_SIZE_PIXELS = 256


def get_quadtree_depth(chromsizes, binsize):
    """
    Depth of quad tree necessary to tesselate the concatenated genome with quad
    tiles such that linear dimension of the tiles is a preset multiple of the
    genomic resolution.

    """
    tile_size_bp = QUAD_TILE_SIZE_PIXELS * binsize
    min_tile_cover = math.ceil(sum(chromsizes) / tile_size_bp)
    return int(math.ceil(np.log2(min_tile_cover)))


def multires_aggregate(input_uri, outfile, nproc, chunksize, lock=None):
    """
    Quad-tree tiling for HiGlass

    """
    from .api import Cooler
    infile, ingroup = parse_cooler_uri(input_uri)

    clr = Cooler(infile, ingroup)
    n_zooms = get_quadtree_depth(clr.chromsizes, clr.binsize)
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

        coarsen(
            outfile + '::' + str(prevLevel),
            outfile + '::' + str(zoomLevel),
            factor,
            nproc,
            chunksize,
            lock
        )
        zoom_levels[zoomLevel] = binsize

    with h5py.File(outfile, 'r+') as fw:
        fw.attrs.update({'max-zoom': n_zooms})
        #grp = fw.require_group('.zooms')
        fw.attrs['max-zooms'] = n_zooms
        fw.attrs.update(zoom_levels)

    return n_zooms, zoom_levels


def merge(output_uri, input_uris, chunksize):
    from .api import Cooler
    logger.info("Merging:\n{}".format('\n'.join(input_uris)))
    clrs = [Cooler(path) for path in input_uris]
    chromsizes = clrs[0].chromsizes
    bins = clrs[0].bins()[['chrom', 'start', 'end']][:]
    assembly = clrs[0].info.get('genome-assembly', None)
    iterator = CoolerMerger(clrs, maxbuf=chunksize)

    is_symm = [clr.symmetric_storage_mode == u'upper' for clr in clrs]
    if all(is_symm):
        symmetric = True
    elif not any(is_symm):
        symmetric = False
    else:
        ValueError("Cannot merge symmetric and asymmetric coolers.")

    create(
        output_uri,
        bins,
        iterator,
        assembly=assembly,
        symmetric=symmetric
    )


def coarsen(input_uri, output_uri, factor, nproc, chunksize, lock=None):
    from .api import Cooler
    c = Cooler(input_uri)
    chromsizes = c.chromsizes
    new_binsize = c.binsize * factor
    new_bins = binnify(chromsizes, new_binsize)
    dtypes = dict(c.pixels()[0:0].dtypes.drop(['bin1_id', 'bin2_id']))

    try:
        # Note: fork before opening to prevent inconsistent global HDF5 state
        if nproc > 1:
            pool = mp.Pool(nproc)

        iterator = CoolerCoarsener(
            input_uri,
            new_bins,
            chunksize,
            batchsize=nproc,
            map=pool.map if nproc > 1 else map)

        create(
            output_uri,
            new_bins,
            iterator,
            dtypes=dtypes,
            symmetric=c.symmetric_storage_mode == u'upper',
            lock=lock,
            append=True)

    finally:
        if nproc > 1:
            pool.close()


def zoomify(input_uris, outfile, resolutions, nproc, chunksize, lock=None):
    from .api import Cooler
    uris = {}
    bases = set()
    for input_uri in input_uris:
        infile, ingroup = parse_cooler_uri(input_uri)
        base_binsize = Cooler(infile, ingroup).binsize
        uris[base_binsize] = (infile, ingroup)
        bases.add(base_binsize)

    resn, pred, mult = get_multiplier_sequence(resolutions, bases)
    n_zooms = len(resn)

    logger.info(
        "Copying base matrices and producing {} new zoom levels.".format(n_zooms)
    )

    # Copy base matrix
    for base_binsize in bases:
        logger.info("Bin size: " + str(base_binsize))
        infile, ingroup = uris[base_binsize]
        with h5py.File(infile, 'r') as src, \
            h5py.File(outfile, 'w') as dest:
            src.copy(ingroup, dest, '/resolutions/{}'.format(base_binsize))

    # Aggregate
    # Use lock to sync read/write ops on same file
    for i in range(n_zooms):
        if pred[i] == -1:
            continue
        prev_binsize = resn[pred[i]]
        binsize = prev_binsize * mult[i]
        logger.info(
            "Aggregating from {} to {}.".format(prev_binsize, binsize))
        coarsen(
            outfile + '::resolutions/{}'.format(prev_binsize),
            outfile + '::resolutions/{}'.format(binsize),
            mult[i],
            nproc,
            chunksize,
            lock
        )

    with h5py.File(outfile, 'r+') as fw:
        fw.attrs.update({
            'format': u'HDF5::MCOOL',
            'format-version': 2,
        })
