import math
import warnings
from bisect import bisect_right
from collections import OrderedDict, defaultdict

import h5py
import multiprocess as mp
import numpy as np
import pandas as pd

from ._logging import get_logger
from ._version import __format_version_mcool__
from .create import ContactBinner, create
from .parallel import lock
from .util import GenomeSegmentation, parse_cooler_uri

__all__ = ["merge_coolers", "coarsen_cooler", "zoomify_cooler"]


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
    10000000,
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
    # k = len(indexes)

    # the virtual cumulative index if no pixels were merged
    cumindex = np.vstack(indexes).sum(axis=0)
    cum_start = 0
    cum_nnz = cumindex[-1]
    # n = len(cumindex)

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
        self.columns = ["count"] if columns is None else columns
        self.agg = {col: "sum" for col in self.columns}
        if agg is not None:
            self.agg.update(agg)

        # check compatibility between input coolers
        binsize = coolers[0].binsize
        if binsize is not None:
            if len({c.binsize for c in coolers}) > 1:
                raise ValueError("Coolers must have the same resolution")
            chromsizes = coolers[0].chromsizes
            for i in range(1, len(coolers)):
                if not np.all(coolers[i].chromsizes == chromsizes):
                    raise ValueError("Coolers must have the same chromosomes")
        else:
            bins = coolers[0].bins()[["chrom", "start", "end"]][:]
            for i in range(1, len(coolers)):
                bins2 = coolers[i].bins()[["chrom", "start", "end"]][:]
                if (len(bins2) != len(bins)) or not np.all(bins2 == bins):
                    raise ValueError("Coolers must have same bin structure")

    def __iter__(self):
        indexes = [c._load_dset("indexes/bin1_offset") for c in self.coolers]
        breakpoints, cum_offsets = merge_breakpoints(indexes, self.maxbuf)
        chunksizes = np.diff(cum_offsets)
        if chunksizes.max() > self.maxbuf:
            warnings.warn(f"Some merge passes will use more than {self.maxbuf} pixels")
        nnzs = [len(c.pixels()) for c in self.coolers]
        logger.info(f"nnzs: {nnzs}")

        starts = [0] * len(self.coolers)
        for bp in breakpoints[1:]:
            stops = [index[bp] for index in indexes]
            logger.info(f"current: {stops}")

            # extract, concat
            combined = pd.concat(
                [
                    c.pixels()[start:stop]
                    for c, start, stop in zip(self.coolers, starts, stops)
                    if (stop - start) > 0
                ],
                axis=0,
                ignore_index=True,
            )

            # sort and aggregate
            df = (
                combined.groupby(["bin1_id", "bin2_id"], sort=True)
                .aggregate(self.agg)
                .reset_index()
            )

            yield {k: v.values for k, v in df.items()}

            starts = stops


def merge_coolers(
    output_uri, input_uris, mergebuf, columns=None, dtypes=None, agg=None, **kwargs
):
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

    Notes
    -----
    The default output file mode is 'w'. If appending output to an existing
    file, pass `mode='a'`.

    See also
    --------
    cooler.coarsen_cooler
    cooler.zoomify_cooler

    """
    # TODO: combine metadata from inputs
    from .api import Cooler

    logger.info("Merging:\n{}".format("\n".join(input_uris)))

    clrs = [Cooler(path) for path in input_uris]

    is_symm = [clr.storage_mode == "symmetric-upper" for clr in clrs]
    if all(is_symm):
        symmetric_upper = True
    elif not any(is_symm):
        symmetric_upper = False
    else:
        raise ValueError("Cannot merge symmetric and non-symmetric coolers.")

    if columns is None:
        columns = ["count"]

    dtype_map = defaultdict(list)
    for clr in clrs:
        pixel_dtypes = clr.pixels().dtypes
        for col in columns:
            if col not in pixel_dtypes:
                raise ValueError(
                    f"Pixel value column '{col}' not found in "
                    f"input '{clr.filename}'."
                )
            else:
                dtype_map[col].append(pixel_dtypes[col])

    if dtypes is None:
        dtypes = {}
    for col in columns:
        if col not in dtypes:
            dtypes[col] = np.find_common_type(dtype_map[col], [])

    bins = clrs[0].bins()[["chrom", "start", "end"]][:]
    assembly = clrs[0].info.get("genome-assembly", None)
    iterator = CoolerMerger(clrs, maxbuf=mergebuf, columns=columns, agg=agg)

    create(
        output_uri,
        bins,
        iterator,
        columns=columns,
        dtypes=dtypes,
        assembly=assembly,
        symmetric_upper=symmetric_upper,
        **kwargs,
    )


def _optimal_prune_partition(edges, maxlen):  # pragma: no cover
    """Given an integer interval partition ``edges``, find the coarsened
    partition with the longest subintervals such that no new subinterval
    created by removing edges exceeds ``maxlen``.
    """
    n = len(edges)
    if n < 2:
        raise ValueError("Partition must have 2 or more edges.")
    opt = np.zeros(n, dtype=int)
    pred = np.zeros(n, dtype=int)

    opt[0] = 0
    for i in range(1, n):
        # default to immediate predecessor edge
        opt[i] = opt[i - 1] + min(maxlen, edges[i] - edges[i - 1])
        pred[i] = i - 1
        # try earlier predecessors until we exceed maxlen
        for k in range(i - 2, -1, -1):
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
    cuts = [maxlen * i for i in range(0, int(np.ceil(cumlen[-1] / maxlen)))]
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


def geomprog(start, mul):
    """
    Generate a geometric progression of integers.

    Beginning with integer ``start``, yield an unbounded geometric progression
    with integer ratio ``mul``.

    """
    start, mul = int(start), int(mul)
    yield start
    while True:
        start *= mul
        yield start


def niceprog(start):
    """
    Generate a nice progression of integers.

    Beginning with integer ``start``, yield a sequence of "nicely" spaced
    integers: an unbounded geometric progression with ratio 10, interspersed
    with steps of ratios 2 and 5.

    """
    start = int(start)
    yield start
    while True:
        for mul in (2, 5, 10):
            yield start * mul
        start *= 10


def preferred_sequence(start, stop, style="nice"):
    """
    Return a sequence of integers with a "preferred" stepping pattern.

    Parameters
    ----------
    start : int
        Starting value in the progression.
    stop : int
        Upper bound of progression, inclusive. Values will not exceed this.
    style : {'nice', 'binary'}
        Style of progression. 'nice' gives geometric steps of 10 with 2 and 5
        in between. 'binary' gives geometric steps of 2.

    Returns
    ------
    list of int

    Examples
    --------
    For certain values of `start` (n * 10^i), nice stepping produces familiar
    "preferred" sequences [1]_:

    Note denominations in Dollars (1-2-5)

        >>> preferred_sequence(1, 100, 'nice')
        [1, 2, 5, 10, 20, 50, 100]


    Coin denominations in Cents

        >>> preferred_sequence(5, 100, 'nice')
        [5, 10, 25, 50, 100]

    .. [1] https://en.wikipedia.org/wiki/Preferred_number#1-2-5_series

    """
    if start > stop:
        return []

    if style == "binary":
        gen = geomprog(start, 2)
    elif style == "nice":
        gen = niceprog(start)
    else:
        ValueError(f"Expected style value of 'binary' or 'nice'; got '{style}'.")

    seq = [next(gen)]
    while True:
        n = next(gen)
        if n > stop:
            break
        seq.append(n)

    return seq


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
                f"Resolution {resn[i]} cannot be derived from "
                f"the base resolutions: {bases}."
            )

    return resn, pred, mult


class CoolerCoarsener(ContactBinner):
    """
    Implementation of cooler coarsening.

    """

    def __init__(self, source_uri, factor, chunksize, columns, agg, batchsize, map=map):
        from .api import Cooler

        self._map = map

        self.source_uri = source_uri
        self.batchsize = batchsize

        assert isinstance(factor, int) and factor > 1
        self.factor = factor
        self.chunksize = int(chunksize)

        self.index_columns = ["bin1_id", "bin2_id"]
        self.value_columns = list(columns)
        self.columns = self.index_columns + self.value_columns
        self.agg = {col: "sum" for col in self.value_columns}
        if agg is not None:
            self.agg.update(agg)

        clr = Cooler(source_uri)
        chromsizes = clr.chromsizes

        # Info for the old bin segmentation
        self.old_binsize = clr.binsize
        self.old_chrom_offset = clr._load_dset("indexes/chrom_offset")
        self.old_bin1_offset = clr._load_dset("indexes/bin1_offset")

        # Calculate the new bin segmentation
        if self.old_binsize is None:
            self.new_binsize = None
        else:
            self.new_binsize = self.old_binsize * factor
        old_bins = clr.bins()[["chrom", "start", "end"]][:]
        self.new_bins = self.coarsen_bins(old_bins, chromsizes, factor)
        self.gs = GenomeSegmentation(chromsizes, self.new_bins)

        # Pre-compute the partition of bin1 offsets that groups the pixels into
        # coarsened bins along the i-axis. Then remove some of the internal
        # edges of this partition to make bigger groups of pixels. This way
        # we ensure that none of the original groups gets split.
        edges = []
        for _chrom, i in self.gs.idmap.items():
            # Respect chrom1 boundaries
            c0 = self.old_chrom_offset[i]
            c1 = self.old_chrom_offset[i + 1]
            edges.extend(self.old_bin1_offset[c0:c1:factor])
        edges.append(self.old_bin1_offset[-1])
        self.edges = _greedy_prune_partition(edges, self.chunksize)

    @staticmethod
    def coarsen_bins(old_bins, chromsizes, factor):
        def _each(group):
            out = group[["chrom", "start"]].copy().iloc[::factor]
            end = group["end"].iloc[factor - 1 :: factor].values
            if len(end) < len(out):
                end = np.r_[end, chromsizes[group.name]]
            out["end"] = end
            return out

        return old_bins.groupby("chrom").apply(_each).reset_index(drop=True)

    def _aggregate(self, span):
        from .api import Cooler

        lo, hi = span

        clr = Cooler(self.source_uri)
        # convert_enum=False returns chroms as raw ints
        table = clr.pixels(join=True, convert_enum=False)[self.columns]
        chunk = table[lo:hi]
        logger.info(f"{lo} {hi}")

        # use the "start" point as anchor for re-binning
        binsize = self.gs.binsize
        chrom_binoffset = self.gs.chrom_binoffset
        chrom_abspos = self.gs.chrom_abspos
        start_abspos = self.gs.start_abspos

        chrom_id1 = chunk["chrom1"].values
        chrom_id2 = chunk["chrom2"].values
        start1 = chunk["start1"].values
        start2 = chunk["start2"].values
        if binsize is None:
            abs_start1 = chrom_abspos[chrom_id1] + start1
            abs_start2 = chrom_abspos[chrom_id2] + start2
            chunk["bin1_id"] = (
                np.searchsorted(start_abspos, abs_start1, side="right") - 1
            )
            chunk["bin2_id"] = (
                np.searchsorted(start_abspos, abs_start2, side="right") - 1
            )
        else:
            rel_bin1 = np.floor(start1 / binsize).astype(int)
            rel_bin2 = np.floor(start2 / binsize).astype(int)
            chunk["bin1_id"] = chrom_binoffset[chrom_id1] + rel_bin1
            chunk["bin2_id"] = chrom_binoffset[chrom_id2] + rel_bin2

        return (
            chunk.groupby(self.index_columns, sort=True)
            .aggregate(self.agg)
            .reset_index()
        )

    def aggregate(self, span):
        try:
            chunk = self._aggregate(span)
        except MemoryError as e:  # pragma: no cover
            raise RuntimeError(str(e)) from e
        return chunk

    def __iter__(self):
        # Distribute batches of `batchsize` pixel spans at once.
        batchsize = self.batchsize
        spans = list(zip(self.edges[:-1], self.edges[1:]))
        for i in range(0, len(spans), batchsize):
            try:
                if batchsize > 1:
                    lock.acquire()
                results = self._map(self.aggregate, spans[i : i + batchsize])
            finally:
                if batchsize > 1:
                    lock.release()
            for df in results:
                yield {k: v.values for k, v in df.items()}


def coarsen_cooler(
    base_uri,
    output_uri,
    factor,
    chunksize,
    nproc=1,
    columns=None,
    dtypes=None,
    agg=None,
    **kwargs,
):
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
    # TODO: decide whether to default to 'count' or whatever is there besides
    # bin1_id, bin2_id dtypes = dict(clr.pixels().dtypes.drop(['bin1_id', 'bin2_id']))

    from .api import Cooler

    clr = Cooler(base_uri)

    factor = int(factor)

    if columns is None:
        columns = ["count"]

    if dtypes is None:
        dtypes = {}

    input_dtypes = clr.pixels().dtypes
    for col in columns:
        if col not in input_dtypes:
            raise ValueError(
                f"Pixel value column '{col}' not found in "
                f"input '{clr.filename}'."
            )
        else:
            dtypes.setdefault(col, input_dtypes[col])

    try:
        # Note: fork before opening to prevent inconsistent global HDF5 state
        if nproc > 1:
            pool = mp.Pool(nproc)
            kwargs.setdefault("lock", lock)

        iterator = CoolerCoarsener(
            base_uri,
            factor,
            chunksize,
            columns=columns,
            agg=agg,
            batchsize=nproc,
            map=pool.map if nproc > 1 else map,
        )

        new_bins = iterator.new_bins

        kwargs.setdefault("append", True)

        create(
            output_uri,
            new_bins,
            iterator,
            dtypes=dtypes,
            symmetric_upper=clr.storage_mode == "symmetric-upper",
            **kwargs,
        )

    finally:
        if nproc > 1:
            pool.close()


def zoomify_cooler(
    base_uris,
    outfile,
    resolutions,
    chunksize,
    nproc=1,
    columns=None,
    dtypes=None,
    agg=None,
    **kwargs,
):
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

    if isinstance(base_uris, str):
        base_uris = [base_uris]

    parsed_uris = {}
    n_bins_longest_chrom = {}
    base_resolutions = set()
    for input_uri in base_uris:
        infile, ingroup = parse_cooler_uri(input_uri)
        clr = Cooler(infile, ingroup)
        base_binsize = 1 if clr.binsize is None else clr.binsize
        parsed_uris[base_binsize] = (infile, ingroup)
        n_bins_longest_chrom[base_binsize] = clr.bins()[:].groupby("chrom").size().max()
        base_resolutions.add(base_binsize)

    # Determine the sequence of reductions.
    resn, pred, mult = get_multiplier_sequence(resolutions, base_resolutions)

    n_zooms = len(resn)

    logger.info(f"Copying base matrices and producing {n_zooms} new zoom levels.")

    if columns is None:
        columns = ["count"]

    # Copy base matrix
    for base_binsize in base_resolutions:
        logger.info("Bin size: " + str(base_binsize))
        infile, ingroup = parsed_uris[base_binsize]
        with h5py.File(infile, "r") as src, \
             h5py.File(outfile, "w") as dest:  # fmt: skip
            prefix = f"/resolutions/{base_binsize}"

            src.copy(ingroup + "/chroms", dest, prefix + "/chroms")
            src.copy(ingroup + "/bins", dest, prefix + "/bins")
            for col in ["bin1_id", "bin2_id", *list(columns)]:
                src.copy(
                    ingroup + f"/pixels/{col}",
                    dest,
                    prefix + f"/pixels/{col}",
                )
            src.copy(ingroup + "/indexes", dest, prefix + "/indexes")
            dest[prefix].attrs.update(src[ingroup].attrs)

    # Aggregate
    # Use lock to sync read/write ops on same file
    for i in range(n_zooms):
        if pred[i] == -1:
            continue
        prev_binsize = resn[pred[i]]
        binsize = prev_binsize * mult[i]
        logger.info(f"Aggregating from {prev_binsize} to {binsize}.")
        coarsen_cooler(
            outfile + f"::resolutions/{prev_binsize}",
            outfile + f"::resolutions/{binsize}",
            mult[i],
            chunksize,
            nproc=nproc,
            columns=columns,
            dtypes=dtypes,
            agg=agg,
            mode="r+",
            **kwargs,
        )

    with h5py.File(outfile, "r+") as fw:
        fw.attrs.update(
            {"format": "HDF5::MCOOL", "format-version": __format_version_mcool__}
        )


def legacy_zoomify(input_uri, outfile, nproc, chunksize, lock=None):
    """
    Quad-tree tiling using legacy MCOOL layout (::0, ::1, ::2, etc.).

    """
    from .api import Cooler

    infile, ingroup = parse_cooler_uri(input_uri)

    clr = Cooler(infile, ingroup)
    n_zooms = get_quadtree_depth(clr.chromsizes, clr.binsize, HIGLASS_TILE_DIM)
    factor = 2

    logger.info(f"total_length (bp): {np.sum(clr.chromsizes)}")
    logger.info(f"binsize: {clr.binsize}")
    logger.info(f"n_zooms: {n_zooms}")
    logger.info(f"quad tile cover: {2 ** n_zooms}")
    logger.info(
        "Copying base matrix to level "
        + f"{n_zooms} and producing {n_zooms} new zoom levels "
        + "counting down to 0..."
    )

    zoom_levels = OrderedDict()
    zoomLevel = str(n_zooms)
    binsize = clr.binsize
    logger.info("Zoom level: " + str(zoomLevel) + " bin size: " + str(binsize))

    # Copy base matrix
    with h5py.File(infile, "r") as src, \
         h5py.File(outfile, "w") as dest:  # fmt: skip

        src.copy(ingroup, dest, str(zoomLevel))
        zoom_levels[zoomLevel] = binsize

    # Aggregate
    # Use lock to sync read/write ops on same file
    for i in range(n_zooms - 1, -1, -1):
        # prev_binsize = binsize
        binsize *= factor
        prevLevel = str(i + 1)
        zoomLevel = str(i)
        logger.info(
            "Aggregating at zoom level: "
            + str(zoomLevel)
            + " bin size: "
            + str(binsize)
        )

        coarsen_cooler(
            outfile + "::" + str(prevLevel),
            outfile + "::" + str(zoomLevel),
            factor,
            chunksize=chunksize,
            nproc=nproc,
            lock=lock,
        )
        zoom_levels[zoomLevel] = binsize

    with h5py.File(outfile, "r+") as fw:
        fw.attrs.update({"max-zoom": n_zooms})
        # grp = fw.require_group('.zooms')
        fw.attrs["max-zooms"] = n_zooms
        fw.attrs.update(zoom_levels)

    return n_zooms, zoom_levels
