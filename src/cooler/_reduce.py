from __future__ import annotations

import os.path as op
import warnings
from bisect import bisect_right
from collections import defaultdict
from collections.abc import Iterator
from typing import Any, Literal

import h5py
import multiprocess as mp
import numpy as np
import pandas as pd

from ._logging import get_logger
from ._typing import MapFunctor
from ._version import __format_version_mcool__
from .api import Cooler
from .create import ContactBinner, create
from .create._writer import _COMBINE_AGGS, _value_aggs
from .parallel import get_lock
from .util import GenomeSegmentation, parse_cooler_uri

__all__ = ["coarsen_cooler", "merge_coolers", "zoomify_cooler"]


logger = get_logger("cooler.reduce")

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


def merge_breakpoints(
    indexes: list[np.ndarray | h5py.Dataset],
    bufsize: int
) -> tuple[np.ndarray, np.ndarray]:
    """
    Given ``k`` bin1_offset indexes, determine how to partition the data from
    ``k`` corresponding pixel tables for a k-way merge.

    The paritition is a subsequence of bin1 IDs, defining the bounds of chunks
    of data that will be loaded into memory from each table in a single "epoch"
    of merging data. The bounds are calculated such that no single epoch will
    load more than ``bufsize`` records into memory.

    However, the ``bufsize`` condition is not guaranteed (see Notes): check
    ``cum_nrecords``.

    Parameters
    ----------
    indexes : sequence of 1D array-like of equal length
        Offset arrays that map bin1 IDs to their offset locations in a
        corresponding pixel table.
    bufsize : int
        Maximum number of pixel records loaded into memory in a single merge
        epoch.

    Returns
    -------
    bin1_partition : 1D array
        Bin1 IDs defining where to partition all the tables for merging.
    cum_nrecords : 1D array
        Cumulative number of records (from all pixel tables combined) that will
        be processed at each epoch.

    Notes
    -----
    The one exception to the post-condition is when a single bin1 increment in
    a table contains more than ``bufsize`` records.
    """
    # This is a "virtual" cumulative index if all the tables were concatenated
    # and sorted but no pixel records were aggregated. It helps us track how
    # many records would be processed at each merge epoch.
    # NOTE: We sum these incrementally in case the indexes are lazy to avoid
    # loading all indexes into memory at once.
    combined_index = np.zeros(indexes[0].shape, dtype=np.int64)
    for i in range(len(indexes)):
        combined_index += indexes[i]
    combined_start = 0
    combined_nnz = combined_index[-1]

    bin1_partition = [0]
    cum_nrecords = [0]
    lo = 0
    while True:
        # Find the last bin1 ID whose records still fit in the buffer
        target = min(combined_start + bufsize, combined_nnz)
        hi = bisect_right(combined_index, target, lo=lo) - 1
        if combined_index[hi] <= combined_start:
            # A single bin1 has more than bufsize records: take it anyway
            hi = bisect_right(combined_index, combined_start, lo=lo)

        hi = min(hi, len(combined_index) - 1)
        bin1_partition.append(hi)
        cum_nrecords.append(combined_index[hi])

        if combined_index[hi] == combined_nnz:
            break

        lo = hi
        combined_start = combined_index[hi]

    bin1_partition = np.array(bin1_partition)
    cum_nrecords = np.array(cum_nrecords)
    return bin1_partition, cum_nrecords


# Aggregations of value columns that give the same result however the pixels
# are grouped: in one step or level by level when coarsening, and input by
# input when merging
_AGGS = _COMBINE_AGGS


class CoolerMerger(ContactBinner):
    """
    Implementation of cooler merging.

    """

    def __init__(
        self,
        coolers: list[Cooler],
        mergebuf: int,
        columns: list[str] | None = None,
        agg: dict[str, Any] | None = None
    ):
        self.coolers = list(coolers)
        self.mergebuf = mergebuf
        self.columns = ["count"] if columns is None else columns
        self.agg = _value_aggs(self.columns, agg)

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

    def __iter__(self) -> Iterator[dict[str, np.ndarray]]:
        # Load the bin1_offset indexes, without keeping the files open.
        indexes = []
        for c in self.coolers:
            with c.open("r") as grp:
                indexes.append(grp["indexes/bin1_offset"][:])

        # Calculate the common partition of bin1 offsets that define the epochs
        # of merging data.
        bin1_partition, cum_nrecords = merge_breakpoints(indexes, self.mergebuf)
        nrecords_per_epoch = np.diff(cum_nrecords)
        if len(nrecords_per_epoch) and nrecords_per_epoch.max() > self.mergebuf:
            warnings.warn(
                f"Merging will buffer up to {nrecords_per_epoch.max():,} pixels at "
                f"once, more than mergebuf={self.mergebuf:,}, because a row of "
                "the matrix can't be split across merge epochs.",
                stacklevel=2,
            )
        nnzs = [len(c.pixels()) for c in self.coolers]
        logger.info(
            f"Merging {sum(nnzs):,} pixels from {len(self.coolers)} coolers in "
            f"{len(nrecords_per_epoch)} epochs"
        )
        logger.debug(f"Pixels per cooler: {nnzs}")
        logger.debug(f"Epoch boundaries (bin1 IDs): {bin1_partition}")
        logger.debug(f"Pixels per epoch: {nrecords_per_epoch}")

        starts = [0] * len(self.coolers)
        for bin1_id in bin1_partition[1:]:
            stops = [index[bin1_id] for index in indexes]

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

            logger.debug(f"Pixels merged up to: {stops}")
            starts = stops


def merge_coolers(
    output_uri: str,
    input_uris: list[str],
    mergebuf: int,
    columns: list[str] | None = None,
    dtypes: dict[str, Any] | None = None,
    agg: dict[str, Any] | None = None,
    **kwargs
) -> None:
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
        How to combine the values of a pixel across the input coolers: a
        mapping of value column names to "sum" (default), "min" or "max". An
        input without the pixel doesn't contribute a value.
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
    logger.info(f"Merging {len(input_uris)} coolers into {output_uri}")
    for uri in input_uris:
        logger.debug(f"Input: {uri}")

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

    dtypes = {} if dtypes is None else dict(dtypes)
    for col in columns:
        if col not in dtypes:
            # Find a common dtype to accomodate all the inputs being merged.
            dtypes[col] = np.result_type(*dtype_map[col])

    bins = clrs[0].bins()[["chrom", "start", "end"]][:]
    assembly = clrs[0].info.get("genome-assembly", None)
    iterator = CoolerMerger(clrs, mergebuf=mergebuf, columns=columns, agg=agg)

    create(
        output_uri,
        bins,
        iterator,
        columns=columns,
        dtypes=dtypes,
        assembly=assembly,
        symmetric_upper=symmetric_upper,
        binsize=clrs[0].binsize,
        **kwargs,
    )


def _greedy_prune_partition(edges: np.ndarray, maxlen: int) -> np.ndarray:
    """Given an integer interval partition ``edges`` from 0..nnz, prune the
    edges to make the new subintervals roughly ``maxlen`` in length.
    """
    edges = np.asarray(edges)
    assert len(edges) >= 2 and edges[0] == 0
    cumlen = np.r_[0, np.cumsum(np.diff(edges))]
    cuts = [maxlen * i for i in range(int(np.ceil(cumlen[-1] / maxlen)))]
    cuts.append(cumlen[-1])
    idx = np.unique(np.searchsorted(cumlen, cuts))
    return edges[idx]


def geomprog(start: int, mul: int) -> Iterator[int]:
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


def niceprog(start: int) -> Iterator[int]:
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


def preferred_sequence(
    start: int,
    stop: int,
    style: Literal["binary", "nice"] = "nice"
) -> list[int]:
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
        raise ValueError(f"Expected style value of 'binary' or 'nice'; got '{style}'.")

    seq = [next(gen)]
    while True:
        n = next(gen)
        if n > stop:
            break
        seq.append(n)

    return seq


def get_multiplier_sequence(
    resolutions: list[int],
    bases: list[int] | None = None
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
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
        if target in bases:
            continue  # data for base resolutions already exists
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

    def __init__(
        self,
        source_uri: str,
        factor: int,
        chunksize: int,
        columns: list[str],
        agg: dict[str, Any] | None,
        batchsize: int,
        map: MapFunctor = map,
    ):
        self._map = map

        self.source_uri = source_uri
        self.batchsize = batchsize

        assert isinstance(factor, int) and factor > 1
        self.factor = factor
        self.chunksize = int(chunksize)

        self.index_columns = ["bin1_id", "bin2_id"]
        self.value_columns = list(columns)
        self.columns = self.index_columns + self.value_columns
        self.agg = _value_aggs(self.value_columns, agg)

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
        for i in self.gs.idmap.to_numpy():
            # Respect chrom1 boundaries
            c0 = self.old_chrom_offset[i]
            c1 = self.old_chrom_offset[i + 1]
            edges.extend(self.old_bin1_offset[c0:c1:factor])
        edges.append(self.old_bin1_offset[-1])
        self.edges = _greedy_prune_partition(edges, self.chunksize)

    @staticmethod
    def coarsen_bins(
        old_bins: pd.DataFrame, chromsizes: pd.Series, factor: int
    ) -> pd.DataFrame:
        def _each(group):
            out = group[["chrom", "start"]].copy().iloc[::factor]
            end = group["end"].iloc[factor - 1 :: factor].values
            if len(end) < len(out):
                end = np.r_[end, chromsizes[group.name]]
            out["end"] = end
            return out

        return (
            old_bins
            .groupby("chrom", observed=True)[["chrom", "start", "end"]]
            .apply(_each)
            .reset_index(drop=True)
        )

    def aggregate(self, span: tuple[int, int]) -> pd.DataFrame:
        lo, hi = span

        clr = Cooler(self.source_uri)
        # convert_enum=False returns chroms as raw ints
        table = clr.pixels(join=True, convert_enum=False)[self.columns]
        chunk = table[lo:hi]
        logger.debug(f"Coarsening pixels {lo:,}-{hi:,}")

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

    def __iter__(self) -> Iterator[dict[str, np.ndarray]]:
        # Distribute batches of `batchsize` pixel spans at once.
        batchsize = self.batchsize
        spans = list(zip(self.edges[:-1], self.edges[1:]))
        for i in range(0, len(spans), batchsize):
            try:
                if batchsize > 1:
                    get_lock().acquire()
                results = self._map(self.aggregate, spans[i : i + batchsize])
            finally:
                if batchsize > 1:
                    get_lock().release()
            for df in results:
                yield {k: v.values for k, v in df.items()}


def coarsen_cooler(
    base_uri: str,
    output_uri: str,
    factor: int,
    chunksize: int,
    nproc: int = 1,
    columns: list[str] | None = None,
    dtypes: dict[str, Any] | None = None,
    agg: dict[str, Any] | None = None,
    **kwargs,
) -> None:
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
        How to combine the pixels of each *k*-by-*k* block: a mapping of value
        column names to "sum" (default), "min" or "max". Pixels that aren't
        stored don't contribute a value.
    kwargs
        Passed to ``cooler.create``.

    See also
    --------
    cooler.zoomify_cooler
    cooler.merge_coolers

    """
    # TODO: decide whether to default to 'count' or whatever is there besides
    # bin1_id, bin2_id dtypes = dict(clr.pixels().dtypes.drop(['bin1_id', 'bin2_id']))
    clr = Cooler(base_uri)

    factor = int(factor)
    if factor < 2:
        raise ValueError(f"The coarsening factor must be at least 2, not {factor}.")

    if columns is None:
        columns = ["count"]

    dtypes = {} if dtypes is None else dict(dtypes)

    input_dtypes = clr.pixels().dtypes
    for col in columns:
        if col not in input_dtypes:
            raise ValueError(
                f"Pixel value column '{col}' not found in " f"input '{clr.filename}'."
            )
        else:
            dtypes.setdefault(col, input_dtypes[col])

    try:
        # Note: fork before opening to prevent inconsistent global HDF5 state
        if nproc > 1:
            kwargs.setdefault("lock", get_lock())
            pool = mp.Pool(nproc)

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

        kwargs.setdefault("mode", "a")
        kwargs.setdefault("assembly", clr.info.get("genome-assembly"))
        kwargs.setdefault("metadata", clr.info.get("metadata"))

        create(
            output_uri,
            new_bins,
            iterator,
            columns=columns,
            dtypes=dtypes,
            symmetric_upper=clr.storage_mode == "symmetric-upper",
            binsize=iterator.new_binsize,
            **kwargs,
        )

    finally:
        if nproc > 1:
            pool.close()
            pool.join()


def zoomify_cooler(
    base_uris: str | list[str],
    outfile: str,
    resolutions: list[int],
    chunksize: int,
    nproc: int = 1,
    columns: list[str] | None = None,
    dtypes: dict[str, Any] | None = None,
    agg: dict[str, Any] | None = None,
    **kwargs,
) -> None:
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
        How to combine the pixels of each block when coarsening: a mapping of
        value column names to "sum" (default), "min" or "max". These give the
        same result whether a resolution is coarsened from the base or from
        the previous resolution.
    kwargs
        Passed to ``cooler.create``.

    See also
    --------
    cooler.coarsen_cooler
    cooler.merge_coolers

    """
    # TODO: provide presets? {'pow2', '4dn'}
    if isinstance(base_uris, str):
        base_uris = [base_uris]

    parsed_uris = {}
    base_resolutions = set()
    for input_uri in base_uris:
        infile, ingroup = parse_cooler_uri(input_uri)
        if op.exists(outfile) and op.samefile(infile, outfile):
            # The output file is truncated before the inputs are read
            raise ValueError(
                f"The output file '{outfile}' contains the input '{input_uri}'."
            )
        clr = Cooler(infile, ingroup)
        base_binsize = 1 if clr.binsize is None else clr.binsize
        if base_binsize in parsed_uris:
            raise ValueError(
                "Base coolers must have different resolutions, and at most one "
                f"can have variable-size bins: '{input_uri}'."
            )
        parsed_uris[base_binsize] = (infile, ingroup)
        base_resolutions.add(base_binsize)

    # Determine the sequence of reductions.
    resn, pred, mult = get_multiplier_sequence(resolutions, base_resolutions)

    n_zooms = len(resn)
    n_new = int(np.sum(pred != -1))
    bases = ", ".join(str(b) for b in sorted(base_resolutions))
    logger.info(f"Copying base resolutions ({bases}) and producing {n_new} new ones")

    if columns is None:
        columns = ["count"]

    # Copy base matrices into a new output file
    with h5py.File(outfile, "w"):
        pass
    for base_binsize in sorted(base_resolutions):
        logger.info(f"Copying base resolution {base_binsize}")
        infile, ingroup = parsed_uris[base_binsize]
        with h5py.File(infile, "r") as src, \
             h5py.File(outfile, "r+") as dest:  # fmt: skip
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
        logger.info(f"Aggregating from {prev_binsize} to {binsize}")
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
