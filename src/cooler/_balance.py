from __future__ import annotations

import warnings
from functools import partial
from operator import add
from typing import Literal

import numpy as np

from ._logging import get_logger
from ._typing import MapFunctor
from .api import Cooler
from .parallel import partition, split
from .util import mad

__all__ = ["balance_cooler"]

logger = get_logger(__name__)


class ConvergenceWarning(UserWarning):
    pass


def _init(chunk):
    return np.copy(chunk["pixels"]["count"])


def _binarize(chunk, data):
    data[data != 0] = 1
    return data


def _zero_diags(n_diags, chunk, data):
    # Only within chromosomes: neighboring bins across a chromosome boundary
    # are trans (#211)
    chrom_ids = chunk["bins"]["chrom"]
    pixels = chunk["pixels"]
    bin1, bin2 = pixels["bin1_id"], pixels["bin2_id"]
    mask = (np.abs(bin1 - bin2) < n_diags) & (chrom_ids[bin1] == chrom_ids[bin2])
    data[mask] = 0
    return data


def _zero_dist(dist, chunk, data):
    # Distance between bin midpoints, which is the same as between starts for
    # uniform bins but doesn't depend on which of two bins is larger
    bins = chunk["bins"]
    chrom_ids = bins["chrom"]
    mids2 = bins["start"] + bins["end"]  # twice the midpoint, to stay integer
    pixels = chunk["pixels"]
    bin1, bin2 = pixels["bin1_id"], pixels["bin2_id"]
    mask = (np.abs(mids2[bin2] - mids2[bin1]) < 2 * dist) & (
        chrom_ids[bin1] == chrom_ids[bin2]
    )
    data[mask] = 0
    return data


def _zero_trans(chunk, data):
    chrom_ids = chunk["bins"]["chrom"]
    pixels = chunk["pixels"]
    mask = chrom_ids[pixels["bin1_id"]] != chrom_ids[pixels["bin2_id"]]
    data[mask] = 0
    return data


def _zero_cis(chunk, data):
    chrom_ids = chunk["bins"]["chrom"]
    pixels = chunk["pixels"]
    mask = chrom_ids[pixels["bin1_id"]] == chrom_ids[pixels["bin2_id"]]
    data[mask] = 0
    return data


def _timesouterproduct(vec, chunk, data):
    pixels = chunk["pixels"]
    data = vec[pixels["bin1_id"]] * vec[pixels["bin2_id"]] * data
    return data


def _marginalize(chunk, data):
    n = len(chunk["bins"]["chrom"])
    pixels = chunk["pixels"]
    # Count main diagonal pixels once, not once per bin1 and bin2
    offdiag = pixels["bin1_id"] != pixels["bin2_id"]
    marg = np.bincount(
        pixels["bin1_id"][offdiag], weights=data[offdiag], minlength=n
    ) + np.bincount(pixels["bin2_id"], weights=data, minlength=n)
    return marg


def _marginals(clr, spans, filters, map, use_lock, bias):
    return (
        split(clr, spans=spans, map=map, use_lock=use_lock)
        .prepare(_init)
        .pipe(filters)
        .pipe(_timesouterproduct, bias)
        .pipe(_marginalize)
        .reduce(add, np.zeros(len(bias)))
    )


def _trans_marginals(clr, spans, filters, map, use_lock, cweights, bias):
    return _marginals(
        clr, spans, [*filters, _zero_cis], map, use_lock, bias * cweights
    )


def _iterate(bias, sel, get_marginals, tol, max_iters, rescale_marginals, where=""):
    """
    Balance the bins ``sel`` of ``bias`` in place.

    Returns the mean and variance of the marginals at the last iteration, which
    are NaN if no iteration was done or if every bin is filtered out. Bins with
    no data left after filtering get NaN weights.

    """
    scale = var = np.nan
    marg = None
    for _ in range(max_iters):
        marg = get_marginals(bias)[sel]
        nzmarg = marg[marg != 0]
        if not len(nzmarg):
            break

        scale = nzmarg.mean()
        var = nzmarg.var()
        factor = marg / scale
        factor[factor == 0] = 1
        bias[sel] /= factor

        logger.info(f"variance is {var}")
        if var < tol:
            break
    else:
        warnings.warn(
            f"Iteration limit reached without convergence{where}.",
            ConvergenceWarning,
            stacklevel=1,
        )

    b = bias[sel]
    if marg is not None:
        b[marg == 0] = 0
    b[b == 0] = np.nan
    if rescale_marginals and np.isfinite(scale):
        b /= np.sqrt(scale)
    bias[sel] = b
    return scale, var


def _balance_genomewide(
    bias,
    clr,
    spans,
    filters,
    chunksize,
    map,
    tol,
    max_iters,
    rescale_marginals,
    use_lock,
):
    get_marginals = partial(
        _marginals, clr, spans, filters, map, use_lock
    )
    scale, var = _iterate(
        bias, slice(None), get_marginals, tol, max_iters, rescale_marginals
    )
    return bias, scale, var


def _balance_cisonly(
    bias,
    clr,
    spans,
    filters,
    chunksize,
    map,
    tol,
    max_iters,
    rescale_marginals,
    use_lock,
):
    chroms = clr.chroms()["name"][:]
    chrom_ids = np.arange(len(clr.chroms()))
    chrom_offsets = clr._load_dset("indexes/chrom_offset")
    bin1_offsets = clr._load_dset("indexes/bin1_offset")
    scales = np.ones(len(chrom_ids))
    variances = np.full_like(scales, np.nan)

    for cid, lo, hi in zip(chrom_ids, chrom_offsets[:-1], chrom_offsets[1:]):
        logger.info(chroms[cid])

        plo, phi = bin1_offsets[lo], bin1_offsets[hi]
        spans = list(partition(plo, phi, chunksize))
        get_marginals = partial(
            _marginals, clr, spans, filters, map, use_lock
        )
        scales[cid], variances[cid] = _iterate(
            bias,
            slice(lo, hi),
            get_marginals,
            tol,
            max_iters,
            rescale_marginals,
            where=f" on {chroms[cid]}",
        )

    # A chromosome with nothing left to balance counts as converged, unless
    # every chromosome is empty
    empty = np.isnan(scales)
    if max_iters > 0 and not empty.all():
        variances[empty] = 0.0

    return bias, scales, variances


def _balance_transonly(
    bias,
    clr,
    spans,
    filters,
    chunksize,
    map,
    tol,
    max_iters,
    rescale_marginals,
    use_lock,
):
    n_bins = len(bias)
    chrom_offsets = clr._load_dset("indexes/chrom_offset")
    cweights = 1.0 / np.concatenate(
        [
            [(1 - (hi - lo) / n_bins)] * (hi - lo)
            for lo, hi in zip(chrom_offsets[:-1], chrom_offsets[1:])
        ]
    )
    get_marginals = partial(
        _trans_marginals, clr, spans, filters, map, use_lock, cweights
    )
    scale, var = _iterate(
        bias, slice(None), get_marginals, tol, max_iters, rescale_marginals
    )
    return bias, scale, var


def balance_cooler(
    clr: Cooler,
    *,
    cis_only: bool = False,
    trans_only: bool = False,
    ignore_diags: int | Literal[False] = 2,
    ignore_dist: int | None = None,
    mad_max: int = 5,
    min_nnz: int = 10,
    min_count: int = 0,
    blacklist: str | None = None,
    rescale_marginals: bool = True,
    x0: np.ndarray | None = None,
    tol: float = 1e-5,
    max_iters: int = 200,
    chunksize: int = 10_000_000,
    map: MapFunctor = map,
    use_lock: bool = False,
    store: bool = False,
    store_name: str = "weight",
) -> tuple[np.ndarray, dict]:
    """
    Iterative correction or matrix balancing of a sparse Hi-C contact map in
    Cooler HDF5 format.

    Parameters
    ----------
    clr : cooler.Cooler
        Cooler object
    cis_only : bool, optional
        Do iterative correction on intra-chromosomal data only.
        Inter-chromosomal data is ignored.
    trans_only : bool, optional
        Do iterative correction on inter-chromosomal data only.
        Intra-chromosomal data is ignored.
    ignore_diags : int or False, optional
        Drop elements occurring on the first ``ignore_diags`` diagonals of the
        matrix (including the main diagonal). Only intra-chromosomal elements
        are dropped.
    ignore_dist : int, optional
        Drop intra-chromosomal elements whose bin midpoints are less than this
        many base pairs apart. Unlike ``ignore_diags``, this also works with
        variable-size bins. Both can be given.
    chunksize : int or None, optional
        Split the contact matrix pixel records into equally sized chunks to
        save memory and/or parallelize. Set to ``None`` to use all the pixels
        at once.
    mad_max : int, optional
        Pre-processing bin-level filter. Each bin's marginal sum is first
        divided by the median marginal sum of its chromosome. Then bins whose
        log normalized marginal sum is less than ``mad_max`` median absolute
        deviations below the genome-wide median are dropped.
    min_nnz : int, optional
        Pre-processing bin-level filter. Drop bins with fewer nonzero elements
        than this value.
    min_count : int, optional
        Pre-processing bin-level filter. Drop bins with lower marginal sum than
        this value.
    blacklist : list or 1D array, optional
        An explicit list of IDs of bad bins to filter out when performing
        balancing.
    rescale_marginals : bool, optional
        Normalize the balancing weights such that the balanced matrix has rows
        / columns that sum to 1.0. The scale factor is stored in the ``stats``
        output dictionary.
    map : callable, optional
        Map function to dispatch the matrix chunks to workers.
        Default is the builtin ``map``, but alternatives include parallel map
        implementations from a multiprocessing pool.
    x0 : 1D array, optional
        Initial weight vector to use. Default is to start with ones(n_bins).
        It is not modified. NaN values are treated as filtered-out bins.
    tol : float, optional
        Convergence criterion is the variance of the marginal (row/col) sum
        vector.
    max_iters : int, optional
        Iteration limit.
    store : bool, optional
        Whether to store the results in the file when finished. Default is
        False.
    store_name : str, optional
        Name of the column of the bin table to save to. Default name is
        'weight'.

    Returns
    -------
    bias : 1D array, whose shape is the number of bins in ``h5``.
        Vector of bin bias weights to normalize the observed contact map.
        Dropped bins will be assigned the value NaN.
        N[i, j] = O[i, j] * bias[i] * bias[j]
    stats : dict
        Summary of parameters used to perform balancing and the average
        magnitude of the corrected matrix's marginal sum at convergence.

    """
    # Balancing with a single weight vector is only defined for symmetric
    # matrices, and the marginals assume the upper triangle is stored
    if clr.storage_mode != "symmetric-upper":
        raise ValueError(
            "Balancing requires a symmetric matrix in symmetric-upper storage, "
            f"but the cooler has '{clr.storage_mode}' storage."
        )

    # Divide the number of elements into non-overlapping chunks
    nnz = int(clr.info["nnz"])
    if chunksize is None:
        chunksize = nnz
        spans = [(0, nnz)]
    else:
        edges = np.arange(0, nnz + chunksize, chunksize)
        spans = list(zip(edges[:-1], edges[1:]))

    # List of pre-marginalization data transformations
    base_filters = []
    if cis_only:
        base_filters.append(_zero_trans)
    if ignore_diags:
        base_filters.append(partial(_zero_diags, ignore_diags))
    if ignore_dist:
        base_filters.append(partial(_zero_dist, ignore_dist))

    # Initialize the bias weights
    n_bins = int(clr.info["nbins"])
    if x0 is not None:
        bias = np.array(x0, dtype=float)
        if bias.shape != (n_bins,):
            raise ValueError(
                f"x0 has shape {bias.shape}, but the cooler has {n_bins} bins."
            )
        bias[np.isnan(bias)] = 0
    else:
        bias = np.ones(n_bins, dtype=float)

    # Drop bins with too few nonzeros from bias
    if min_nnz > 0:
        filters = [_binarize, *base_filters]
        marg_nnz = (
            split(clr, spans=spans, map=map, use_lock=use_lock)
            .prepare(_init)
            .pipe(filters)
            .pipe(_marginalize)
            .reduce(add, np.zeros(n_bins))
        )
        bias[marg_nnz < min_nnz] = 0

    filters = base_filters
    marg = (
        split(clr, spans=spans, map=map, use_lock=use_lock)
        .prepare(_init)
        .pipe(filters)
        .pipe(_marginalize)
        .reduce(add, np.zeros(n_bins))
    )

    # Drop bins with too few total counts from bias
    if min_count:
        bias[marg < min_count] = 0

    # MAD-max filter on the marginals
    if mad_max > 0:
        offsets = clr._load_dset("indexes/chrom_offset")
        for lo, hi in zip(offsets[:-1], offsets[1:]):
            c_marg = marg[lo:hi]
            if (c_marg > 0).any():  # bins of empty chromosomes stay at 0
                marg[lo:hi] /= np.median(c_marg[c_marg > 0])
        logNzMarg = np.log(marg[marg > 0])
        med_logNzMarg = np.median(logNzMarg)
        dev_logNzMarg = mad(logNzMarg)
        cutoff = np.exp(med_logNzMarg - mad_max * dev_logNzMarg)
        bias[marg < cutoff] = 0

    # Filter out pre-determined bad bins
    if blacklist is not None:
        bias[blacklist] = 0

    # Do balancing
    if cis_only:
        bias, scale, var = _balance_cisonly(
            bias,
            clr,
            spans,
            base_filters,
            chunksize,
            map,
            tol,
            max_iters,
            rescale_marginals,
            use_lock,
        )
    elif trans_only:
        bias, scale, var = _balance_transonly(
            bias,
            clr,
            spans,
            base_filters,
            chunksize,
            map,
            tol,
            max_iters,
            rescale_marginals,
            use_lock,
        )
    else:
        bias, scale, var = _balance_genomewide(
            bias,
            clr,
            spans,
            base_filters,
            chunksize,
            map,
            tol,
            max_iters,
            rescale_marginals,
            use_lock,
        )

    stats = {
        "tol": tol,
        "min_nnz": min_nnz,
        "min_count": min_count,
        "mad_max": mad_max,
        "cis_only": cis_only,
        "trans_only": trans_only,
        "ignore_diags": ignore_diags,
        "scale": scale,
        "converged": var < tol,
        "var": var,
        "divisive_weights": False,
    }
    if ignore_dist:
        stats["ignore_dist"] = ignore_dist

    if store:
        from .create._writer import _add_columns

        _add_columns(
            clr.uri,
            "bins",
            {store_name: bias},
            attrs={store_name: stats},
            replace=True,
            h5opts={"compression": "gzip", "compression_opts": 6},
        )

    return bias, stats


iterative_correction = balance_cooler  # alias
