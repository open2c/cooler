# -*- coding: utf-8 -*-
from __future__ import division, print_function
from multiprocess import Pool, Lock
from functools import partial
from operator import add
import warnings
import six

import numpy as np
import pandas
import h5py

from ._logging import get_logger
from .api import Cooler
from .tools import split, partition
from .util import mad

__all__ = ['balance_cooler']

logger = get_logger(__name__)


class ConvergenceWarning(UserWarning):
    pass


def _init(chunk):
    return np.copy(chunk['pixels']['count'])


def _binarize(chunk, data):
    data[data != 0] = 1
    return data


def _zero_diags(n_diags, chunk, data):
    pixels = chunk['pixels']
    mask = np.abs(pixels['bin1_id'] - pixels['bin2_id']) < n_diags
    data[mask] = 0
    return data


def _zero_trans(chunk, data):
    chrom_ids = chunk['bins']['chrom']
    pixels = chunk['pixels']
    mask = chrom_ids[pixels['bin1_id']] != chrom_ids[pixels['bin2_id']]
    data[mask] = 0
    return data


def _zero_cis(chunk, data):
    chrom_ids = chunk['bins']['chrom']
    pixels = chunk['pixels']
    mask = chrom_ids[pixels['bin1_id']] == chrom_ids[pixels['bin2_id']]
    data[mask] = 0
    return data


def _timesouterproduct(vec, chunk, data):
    pixels = chunk['pixels']
    data = (vec[pixels['bin1_id']]
                * vec[pixels['bin2_id']]
                * data)
    return data


def _marginalize(chunk, data):
    n = len(chunk['bins']['chrom'])
    pixels = chunk['pixels']
    marg = (
          np.bincount(pixels['bin1_id'], weights=data, minlength=n)
        + np.bincount(pixels['bin2_id'], weights=data, minlength=n)
    )
    return marg


def _balance_genomewide(bias, clr, spans, filters, chunksize, map, tol, max_iters,
                        rescale_marginals, use_lock):
    scale = 1.0
    n_bins = len(bias)

    for _ in range(max_iters):
        marg = (
            split(clr, spans=spans, map=map, use_lock=use_lock)
                .prepare(_init)
                .pipe(filters)
                .pipe(_timesouterproduct, bias)
                .pipe(_marginalize)
                .reduce(add, np.zeros(n_bins))
        )

        nzmarg = marg[marg != 0]
        if not len(nzmarg):
            scale = np.nan
            bias[:] = np.nan
            var = 0.0
            break

        marg = marg / nzmarg.mean()
        marg[marg == 0] = 1
        bias /= marg

        var = nzmarg.var()
        logger.info("variance is {}".format(var))
        if var < tol:
            break
    else:
        warnings.warn(
            'Iteration limit reached without convergence.',
            ConvergenceWarning)

    scale = nzmarg.mean()
    bias[bias == 0] = np.nan
    if rescale_marginals:
        bias /= np.sqrt(scale)

    return bias, scale, var


def _balance_cisonly(bias, clr, spans, filters, chunksize, map, tol, max_iters,
                     rescale_marginals, use_lock):
    chroms = clr.chroms()['name'][:]
    chrom_ids = np.arange(len(clr.chroms()))
    chrom_offsets = clr._load_dset('indexes/chrom_offset')
    bin1_offsets = clr._load_dset('indexes/bin1_offset')
    scales = np.ones(len(chrom_ids))
    n_bins = len(bias)

    for cid, lo, hi in zip(chrom_ids, chrom_offsets[:-1], chrom_offsets[1:]):
        logger.info(chroms[cid])

        plo, phi = bin1_offsets[lo], bin1_offsets[hi]
        spans = list(partition(plo, phi, chunksize))
        scale = 1.0
        for _ in range(max_iters):
            marg = (
                split(clr, spans=spans, map=map, use_lock=use_lock)
                    .prepare(_init)
                    .pipe(filters)
                    .pipe(_timesouterproduct, bias)
                    .pipe(_marginalize)
                    .reduce(add, np.zeros(n_bins))
            )

            marg = marg[lo:hi]
            nzmarg = marg[marg != 0]
            if not len(nzmarg):
                scale = np.nan
                bias[lo:hi] = np.nan
                var = 0.0
                break

            marg = marg / nzmarg.mean()
            marg[marg == 0] = 1
            bias[lo:hi] /= marg

            var = nzmarg.var()
            logger.info("variance is {}".format(var))
            if var < tol:
                break

        else:
            warnings.warn(
                'Iteration limit reached without convergence on {}.'.format(
                    chroms[cid]),
                ConvergenceWarning)

        scale = nzmarg.mean()
        b = bias[lo:hi]
        b[b == 0] = np.nan
        scales[cid] = scale
        if rescale_marginals:
            bias[lo:hi] /= np.sqrt(scale)

    return bias, scales, var


def _balance_transonly(bias, clr, spans, filters, chunksize, map, tol, max_iters,
                        rescale_marginals, use_lock):
    scale = 1.0
    n_bins = len(bias)

    chrom_offsets = clr._load_dset('indexes/chrom_offset')
    cweights = 1. / np.concatenate([
        [(1 - (hi - lo)/n_bins)] * (hi - lo) for lo, hi in
            zip(chrom_offsets[:-1], chrom_offsets[1:])
    ])

    for _ in range(max_iters):
        marg = (
            split(clr, spans=spans, map=map, use_lock=use_lock)
                .prepare(_init)
                .pipe(filters)
                .pipe(_zero_cis)
                .pipe(_timesouterproduct, bias * cweights)
                .pipe(_marginalize)
                .reduce(add, np.zeros(n_bins))
        )

        nzmarg = marg[marg != 0]
        if not len(nzmarg):
            scale = np.nan
            bias[:] = np.nan
            var = 0.0
            break

        marg = marg / nzmarg.mean()
        marg[marg == 0] = 1
        bias /= marg

        var = nzmarg.var()
        logger.info("variance is {}".format(var))
        if var < tol:
            break
    else:
        warnings.warn(
            'Iteration limit reached without convergence.',
            ConvergenceWarning)

    scale = nzmarg.mean()
    bias[bias == 0] = np.nan
    if rescale_marginals:
        bias /= np.sqrt(scale)

    return bias, scale, var



def balance_cooler(clr, chunksize=None, map=map, tol=1e-5,
                   min_nnz=0, min_count=0, mad_max=0,
                   cis_only=False, trans_only=False, ignore_diags=False,
                   max_iters=200, rescale_marginals=True,
                   use_lock=False, blacklist=None, x0=None,
                   store=False, store_name='weight'):
    """
    Iterative correction or matrix balancing of a sparse Hi-C contact map in
    Cooler HDF5 format.

    Parameters
    ----------
    clr : cooler.Cooler
        Cooler object
    chunksize : int, optional
        Split the contact matrix pixel records into equally sized chunks to
        save memory and/or parallelize. Default is to use all the pixels at
        once.
    map : callable, optional
        Map function to dispatch the matrix chunks to workers.
        Default is the builtin ``map``, but alternatives include parallel map
        implementations from a multiprocessing pool.
    tol : float, optional
        Convergence criterion is the variance of the marginal (row/col) sum
        vector.
    min_nnz : int, optional
        Pre-processing bin-level filter. Drop bins with fewer nonzero elements
        than this value.
    min_count : int, optional
        Pre-processing bin-level filter. Drop bins with lower marginal sum than
        this value.
    mad_max : int, optional
        Pre-processing bin-level filter. Drop bins whose log marginal sum is
        less than ``mad_max`` median absolute deviations below the median log
        marginal sum.
    cis_only : bool, optional
        Do iterative correction on intra-chromosomal data only.
        Inter-chromosomal data is ignored.
    trans_only : bool, optional
        Do iterative correction on inter-chromosomal data only.
        Intra-chromosomal data is ignored.
    blacklist : list or 1D array, optional
        An explicit list of IDs of bad bins to filter out when performing
        balancing.
    ignore_diags : int or False, optional
        Drop elements occurring on the first ``ignore_diags`` diagonals of the
        matrix (including the main diagonal).
    max_iters : int, optional
        Iteration limit.
    rescale_marginals : bool, optional
        Normalize the balancing weights such that the balanced matrix has rows /
        columns that sum to 1.0. The scale factor is stored in the ``stats``
        output dictionary.
    x0 : 1D array, optional
        Initial weight vector to use. Default is to start with ones(n_bins).
    store : bool, optional
        Whether to store the results in the file when finished. Default is False.
    store_name : str, optional
        Name of the column of the bin table to save to. Default name is 'weight'.

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
    # Divide the number of elements into non-overlapping chunks
    nnz = clr.info['nnz']
    if chunksize is None:
        chunksize = nnz
        spans = [(0, nnz)]
    else:
        edges = np.arange(0, nnz+chunksize, chunksize)
        spans = list(zip(edges[:-1], edges[1:]))

    # List of pre-marginalization data transformations
    base_filters = []
    if cis_only:
        base_filters.append(_zero_trans)
    if ignore_diags:
        base_filters.append(partial(_zero_diags, ignore_diags))

    # Initialize the bias weights
    n_bins = clr.info['nbins']
    if x0 is not None:
        bias = x0
        bias[np.isnan(bias)] = 0
    else:
        bias = np.ones(n_bins, dtype=float)

    # Drop bins with too few nonzeros from bias
    if min_nnz > 0:
        filters = [_binarize] + base_filters
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
        offsets = clr._load_dset('indexes/chrom_offset')
        for lo, hi in zip(offsets[:-1], offsets[1:]):
            c_marg = marg[lo:hi]
            marg[lo:hi] /= np.median(c_marg[c_marg > 0])
        logNzMarg = np.log(marg[marg>0])
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
            bias, clr, spans, base_filters, chunksize, map, tol, max_iters,
            rescale_marginals, use_lock)
    elif trans_only:
        bias, scale, var = _balance_transonly(
            bias, clr, spans, base_filters, chunksize, map, tol, max_iters,
            rescale_marginals, use_lock)
    else:
        bias, scale, var = _balance_genomewide(
            bias, clr, spans, base_filters, chunksize, map, tol, max_iters,
            rescale_marginals, use_lock)

    stats = {
        'tol': tol,
        'min_nnz': min_nnz,
        'min_count': min_count,
        'mad_max': mad_max,
        'cis_only': cis_only,
        'ignore_diags': ignore_diags,
        'scale': scale,
        'converged': var < tol,
        'var': var,
    }

    if store:
        with clr.open('r+') as grp:
            if store_name in grp['bins']:
                del grp['bins'][store_name]
            h5opts = dict(compression='gzip', compression_opts=6)
            grp['bins'].create_dataset(store_name, data=bias, **h5opts)
            grp['bins'][store_name].attrs.update(stats)

    return bias, stats


iterative_correction = balance_cooler  # alias
