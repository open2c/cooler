# -*- coding: utf-8 -*-
from __future__ import division, print_function
from multiprocess import Pool, Lock
from functools import partial
import warnings
import six

import numpy as np
import pandas
import h5py

from . import get_logger
from .api import Cooler
from .tools import split, partition
from .util import mad


logger = get_logger()


def _init_transform(chunk):
    return chunk, np.copy(chunk['pixels']['count'])


def _binarize_mask(args):
    chunk, data = args
    data[data != 0] = 1
    return chunk, data


def _dropdiag_mask(n_diags, args):
    chunk, data = args
    pixels = chunk['pixels']
    mask = np.abs(pixels['bin1_id'] - pixels['bin2_id']) < n_diags
    data[mask] = 0
    return chunk, data


def _cisonly_mask(args):
    chunk, data = args
    chrom_ids = chunk['bins']['chrom']
    pixels = chunk['pixels']
    mask = chrom_ids[pixels['bin1_id']] != chrom_ids[pixels['bin2_id']]
    data[mask] = 0
    return chunk, data


def _timesouterproduct_transform(vec, args):
    chunk, data = args
    pixels = chunk['pixels']
    data = (vec[pixels['bin1_id']]
                * vec[pixels['bin2_id']]
                * data)
    return chunk, data


def _marginalize_transform(args):
    chunk, data = args
    n = len(chunk['bins']['chrom'])
    pixels = chunk['pixels']
    marg = (
          np.bincount(pixels['bin1_id'], weights=data, minlength=n)
        + np.bincount(pixels['bin2_id'], weights=data, minlength=n)
    )
    return marg


def _balance_genomewide(bias, c, spans, filters, chunksize, map, tol, max_iters,
                        rescale_marginals, use_lock):
    scale = 1.0
    for _ in range(max_iters):
        marg = np.sum(
            split(c, cooler_root=c.root, spans=spans, map=map, use_lock=use_lock)
                .pipe(_init_transform)
                .pipe(filters)
                .pipe(_timesouterproduct_transform, bias)
                .pipe(_marginalize_transform)
                .combine(),
            axis=0)
        
        nzmarg = marg[marg != 0]
        if not len(nzmarg):
            scale = np.nan
            bias[:] = np.nan
            break

        marg = marg / nzmarg.mean()
        marg[marg == 0] = 1
        bias /= marg

        var = nzmarg.var()
        logger.info("variance is {}".format(var))
        if var < tol:
            scale = nzmarg.mean()
            bias[bias == 0] = np.nan
            break
    else:
        raise RuntimeError('Iteration limit reached without convergence.')

    if rescale_marginals:
        bias /= np.sqrt(scale)

    return bias, scale


def _balance_cisonly(bias, c, spans, filters, chunksize, map, tol, max_iters,
                     rescale_marginals, use_lock):
    chroms = c.chroms()['name'][:]
    chrom_ids = np.arange(len(c.chroms()))
    chrom_offsets = c._load_dset('/indexes/chrom_offset')
    bin1_offsets = c._load_dset('/indexes/bin1_offset')
    scales = np.ones(len(chrom_ids))
    
    for cid, lo, hi in zip(chrom_ids, chrom_offsets[:-1], chrom_offsets[1:]):
        logger.info(chroms[cid])

        plo, phi = bin1_offsets[lo], bin1_offsets[hi]
        spans = list(partition(plo, phi, chunksize))
        scale = 1.0
        for _ in range(max_iters):
            marg = np.sum(
                split(c, cooler_root=c.root,  spans=spans, map=map, use_lock=use_lock)
                    .pipe(_init_transform)
                    .pipe(filters)
                    .pipe(_timesouterproduct_transform, bias)
                    .pipe(_marginalize_transform)
                    .combine(),
                axis=0)

            marg = marg[lo:hi]
            nzmarg = marg[marg != 0]
            if not len(nzmarg):
                scale = np.nan
                bias[lo:hi] = np.nan
                break

            marg = marg / nzmarg.mean()
            marg[marg == 0] = 1
            bias[lo:hi] /= marg

            var = nzmarg.var()
            logger.info("variance is {}".format(var))
            if var < tol:
                scale = nzmarg.mean()
                b = bias[lo:hi]
                b[b == 0] = np.nan
                break

        else:
            raise RuntimeError('Iteration limit reached without convergence.')
        
        scales[cid] = scale
        if rescale_marginals:
            bias[lo:hi] /= np.sqrt(scale)
        
    return bias, scales


def iterative_correction(h5, cooler_root='/', chunksize=None, map=map, tol=1e-5,
                         min_nnz=0, min_count=0, mad_max=0,
                         cis_only=False, ignore_diags=False,
                         max_iters=200, rescale_marginals=True,
                         use_lock=True):
    """
    Iterative correction or matrix balancing of a sparse Hi-C contact map in
    Cooler HDF5 format.

    Parameters
    ----------
    h5 : h5py.File object
        Cooler file
    cooler_root : str, optional
        Path of the root node of a cooler tree. Default is the file root, '/'.
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
        less than ``mad_max`` mean absolute deviations below the median log
        marginal sum.
    cis_only: bool, optional
        Do iterative correction on intra-chromosomal data only.
        Inter-chromosomal data is ignored.
    ignore_diags : int or False, optional
        Drop elements occurring on the first ``ignore_diags`` diagonals of the
        matrix (including the main diagonal).
    max_iters : int, optional
        Iteration limit.
    rescale_marginals : bool, optional
        Normalize the balancing weights such that the balanced matrix has rows /
        columns that sum to 1.0. The scale factor is stored in the ``stats``
        output dictionary.

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
    filepath = h5.file.filename
    c = Cooler(filepath, cooler_root)

    # Divide the number of elements into non-overlapping chunks
    nnz = h5[cooler_root].attrs['nnz']
    if chunksize is None:
        chunksize = nnz
        spans = [(0, nnz)]
    else:
        edges = np.arange(0, nnz+chunksize, chunksize)
        spans = list(zip(edges[:-1], edges[1:]))

    # List of pre-marginalization data transformations
    base_filters = []
    if cis_only:
        base_filters.append(_cisonly_mask)
    if ignore_diags:
        base_filters.append(partial(_dropdiag_mask, ignore_diags))

    # Initialize the bias weights
    n_bins = h5[cooler_root].attrs['nbins']
    bias = np.ones(n_bins, dtype=float)

    # Drop bins with too few nonzeros from bias
    if min_nnz > 0:
        filters = [_binarize_mask] + base_filters
        marg_nnz = np.sum(
            split(c, cooler_root=c.root, spans=spans, map=map, use_lock=use_lock)
                .pipe(_init_transform)
                .pipe(filters)
                .pipe(_marginalize_transform)
                .combine(),
            axis=0)
        bias[marg_nnz < min_nnz] = 0

    filters = base_filters
    marg = np.sum(
        split(c, cooler_root=c.root, spans=spans, map=map, use_lock=use_lock)
            .pipe(_init_transform)
            .pipe(filters)
            .pipe(_marginalize_transform)
            .combine(),
        axis=0)

    # Drop bins with too few total counts from bias
    if min_count:
        bias[marg < min_count] = 0

    # MAD-max filter on the marginals
    if mad_max > 0:
        offsets = h5[cooler_root]['indexes']['chrom_offset'][:]
        for lo, hi in zip(offsets[:-1], offsets[1:]):
            c_marg = marg[lo:hi]
            marg[lo:hi] /= np.median(c_marg[c_marg > 0])
        logNzMarg = np.log(marg[marg>0])
        logMedMarg = np.median(logNzMarg)
        madSigma = mad(logNzMarg) / 0.6745
        cutoff = np.exp(logMedMarg - mad_max * madSigma)
        bias[marg < cutoff] = 0

    # Do balancing
    if cis_only:
        bias, scale = _balance_cisonly(
            bias, c, spans, base_filters, chunksize, map, tol, max_iters,
            rescale_marginals, use_lock)
    else:
        bias, scale = _balance_genomewide(
            bias, c, spans, base_filters, chunksize, map, tol, max_iters,
            rescale_marginals, use_lock)

    stats = {
        'tol': tol,
        'min_nnz': min_nnz,
        'min_count': min_count,
        'mad_max': mad_max,
        'cis_only': cis_only,
        'ignore_diags': ignore_diags,
        'scale': scale,
    }

    return bias, stats
