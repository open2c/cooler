# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals
from multiprocessing import Pool, Lock
import six

import numpy as np
import numexpr
import pandas
import h5py


# global lock because libhdf5 builds are not thread-safe by default :(
lock = Lock()


class Worker(object):
    """
    Worker to do partial marginalization of a sparse heatmap in Cooler format.

    A worker fetches a set of records representing matrix elements (pixels) and
    uses a weighted bincount algorithm to compute the partial marginal sum of
    the matrix coming from that set.

    The initial weight assigned to each pixel in a chunk is its value. The
    initializer accepts a sequence of "filters" or transformations to
    sequentially process the pixel weights before aggregating/marginalizing,
    e.g., to remove unwanted elements by setting their weights to zero.

    Parameters
    ----------
    cooler_path : str
        Path to Cooler HDF5 file. Must be a string and not a file-like object
        when running in parallel because the latter cannot be pickled.

    filters : sequence of callables
        Transformations that the worker will apply on the pixel weights prior
        to bincount. The required signature is ``f(chunk, data_weights)``,
        where ``chunk`` is a dictionary containing the data chunk and bin
        information, and ``data_weights`` is the 1D-array of current pixel
        weights of the chunk. The filter must return the updated array of pixel
        weights.

    Example
    -------
    # Create a worker and register filters f1, f2, f3
    >>> w = Worker('mycooler.coo', [f1, f2, f3])

    # Map to three consecutive chunks
    >>> marg1 = w((0, 1000))
    >>> marg2 = w((1000, 2000))
    >>> marg3 = w((2000, 3000))

    # Reduce the partial marginals into one
    >>> marg = np.sum([marg1, marg2, marg3], axis=0)

    """
    def __init__(self, cooler_path, filters=None):
        self.filepath = cooler_path
        self.filters = filters if filters is not None else []

    def __call__(self, span):
        lo, hi = span
        import numpy as np
        import h5py

        # prepare chunk dict
        chunk = {}
        lock.acquire()
        with h5py.File(self.filepath, 'r') as h5:
            n_bins_total = h5['bins/chrom_id'].shape[0]
            chunk['bin1_id'] = h5['matrix/bin1_id'][lo:hi]
            chunk['bin2_id'] = h5['matrix/bin2_id'][lo:hi]
            chunk['count'] = h5['matrix/count'][lo:hi]
            chunk['bintable'] = {
                'chrom_id': h5['bins/chrom_id'][:],
                'start':    h5['bins/start'][:],
                'end':      h5['bins/end'][:],
            }
        lock.release()

        # apply filters to chunk
        data_weights = chunk['count']
        for filter_ in self.filters:
            data_weights = filter_(chunk, data_weights)

        # marginalize
        marg = np.zeros(n_bins_total, dtype=float)
        marg += np.bincount(
            chunk['bin1_id'], weights=data_weights, minlength=n_bins_total)
        marg += np.bincount(
            chunk['bin2_id'], weights=data_weights, minlength=n_bins_total)
        return marg


class BinarizeFilter(object):
    def __call__(self, chunk, data_weights):
        data_weights[data_weights != 0] = 1
        return data_weights


class DropDiagFilter(object):
    def __init__(self, n_diags):
        self.n_diags = n_diags
    def __call__(self, chunk, data_weights):
        mask = np.abs(chunk['bin1_id'] - chunk['bin2_id']) < self.n_diags
        data_weights[mask] = 0
        return data_weights


class CisOnlyFilter(object):
    def __call__(self, chunk, data_weights):
        chrom_ids = chunk['bintable']['chrom_id']
        mask = chrom_ids[chunk['bin1_id']] != chrom_ids[chunk['bin2_id']]
        data_weights[mask] = 0
        return data_weights


class TimesOuterProductFilter(object):
    def __init__(self, vec):
        self.vec = vec
    def __call__(self, chunk, data_weights):
        data_weights = (self.vec[chunk['bin1_id']]
                            * self.vec[chunk['bin2_id']]
                            * data_weights)
        return data_weights


def iterative_correction(coo, chunksize=None, map=map, tol=1e-5,
                         min_nnz=0, min_count=0,
                         cis_only=False, ignore_diags=False):
    """
    Iterative correction or matrix balancing of a sparse Hi-C contact map in
    Cooler HDF5 format.

    Parameters
    ----------
    coo : h5py.File object
        Cooler file
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
    cis_only: bool, optional
        Do iterative correction on intra-chromosomal data only.
        Inter-chromosomal data is ignored.
    ignore_diags : int or False, optional
        Drop elements occurring on the first ``ignore_diags`` diagonals of the
        matrix.


    Returns
    -------
    bias : 1D array, whose shape is the number of bins in ``coo``.
        Vector of bin bias weights to normalize the observed contact map.
        Dropped bins will be assigned the value NaN.
        N[i, j] = O[i, j] * bias[i] * bias[j]

    """

    # Divide the number of elements into non-overlapping chunks
    nnz = coo.attrs['nnz']
    if chunksize is None:
        spans = [(0, nnz)]
    else:
        edges = np.arange(0, nnz+chunksize, chunksize)
        spans = list(zip(edges[:-1], edges[1:]))

    # List of pre-marginalization data transformations
    base_filters = []
    if cis_only:
        base_filters.append(CisOnlyFilter())
    if ignore_diags:
        base_filters.append(DropDiagFilter(ignore_diags))

    # Initialize the bias weights
    n_bins = coo.attrs['nbins']
    bias = np.ones(n_bins, dtype=float)

    # Drop bins with too few nonzeros from bias
    filters = [BinarizeFilter()] + base_filters
    marg_partials = map(Worker(coo.filename, filters), spans)
    marg_nnz = np.sum(list(marg_partials), axis=0)
    bias[marg_nnz < min_nnz] = 0

    # Drop bins with too few total counts from bias
    filters = base_filters
    marg_partials = map(Worker(coo.filename, filters), spans)
    marg = np.sum(list(marg_partials), axis=0)
    bias[marg < min_count] = 0

    # Do balancing
    while True:
        filters = base_filters + [TimesOuterProductFilter(bias)]
        worker = Worker(coo.filename, filters)
        marg_partials = map(worker, spans)
        marg = np.sum(list(marg_partials), axis=0)

        marg_ = marg[marg != 0]
        m = marg / marg_.mean()
        m[m == 0] = 1
        bias /= m

        var = marg_.var()
        print("variance is", var)
        if var < tol:
            bias[bias==0] = np.nan
            break

    return bias
