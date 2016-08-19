# -*- coding: utf-8 -*-
from __future__ import division, print_function
from multiprocessing import Pool, Lock
import warnings
import six

import numpy as np
import pandas
import h5py


# Lock prevents race condition if HDF5 file is already open before forking.
# See discussion: <https://groups.google.com/forum/#!topic/h5py/bJVtWdFtZQM>
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

    cooler_root : str
        HDF5 path to root group of a Cooler tree.

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
    def __init__(self, cooler_path, cooler_root, filters=None):
        self.filepath = cooler_path
        self.root = cooler_root
        self.filters = filters if filters is not None else []

    def __call__(self, span):
        lo, hi = span
        import numpy as np
        import h5py

        # prepare chunk dict
        chunk = {}
        lock.acquire()
        with h5py.File(self.filepath, 'r') as h5:
            coo = h5[self.root]
            n_bins_total = coo['bins/chrom'].shape[0]
            chunk['bin1_id'] = coo['pixels/bin1_id'][lo:hi]
            chunk['bin2_id'] = coo['pixels/bin2_id'][lo:hi]
            chunk['count']   = coo['pixels/count'][lo:hi]
            chunk['bintable'] = {
                'chrom_id':  coo['bins/chrom'][:],
                'start':     coo['bins/start'][:],
                'end':       coo['bins/end'][:],
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


def mad(data, axis=None):
    return np.median(np.abs(data - np.median(data, axis)), axis)


def iterative_correction(h5, cooler_root='/', chunksize=None, map=map, tol=1e-5,
                         min_nnz=0, min_count=0, mad_max=0,
                         cis_only=False, ignore_diags=False,
                         max_iters=200):
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

    Returns
    -------
    bias : 1D array, whose shape is the number of bins in ``h5``.
        Vector of bin bias weights to normalize the observed contact map.
        Dropped bins will be assigned the value NaN.
        N[i, j] = O[i, j] * bias[i] * bias[j]

    """
    filepath = h5.file.filename

    # Divide the number of elements into non-overlapping chunks
    nnz = h5[cooler_root].attrs['nnz']
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
    n_bins = h5[cooler_root].attrs['nbins']
    bias = np.ones(n_bins, dtype=float)

    # Drop bins with too few nonzeros from bias
    if min_nnz > 0:
        filters = [BinarizeFilter()] + base_filters
        marg_partials = map(Worker(filepath, cooler_root, filters), spans)
        marg_nnz = np.sum(list(marg_partials), axis=0)
        bias[marg_nnz < min_nnz] = 0

    filters = base_filters
    marg_partials = map(Worker(filepath, cooler_root, filters), spans)
    marg = np.sum(list(marg_partials), axis=0)

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
    for _ in range(max_iters):
        filters = base_filters + [TimesOuterProductFilter(bias)]
        worker = Worker(filepath, cooler_root, filters)
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
    else:
        warnings.warn('Iteration limit reached without convergence.')

    return bias
