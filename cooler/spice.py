from __future__ import division, print_function, unicode_literals
from contextlib import contextmanager
from multiprocessing import Pool
from copy import copy
import warnings

from scipy.sparse import coo_matrix
import numpy as np
import numexpr
import pandas
import h5py
import six

from hiclib.hicShared import h5dictBinarySearch, binarySearch
from mirnylib.genome import Genome


def _marginalize(lo, hi, filepath, bintable, weights, filters=None, nonzero=False):
    """
    Actual worker of coverage calculation

    NOTE: Change this to take in a number of arbitrary filters.
    Also, chrom_ids should be in the file's bin table

    e.g.
    nonzeroCounts --> filter the weights to be 1
    cisOnly --> filter weights where bins i and j belong to different chroms
    skipDiag --> filter weights where bins i and j are too close

    params: filename, filters, lo, hi

    """
    import numpy as np
    import h5py

    print(lo, hi)

    with h5py.File(filename, 'r') as h5:
        bin1   = h5["heatmap/bin1"][lo:hi]
        bin2   = h5["heatmap/bin2"][lo:hi]
        counts = h5["heatmap/count"][lo:hi]

    if nonzero:
        weights = np.ones(len(bin1))
    else:
        weights = weights[bin1] * weights[bin2] * counts

    if filters is not None:
        for filter in filters:
            weights = filter(weights, bin1, bin2, counts, bintable)

    marg = np.zeros(len(bintable), dtype=np.float64)
    marg += np.bincount(bin1, weights=weights, minlength=n_bins)
    marg += np.bincount(bin2, weights=weights, minlength=n_bins)
    return marg


def marginalize(filepath, bintable, weights=None, filters=None, nonzero=False, 
                chunksize=100000000, parallel=False):
    if weights is not None:
        weights = np.array(weights, dtype=np.float64)
    else:
        weights = np.ones(len(bintable), dtype=np.float64)

    with open_hdf5(filepath, 'r') as h5:
        nnz = len(h5['bin1'])

    inputs = []
    for lo, hi in _get_random_chunks(nnz, chunksize):
        inputs.append(
            (lo, hi, filepath, bintable, weights, filters, nonzero))

    if parallel:
        map = multiprocessing.Pool().map

    parts = map(_marginalize, inputs)
    marg = np.sum(parts, axis=0)
    return marg


def cisonly_filter(weights, bin1, bin2, counts, bintable):
    chrom_ids = bintable['chrom_id']
    weights[chrom_ids[bin1] != chrom_ids[bin2]] = 0
    return weights

def skipdiag_filter(weights, bin1, bin2, counts, binsize):
    weights[np.abs(bin1 - bin2) < skipDiag] = 0
    return weights




def iterative_correction(hmH5, genome_table, 
                         min_nonzero=20, min_count=40,
                         cisOnly=False, skipDiag=2, 
                         tol=1e-6, chunksize=40000000):
    n_bins = genome_table['chrom_id']
    genome_table['weights'] = np.ones(n_bins, dtype=np.float64)

    # Apply nonzero filter
    nz_marg = nonzero_marginals(hmH5, genome_table)
    genome_table['weights'][nz_marg < min_nonzero] = 0
    setToZero = ((nz_marg < min_nonzero) & (nz_marg > 0)).sum()
    remaining = (nz_marg >= min_nonzero).sum()
    print(("Set to zero {0} bins with less than {1} nonzero interactions;"
           "{2} elements remaining").format(setToZero, min_nonzero, remaining))

    # Apply minimum count filter
    marg = marginals(
        hmH5, genome_table, 
        weights=True, cisOnly=cisOnly, skipDiag=skipDiag, chunksize=chunksize)
    genome_table['weights'][marg < min_count] = 0
    setToZero = ((marg < minCount) * (marg > 0)).sum()
    remaining = (marg >= minCount).sum()
    print(("Set to zero {0} bins with less than {1} reads;"
           "{2} elements remaining").format(setToZero, min_count, remaining))

    # Do correction
    marg /= np.mean(marg[marg > 0])
    while True:
        marg[genome_table['weights'] == 0] = 1
        genome_table["weights"] /= marg
        marg = marginals(
            hmH5, genome_table, 
            weights=True, cisOnly=cisOnly, 
            skipDiag=skipDiag, chunksize=chunksize)
        marg /= np.mean(marg[marg > 0])
        print("variance is", marg[marg > 0].var())

        if marg[marg > 0].var() < tol:
            marg[genome_table['weights'] == 0] = 1
            genome_table['weights'] /= marg
            break
    return genome_table
