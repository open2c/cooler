# -*- coding: utf-8 -*-
from __future__ import division, print_function
from collections import OrderedDict
from six import iteritems
import mock
import os

from scipy import sparse
import numpy as np
import pandas
import h5py

from nose.tools import assert_raises
from nose import with_setup
import cooler

testdir = os.path.dirname(os.path.realpath(__file__))


# class MockCooler(dict):
#     file = mock.Mock(['mode'])


# def generate_mock_cooler(mat):
#     binsize = 100
#     n_bins = mat.shape[0]
#     r = sparse.triu(mat, k=0).tocsr()  # drop lower triangle, keep main diagonal
#     n_chr1 = n_bins //2
#     n_chr2 = n_bins - (n_bins//2)

#     mock_cooler = MockCooler({
#         'chroms': {
#             'name':   np.array(['chr1', 'chr2'], dtype='S'),
#             'length': np.array([1000, 1000], dtype=np.int32),
#         },
#         'bins': {
#             'chrom':    np.array(([0] * n_chr1) + ([1] * n_chr2), dtype=int),
#             'start':    np.r_[np.arange(0, n_chr1*binsize, binsize),
#                               np.arange(0, n_chr2*binsize, binsize)].astype(int),
#             'end':      np.r_[np.arange(binsize, n_chr1*binsize + binsize, binsize),
#                               np.arange(binsize, n_chr2*binsize + binsize, binsize)].astype(int),
#         },
#         'pixels': {
#             'bin1_id':  r.tocoo().row,
#             'bin2_id':  r.indices,
#             'count':    r.data,
#             'mask':     np.ones(r.nnz, dtype=bool),
#         },
#         'indexes': {
#             'chrom_offset': np.array([0, n_chr1, n_bins], dtype=np.int32),
#             'bin1_offset':  r.indptr,
#         },
#     })

#     mock_cooler.attrs = {
#         'bin-size': binsize,
#         'bin-type': 'fixed',
#         'nchroms': 2,
#         'nbins': n_bins,
#         'nnz': r.nnz,
#         'metadata': '{}',
#     }

#     mock_cooler.file.mode = 'r'

#     return mock_cooler


def test_balancing():
    fp = os.path.join(testdir, 'data', 'GM12878-MboI-matrix.2000kb.cool')
    tol = 1e-2

    with h5py.File(fp, 'r') as h5:
        weights = cooler.ice.iterative_correction(h5, ignore_diags=1, min_nnz=10, tol=tol)

    # Extract matrix and apply weights
    mat = cooler.Cooler(fp).matrix()[:, :]
    mat.data = weights[mat.row] * weights[mat.col] * mat.data
    arr = mat.toarray()

    # Re-apply bin level filters
    mask = np.isnan(weights)
    arr[mask, :] = 0
    arr[:, mask] = 0
    
    # Apply diagonal filter
    np.fill_diagonal(arr, 0)

    # Check that the balanced marginal is flat
    marg = np.sum(arr, axis=0)
    var = np.var(marg[marg != 0])
    assert var < tol 



