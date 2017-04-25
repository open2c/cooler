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


def test_balancing():
    fp = os.path.join(testdir, 'data', 'GM12878-MboI-matrix.2000kb.cool')
    tol = 1e-2

    weights, stats = cooler.ice.iterative_correction(
        cooler.Cooler(fp), 
        ignore_diags=1, 
        min_nnz=10, 
        tol=tol)

    # Extract matrix and apply weights
    mat = cooler.Cooler(fp).matrix(balance=False, sparse=True)[:, :]
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

    # Check that the balanced marginal is unity
    conv_marg = marg[~np.isnan(marg)].mean()
    err_marg = marg[~np.isnan(marg)].std()
    assert np.isclose(conv_marg, 1, atol=err_marg)


def test_balancing_cisonly():
    fp = os.path.join(testdir, 'data', 'GM12878-MboI-matrix.2000kb.cool')
    tol = 1e-2

    with h5py.File(fp, 'r') as h5:
        chrom_offsets = h5['indexes/chrom_offset'][:]
        weights, stats = cooler.ice.iterative_correction(
            cooler.Cooler(h5), 
            ignore_diags=1, 
            min_nnz=10, 
            tol=tol, 
            cis_only=True)

    # Extract matrix and apply weights
    mat = cooler.Cooler(fp).matrix(balance=False, sparse=True)[:, :]
    mat.data = weights[mat.row] * weights[mat.col] * mat.data
    arr = mat.toarray()

    # Re-apply bin level filters
    mask = np.isnan(weights)
    arr[mask, :] = 0
    arr[:, mask] = 0
    
    # Apply diagonal filter
    np.fill_diagonal(arr, 0)

    # Filter out trans data
    spans = list(zip(chrom_offsets[:-1], chrom_offsets[1:]))
    from scipy.linalg import block_diag
    blocks = [np.ones((hi -lo,) * 2) for lo, hi in spans]
    mask = block_diag(*blocks).astype(bool)
    arr[~mask] = 0

    # Check that the balanced marginal is flat
    marg = np.sum(arr, axis=0)
    for lo, hi in spans:
        m = marg[lo:hi]
        m = m[m != 0]
        if len(m):
            print(lo, hi)
            var = np.var(m[m != 0])
            assert var < tol 

            # Check that the balanced marginal is unity
            conv_marg = m[~np.isnan(m)].mean()
            err_marg = m[~np.isnan(m)].std()
            assert np.isclose(conv_marg, 1, atol=err_marg)
