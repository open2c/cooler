# -*- coding: utf-8 -*-
from __future__ import division, print_function
from scipy import sparse
import numpy as np
import pandas
import h5py

from nose.tools import assert_raises
import cooler.core


class MockCooler(dict):
    pass

binsize = 100
n_bins = 20
r = sparse.random(n_bins, n_bins, density=1, random_state=1)
r = sparse.triu(r, k=1).tocsr()
r_full = r.toarray() + r.toarray().T

mock_cooler = MockCooler({
    'chroms': {
        'name':   np.array(['chr1', 'chr2'], dtype='S'),
        'length': np.array([1000, 1000], dtype=np.int32),
    },
    'bins': {
        'chrom':    np.array([0,0,0,0,0,0,0,0,0,0,
                              1,1,1,1,1,1,1,1,1,1], dtype=int),
        'start':    np.array([0,100,200,300,400,500,600,700,800,900,
                              0,100,200,300,400,500,600,700,800,900],
                              dtype=int),
        'end':      np.array([100,200,300,400,500,600,700,800,900,1000,
                              100,200,300,400,500,600,700,800,900,1000],
                              dtype=int),
        'mask':     np.array([1,1,1,1,1,1,1,1,1,1,
                              1,1,1,1,1,1,1,1,1,1], dtype=bool),
        'bias':     np.array([1,1,1,1,1,1,1,1,1,1,
                              1,1,1,1,1,1,1,1,1,1], dtype=float),
        'E1':       np.zeros(20, dtype=float),
    },
    'pixels': {
        'bin1_id':  r.tocoo().row,
        'bin2_id':  r.indices,
        'count':    r.data,
        'mask':     np.ones(r.nnz, dtype=bool),
    },
    'indexes': {
        'chrom_offset': np.array([0, 10, 20], dtype=np.int32),  # nchroms + 1
        'bin1_offset':  r.indptr,  # nbins + 1
    },
})

mock_cooler.attrs = {
    'bin-size': binsize,
    'bin-type': 'fixed',
    'nchroms': 2,
    'nbins': n_bins,
    'nnz': r.nnz,
    'metadata': '{}',
}

chromID_lookup = pandas.Series({'chr1': 0, 'chr2': 1})


def test_selector1d():
    slicer = lambda fields, lo, hi: (lo, hi)
    fetcher = lambda x: x
    nmax = 50

    s = cooler.core.RangeSelector1D(None, slicer, fetcher, nmax)
    assert s[30] == (30, 31)
    assert s[10:20] == (10, 20)
    assert s[:20] == (0, 20)
    assert s[10:] == (10, nmax)
    assert s[:] == (0, nmax)
    assert s[:nmax] == (0, nmax)
    assert s[:-10] == (0, nmax-10)
    assert s[1:1] == (1, 1)
    assert_raises(IndexError, lambda : s[:, :])
    assert_raises(ValueError, lambda : s[::2])
    #assert_raises(TypeError, lambda : s['blah'])
    assert s.shape == (nmax,)

    # FIXME - questionable behavior
    assert s[30:20] == (30, 20)  # lo > hi
    assert s[nmax+10:nmax+30] == (nmax+10, nmax+30)  # lo > nmax
    assert s[10.0] == (10, 11)  # accepting floats
    #assert s[10.1] == (10.1, 11.1)  # not casting
    #assert s[nmax+10] == (nmax+10, nmax+11)


def test_selector2d():
    slicer = lambda field, i0, i1, j0, j1: (i0, i1, j0, j1)
    fetcher = lambda x: x
    nmax = 50

    s = cooler.core.RangeSelector2D(None, slicer, fetcher, (nmax, nmax))
    assert s[30] == (30, 31, 0, nmax)
    assert s[10:20, 10:20] == (10, 20, 10, 20)
    assert s[:] == (0, nmax, 0, nmax)
    assert_raises(IndexError, lambda : s[:, :, :])
    assert_raises(ValueError, lambda : s[::2, :])
    assert s.shape == (nmax, nmax)


def test_region_to_extent():
    region = ('chr1', 159, 402)
    first, last = 1, 4
    assert cooler.api.region_to_extent(
        mock_cooler, chromID_lookup, region, binsize) == (first, last+1)
    assert cooler.api.region_to_extent(
        mock_cooler, chromID_lookup, region, None) == (first, last+1)

    region = ('chr1', 159, 400)
    first, last = 1, 3
    assert cooler.api.region_to_extent(
        mock_cooler, chromID_lookup, region, binsize) == (first, last+1)
    assert cooler.api.region_to_extent(
        mock_cooler, chromID_lookup, region, None) == (first, last+1)


def test_slice_matrix():
    slices = [
        (0, 10, 0, 10),
        (0, 10, 10, 20),
        (5, 15, 10, 20),
        (10, 20, 5, 15),
        (1, 1, 5, 15),
        (1, 1, 1, 1),
    ]
    for i0, i1, j0, j1 in slices:
        triu_reader = cooler.core.TriuReader(mock_cooler, 'count', max_chunk=10)

        # triangular query
        index = triu_reader.index_col(i0, i1, j0, j1)
        i, j, v = triu_reader.query(i0, i1, j0, j1)
        assert len(index) == len(i)

        # rectangular query
        i, j, v = cooler.core.query_rect(triu_reader.query, i0, i1, j0, j1)
        mat = sparse.coo_matrix((v, (i-i0, j-j0)), (i1-i0, j1-j0)).toarray()
        assert np.allclose(r_full[i0:i1, j0:j1], mat)


