# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals
from scipy import sparse
import numpy as np
import pandas
import h5py

from nose.tools import assert_raises
import cooler.api


class MockCooler(dict):
    pass

binsize = 100
n_bins = 20
r = sparse.random(n_bins, n_bins, density=1, random_state=1)
r = sparse.triu(r, k=1).tocsr()
r_full = r.toarray() + r.toarray().T

mock_cooler = MockCooler({
    'scaffolds': {
        'name':   np.array(['chr1', 'chr2'], dtype='S'),
        'length': np.array([1000, 1000], dtype=np.int32),
    },
    'bins': {
        'chrom_id': np.array([0,0,0,0,0,0,0,0,0,0,
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
    'matrix': {
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


def test_get():
    table = cooler.api.get(mock_cooler, 'scaffolds')
    assert np.all(table['length'] == mock_cooler['scaffolds']['length'])


def test_chromtable():
    table = cooler.api.chromtable(mock_cooler)
    assert np.all(table['length'] == mock_cooler['scaffolds']['length'])


def test_bintable():
    lo, hi = 2, 10
    table = cooler.api.bintable(mock_cooler, lo, hi)
    assert np.all(chromID_lookup[table['chrom']] == mock_cooler['bins']['chrom_id'][lo:hi])
    assert np.all(table['start'] == mock_cooler['bins']['start'][lo:hi])
    assert np.all(table['end'] == mock_cooler['bins']['end'][lo:hi])


def test_pixeltable():
    lo, hi = 2, 10
    table = cooler.api.pixeltable(mock_cooler, lo, hi, join=False)
    assert np.all(table['bin1_id'] == mock_cooler['matrix']['bin1_id'][lo:hi])
    assert np.all(table['bin2_id'] == mock_cooler['matrix']['bin2_id'][lo:hi])


def test_info():
    pass


def test_cooler():
    pass
