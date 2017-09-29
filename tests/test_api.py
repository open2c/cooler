# -*- coding: utf-8 -*-
from __future__ import division, print_function
import os.path as op

from scipy import sparse
import numpy as np
import pandas
import h5py

from nose.tools import assert_raises
import cooler.api
import mock

testdir = op.realpath(op.dirname(__file__))


class MockHDF5(dict):
    file = mock.Mock(['mode'])

binsize = 100
n_bins = 20
r = sparse.random(n_bins, n_bins, density=1, random_state=1)
r = sparse.triu(r, k=1).tocsr()
r_full = r.toarray() + r.toarray().T

mock_cooler = MockHDF5({
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

mock_cooler.file = mock_cooler
mock_cooler.file.mode = 'r'
mock_cooler.file.filename = 'mock.cool'
mock_cooler.name = '/'
mock_cooler['/'] = mock_cooler

chromID_lookup = pandas.Series({'chr1': 0, 'chr2': 1})


def test_get():
    table = cooler.api.get(mock_cooler['chroms'])
    assert np.all(table['length'] == mock_cooler['chroms']['length'])


def test_chromtable():
    table = cooler.api.chroms(mock_cooler)
    assert np.all(table['length'] == mock_cooler['chroms']['length'])


def test_bintable():
    lo, hi = 2, 10
    table = cooler.api.bins(mock_cooler, lo, hi)
    assert np.all(chromID_lookup[table['chrom']] == mock_cooler['bins']['chrom'][lo:hi])
    assert np.all(table['start'] == mock_cooler['bins']['start'][lo:hi])
    assert np.all(table['end'] == mock_cooler['bins']['end'][lo:hi])


def test_pixeltable():
    lo, hi = 2, 10
    table = cooler.api.pixels(mock_cooler, lo, hi, join=False)
    assert np.all(table['bin1_id'] == mock_cooler['pixels']['bin1_id'][lo:hi])
    assert np.all(table['bin2_id'] == mock_cooler['pixels']['bin2_id'][lo:hi])

    table = cooler.api.pixels(mock_cooler, lo, hi, join=True)
    assert table.shape == (hi-lo, len(mock_cooler['pixels']) + 4)


def test_info():
    pass


def test_cooler():
    c = cooler.Cooler(mock_cooler)

    # bin table
    table = c.bins().fetch('chr1')
    assert np.all(table['start'] == mock_cooler['bins']['start'][0:10])
    assert np.all(table['end'] == mock_cooler['bins']['end'][0:10])

    # offsets
    assert c.offset('chr1') == 0
    assert c.extent('chr1') == (0, 10)

    # 2D range queries as rectangular or triangular
    A1 = np.triu(c.matrix(balance=False).fetch('chr2'))
    df = c.matrix(as_pixels=True, join=False, balance=False).fetch('chr2')
    i0 = c.offset('chr2')
    i, j, v = df['bin1_id'], df['bin2_id'], df['count']
    mat = sparse.coo_matrix((v, (i-i0, j-i0)), (A1.shape))
    A2 = np.triu(mat.toarray())
    assert np.all(A1 == A2)


def test_annotate():
    c = cooler.Cooler(mock_cooler)

    # works with full bin table / view or only required bins
    df = c.matrix(as_pixels=True, balance=False).fetch('chr1')
    df1 = cooler.annotate(df, c.bins()[:])
    df2 = cooler.annotate(df, c.bins())
    df3 = cooler.annotate(df, c.bins().fetch('chr1'))
    assert np.all(df1 == df2)
    assert np.all(df1 == df3)

    # works on empty dataframe
    df4 = cooler.annotate(df[0:0], c.bins()[:])
    assert np.all(df4.columns == df3.columns)
    assert len(df4) == 0


def test_matrix_as_pixels():
    c = cooler.Cooler(op.join(
        testdir,
        'data',
        'dixon2012-h1hesc-hindiii-allreps-filtered.1000kb.multires.cool::4'))
    df = c.matrix(as_pixels=True, join=True, balance=True).fetch(
        "chr10:6052652-6104288", 
        "chr10:7052652-8104288")
    print(df)
    assert len(df.dropna()) == 2
