from __future__ import absolute_import, print_function, division
import scipy.sparse as sps
import os.path as op
import pandas as pd
import numpy as np
import pandas
import h5py

import cooler.api

testdir = op.realpath(op.dirname(__file__))


def test_info():
    pass


def test_get(mock_cooler):
    table = cooler.api.get(mock_cooler['chroms'])
    assert np.all(table['length'] == mock_cooler['chroms']['length'])


def test_chromtable(mock_cooler):
    table = cooler.api.chroms(mock_cooler)
    assert np.all(table['length'] == mock_cooler['chroms']['length'])


def test_bintable(mock_cooler):
    chromID_lookup = pd.Series({'chr1': 0, 'chr2': 1})
    lo, hi = 2, 10
    table = cooler.api.bins(mock_cooler, lo, hi)
    assert np.all(chromID_lookup[table['chrom']] == mock_cooler['bins']['chrom'][lo:hi])
    assert np.all(table['start'] == mock_cooler['bins']['start'][lo:hi])
    assert np.all(table['end'] == mock_cooler['bins']['end'][lo:hi])


def test_bintable_many_contigs():
    # In a file with many contigs, bins/chrom does not have an ENUM header,
    # so chromosome names are taken from the chroms/name
    c = cooler.api.Cooler(op.join(testdir, 'data', 'manycontigs.1.cool'))
    bins = c.bins()[:10]
    assert pd.api.types.is_categorical_dtype(bins['chrom'].dtype)


def test_pixeltable(mock_cooler):
    lo, hi = 2, 10
    table = cooler.api.pixels(mock_cooler, lo, hi, join=False)
    assert np.all(table['bin1_id'] == mock_cooler['pixels']['bin1_id'][lo:hi])
    assert np.all(table['bin2_id'] == mock_cooler['pixels']['bin2_id'][lo:hi])

    table = cooler.api.pixels(mock_cooler, lo, hi, join=True)
    assert table.shape == (hi-lo, len(mock_cooler['pixels']) + 4)


def test_cooler(mock_cooler):
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
    mat = sps.coo_matrix((v, (i-i0, j-i0)), (A1.shape))
    A2 = np.triu(mat.toarray())
    assert np.all(A1 == A2)


def test_annotate(mock_cooler):
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


# def test_matrix_as_pixels():
#     c = cooler.Cooler(op.join(
#         testdir,
#         'data',
#         'dixon2012-h1hesc-hindiii-allreps-filtered.1000kb.multires.cool::4'))
#     df = c.matrix(as_pixels=True, join=True, balance=True).fetch(
#         "chr10:6052652-6104288",
#         "chr10:7052652-8104288")
#     print(df)
#     assert len(df.dropna()) == 2
