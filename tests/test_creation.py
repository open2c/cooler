# -*- coding: utf-8 -*-
from __future__ import division, print_function
from collections import OrderedDict
from six import iteritems
import tempfile
import os.path as op
import os
import numpy as np
import pandas as pd
import h5py

from _common import isolated_filesystem
import cooler.create
import cooler
import pytest


tmp = tempfile.gettempdir()
testdir = op.dirname(op.realpath(__file__))


@pytest.mark.parametrize("fp", [
    op.join(testdir, 'data', 'hg19.GM12878-MboI.matrix.2000kb.cool')
])
def test_create_append(fp):
    import dask.dataframe as dd

    c = cooler.Cooler(fp)
    chromsizes = c.chromsizes
    bins = c.bins()[:]
    pixels = c.pixels()[:]

    # create
    cooler.create.create(op.join(tmp, 'test.df.2000kb.cool'), bins, pixels)
    cooler.create.create(op.join(tmp, 'test.dict.2000kb.cool'), bins, {k:v for k,v in iteritems(pixels)})
    cooler.create.create(op.join(tmp, 'test.iter_df.2000kb.cool'), bins, [pixels])
    cooler.create.create(op.join(tmp, 'test.iter_dict.2000kb.cool'), bins, [{k:v for k,v in iteritems(pixels)}])
    ddf = dd.from_pandas(pixels, npartitions=3)
    cooler.create.create(op.join(tmp, 'test.ddf.2000kb.cool'), bins, ddf)

    # Append
    cooler.create.append(op.join(tmp, 'test.df.2000kb.cool'),
                     'bins',
                     {'start_1based': bins.apply(lambda x: x.start + 1, axis=1)})
    cooler.create.append(op.join(tmp, 'test.df.2000kb.cool'),
                     'bins',
                     {'ones': 1})
    series = (ddf['count'] / ddf['count'].sum())
    series.name = 'normed'
    cooler.create.append(op.join(tmp, 'test.df.2000kb.cool'),
                     'pixels',
                     series)
    cooler.create.append(op.join(tmp, 'test.df.2000kb.cool'),
                     'pixels',
                     series, force=True)
    cooler.create.append(op.join(tmp, 'test.df.2000kb.cool'),
                     'bins',
                     {'twos': [np.ones(1000, dtype=int)*2,
                               np.ones(561, dtype=int)*2]},
                     chunked=True, force=True)
    c2 = cooler.Cooler(op.join(tmp, 'test.df.2000kb.cool'))
    assert (len(c2.bins().columns) == 6)
    assert (len(c2.pixels().columns) == 4)


@pytest.mark.parametrize("f_hm,f_cool", [(
    op.join(testdir, 'data', 'hg19.IMR90-MboI.matrix.2000kb.npy'),
    op.join(tmp, 'test.cool'),
)])
def test_roundtrip(f_hm, f_cool):
    chromsizes = cooler.read_chromsizes(
        'https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes',
        name_patterns=(r'^chr[0-9]+$', r'chrX$'))
    binsize = 2000000
    bintable = cooler.binnify(chromsizes, binsize)

    heatmap = np.load(f_hm)
    reader = cooler.create.ArrayLoader(bintable, heatmap, 100000)
    cooler.create.create(f_cool, bintable, reader, assembly='hg19')

    h5 = h5py.File(f_cool, 'r')
    new_chromtable = cooler.api.chroms(h5)
    assert np.all(chromsizes.index == new_chromtable['name'])

    new_bintable = cooler.api.bins(h5)
    assert np.all(bintable == new_bintable)

    info = cooler.api.info(h5)
    assert info['genome-assembly'] == 'hg19'
    assert info['bin-type'] == 'fixed'
    assert info['bin-size'] == binsize

    mat = cooler.api.matrix(h5, 0, 100, 0, 100, 'count', balance=False)
    assert mat.shape == (100, 100)
    assert np.allclose(heatmap[:100,:100], mat)

    mat = cooler.Cooler(h5).matrix('count', balance=False)[:100, :100]
    assert mat.shape == (100, 100)
    assert np.allclose(heatmap[:100,:100], mat)

    mat = cooler.api.matrix(h5, 100, 200, 100, 200, 'count', balance=False)
    assert mat.shape == (100, 100)
    assert np.allclose(heatmap[100:200,100:200], mat)

    mat = cooler.Cooler(h5).matrix('count', balance=False)[100:200, 100:200]
    assert mat.shape == (100, 100)
    assert np.allclose(heatmap[100:200,100:200], mat)

    try:
        os.remove(f_cool)
    except OSError:
        pass


def test_rename_chroms():
    from shutil import copyfile

    with isolated_filesystem() as fs:
        copyfile(op.join(testdir, 'data', 'toy.asymm.4.cool'), 'toy.asymm.4.cool')
        clr = cooler.Cooler('toy.asymm.4.cool')
        assert clr.chromnames == ['chr1', 'chr2']
        cooler.rename_chroms(clr, {'chr1': '1', 'chr2': '2'})
        assert clr.chromnames == ['1', '2']  # the Cooler object is refreshed


def test_create_custom_cols():

    with isolated_filesystem() as fs:
        df = pd.DataFrame({
            'bin1_id': [0, 1, 1, 1, 2, 2, 3, 4, 5],
            'bin2_id': [1, 1, 3, 4, 5, 6, 7, 8, 9],
            'foo':     [1, 1, 1, 1, 1, 2, 2, 2, 2],
            'bar':     [.1, .2, .3, .4, .5, .6, .7, .8, .9],
        }, columns=['bin1_id', 'bin2_id', 'foo', 'bar'])
        bins = pd.DataFrame({
            'chrom': ['chr1']*5 + ['chr2']*5,
            'start': list(range(5))*2,
            'end': list(range(1,6))*2,
        })
        # works in unordered mode
        cooler.create_cooler('test.cool', bins, df, columns=['foo', 'bar'])
        clr = cooler.Cooler('test.cool')
        assert len(clr.pixels().columns) == 4
        assert np.allclose(df, clr.pixels()[['bin1_id', 'bin2_id', 'foo', 'bar']][:])

        # works in ordered mode
        cooler.create_cooler('test.cool', bins, df, columns=['foo', 'bar'], ordered=True)
        clr = cooler.Cooler('test.cool')
        assert len(clr.pixels().columns) == 4
        assert np.allclose(df, clr.pixels()[['bin1_id', 'bin2_id', 'foo', 'bar']][:])

        # raises if no custom columns specified and 'count' does not exist
        with pytest.raises(ValueError):
            cooler.create_cooler('test.cool', bins, df, columns=None, ordered=True)

