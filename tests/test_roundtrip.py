# -*- coding: utf-8 -*-
from __future__ import division, print_function
from six import BytesIO, iteritems
import tempfile
import os

import numpy as np
import h5py

from nose.tools import assert_raises
from nose import with_setup
import cooler

testdir = os.path.dirname(os.path.realpath(__file__))
tmp = tempfile.gettempdir()
testfile_path = os.path.join(tmp, 'test.cool')


def teardown_func():
    try:
        os.remove(testfile_path)
    except OSError:
        pass


@with_setup(teardown=teardown_func)
def test_roundtrip():
    chromsizes = cooler.read_chromsizes(
        'https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes',
        name_patterns=(r'^chr[0-9]+$', r'chrX$'))
    binsize = 2000000
    bintable = cooler.binnify(chromsizes, binsize)

    heatmap = np.load(os.path.join(testdir, 'data', 'IMR90-MboI-matrix.2000kb.npy'))
    reader = cooler.io.ArrayLoader(bintable, heatmap, 100000)
    cooler.io.create(testfile_path, bintable, reader, assembly='hg19')

    h5 = h5py.File(testfile_path, 'r')
    new_chromtable = cooler.chroms(h5)
    assert np.all(chromsizes.index == new_chromtable['name'])

    new_bintable = cooler.bins(h5)
    assert np.all(bintable == new_bintable)

    info = cooler.info(h5)
    assert info['genome-assembly'] == 'hg19'
    assert info['bin-type'] == 'fixed'
    assert info['bin-size'] == binsize

    mat = cooler.matrix(h5, 0, 100, 0, 100, 'count', balance=False)
    assert mat.shape == (100, 100)
    assert np.allclose(heatmap[:100,:100], mat)

    mat = cooler.Cooler(h5).matrix('count', balance=False)[:100, :100]
    assert mat.shape == (100, 100)
    assert np.allclose(heatmap[:100,:100], mat)

    mat = cooler.matrix(h5, 100, 200, 100, 200, 'count', balance=False)
    assert mat.shape == (100, 100)
    assert np.allclose(heatmap[100:200,100:200], mat)

    mat = cooler.Cooler(h5).matrix('count', balance=False)[100:200, 100:200]
    assert mat.shape == (100, 100)
    assert np.allclose(heatmap[100:200,100:200], mat)
