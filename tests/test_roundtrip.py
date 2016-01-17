# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals
from six import BytesIO
import os

import numpy as np
import h5py

from nose.tools import assert_raises
import cooler

testdir = os.path.dirname(os.path.realpath(__file__))
testfile_path = os.path.join(testdir, 'test.cool')


def teardown_func():
    try:
        os.remove(testfile_path)
    except OSError:
        pass


def test_roundtrip():
    chromtable = cooler.read_chrominfo(
        'https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes',
        name_patterns=(r'^chr[0-9]+$', r'chrX$'))

    binsize = 2000000
    bintable = cooler.make_bintable(chromtable['length'], binsize)

    heatmap = np.load(os.path.join(testdir, 'IMR90_inSitu-all-MboI-2000k.npy'))
    with h5py.File(testfile_path, 'w') as h5:
        cooler.io.from_dense(h5, chromtable, bintable, heatmap,
                             binsize=binsize,
                             info={'genome-assembly': 'hg19'})

    h5 = h5py.File(testfile_path, 'r')
    new_chromtable = cooler.chromtable(h5)
    assert np.all(chromtable['name'] == new_chromtable['name'])

    new_bintable = cooler.bintable(h5)
    assert np.all(bintable == new_bintable)

    info = cooler.info(h5)
    assert info['genome-assembly'] == 'hg19'
    assert info['bin-type'] == 'fixed'
    assert info['bin-size'] == binsize

    mat = cooler.matrix(h5, 0, 100, 0, 100, 'count')
    assert mat.shape == (100, 100)
    assert np.allclose(heatmap[:100,:100], mat.toarray())

    mat = cooler.Cooler(h5).matrix('count')[:100, :100]
    assert mat.shape == (100, 100)
    assert np.allclose(heatmap[:100,:100], mat.toarray())

test_roundtrip.teardown = teardown_func
