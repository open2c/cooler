# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals
from six import BytesIO
import os

import numpy as np
import requests
import h5py

import cooler


def teardown_func():
    try:
        os.remove('test.cool')
    except OSError:
        pass


def test_roundtrip():
    # *.chrom.sizes or *.chromInfo.txt files come from UCSC GoldenPath
    r = requests.get('https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes')
    chrom_table = cooler.genome.read_chromsizes(BytesIO(r.content))
    # exclude Y and M
    chrom_table = chrom_table.loc['chr1':'chrX']

    binsize = 2000000
    bin_table = cooler.genome.binnify(chrom_table, binsize)

    store = np.load('IMR90_inSitu-all-MboI-2000k.npy')
    with h5py.File('test.cool', 'w') as h5:
        cooler.from_dense(h5, chrom_table, bin_table, store['heatmap'][:],
                          bintype='fixed',
                          metadata={'genome-assembly': 'hg19'},
                          h5opts={})

    clr = h5py.File('test.cool', 'r')
    new_chrom_table = cooler.get_scaffolds(clr)
    assert np.all(chrom_table.index == new_chrom_table.index)

    new_bin_table = cooler.get_bins(clr)
    assert np.all(new_bin_table == bin_table)

    info = cooler.get_metadata(clr)
    assert info['genome-assembly'] == 'hg19'
    assert info['bin-type'] == 'fixed'
    assert info['bin-size'] == binsize

    mat = cooler.get_matrix(clr, ('chr1', 0, 200000000), dense=False)
    assert mat.shape == (1000, 1000)

test_roundtrip.teardown = teardown_func

