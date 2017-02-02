from __future__ import print_function

import cooler.contrib.higlass as cch
import h5py
import os.path as op

testdir = op.realpath(op.dirname(__file__))

def test_data_retrieval():
    data_file = op.join(testdir, 'data', 'dixon2012-h1hesc-hindiii-allreps-filtered.1000kb.multires.cool')
    
    f = h5py.File(data_file, 'r')
    
    data = cch.get_data(f, 0, 0, 3276799999, 0, 3276799999)

    assert(data['genome_start1'].iloc[0] == 0.)
    assert(data['genome_start2'].iloc[0] == 0.)
