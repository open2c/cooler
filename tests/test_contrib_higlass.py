from __future__ import print_function

import cooler.contrib.higlass as cch
import h5py
import os.path as op

testdir = op.realpath(op.dirname(__file__))


def test_data_retrieval():
    data_file = op.join(testdir, 'data', 'dixon2012-h1hesc-hindiii-allreps-filtered.1000kb.multires.cool')
    
    f = h5py.File(data_file, 'r')
    
    data = cch.get_data(f['0'], 0, 3276799999, 0, 3276799999)

    assert(data['genome_start1'].iloc[0] == 0.)
    assert(data['genome_start2'].iloc[0] == 0.)

    data = cch.get_data(f['4'], 0, 256000000, 0, 256000000)

    assert(data['genome_start1'].iloc[-1] > 255000000)
    assert(data['genome_start1'].iloc[-1] < 256000000)
    #print("ge1", data['genome_end1'])

    data = cch.get_data(f['4'], 0, 256000000, 0, 256000000, transform='weight')
    assert 'balanced' in data
    data = cch.get_data(f['4'], 0, 256000000, 0, 256000000, transform='none')
    assert 'count' in data


def test_get_info():
    data_file = op.join(testdir, 'data', 'dixon2012-h1hesc-hindiii-allreps-filtered.1000kb.multires.cool')

    info = cch.get_info(data_file)
    assert('transforms' in info)
