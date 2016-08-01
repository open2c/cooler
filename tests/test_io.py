# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals
import os

import numpy as np
import pandas
import h5py

from nose.tools import assert_raises
from nose import with_setup
import cooler.io
import cooler

testdir = os.path.dirname(os.path.realpath(__file__))
testfile_path = os.path.join(testdir, 'test.cool')


# Generate mock data
n_chroms = 2
clen = 2000
chromtable = pandas.DataFrame({
    'name': ['chr1', 'chr2'],
    'length': [clen, clen],
}, columns=['name', 'length'])
chromsizes = chromtable.set_index('name')['length']

class MockReads(dict):
    pass

n_records = 3000
np.random.seed(1)
chrms = np.random.randint(0, n_chroms, n_records * 2)
cuts = np.random.randint(0, clen, n_records * 2)
abs_cuts = np.array([clen * chrm + cut for chrm, cut in zip(chrms, cuts)])
abs_cuts1, abs_cuts2 = abs_cuts[:n_records], abs_cuts[n_records:]
mock_reads = MockReads({
    'chrms1': chrms[:n_records],
    'cuts1':  cuts[:n_records],
    'chrms2': chrms[n_records:],
    'cuts2':  cuts[n_records:],
})
# Triu-sort
mask = abs_cuts1 > abs_cuts2
mock_reads['chrms1'][mask], mock_reads['chrms2'][mask] = mock_reads['chrms2'][mask], mock_reads['chrms1'][mask]
mock_reads['cuts1'][mask], mock_reads['cuts2'][mask] = mock_reads['cuts2'][mask], mock_reads['cuts1'][mask]
abs_cuts1[mask], abs_cuts2[mask] = abs_cuts2[mask], abs_cuts1[mask]
idx = np.lexsort([abs_cuts2, abs_cuts1])
for key in mock_reads:
    mock_reads[key] = mock_reads[key][idx]


def teardown_func():
    try:
        os.remove(testfile_path)
    except OSError:
        pass


def should_not_depend_on_chunksize(bintable, binsize):
    # try different chunk sizes
    with h5py.File(testfile_path, 'w') as h5:
        reader = cooler.io.HDF5Aggregator(
            mock_reads, chromsizes, bintable,
            binsize, chunksize=66)
        cooler.io.create(h5, chromsizes, bintable, reader, binsize)
        oc1 = h5['indexes']['chrom_offset'][:]
        ob1 = h5['indexes']['bin1_offset'][:]
        p1 = cooler.pixeltable(h5, join=False)

    with h5py.File(testfile_path, 'w') as h5:
        reader = cooler.io.HDF5Aggregator(
            mock_reads, chromsizes, bintable,
            binsize, chunksize=666)
        cooler.io.create(h5, chromsizes, bintable, reader, binsize)
        oc2 = h5['indexes']['chrom_offset'][:]
        ob2 = h5['indexes']['bin1_offset'][:]
        p2 = cooler.pixeltable(h5, join=False)

    assert np.all(oc1 == oc2)
    assert np.all(ob1 == ob2)
    assert np.all(p1.values == p2.values)


def should_raise_if_input_not_sorted(bintable, binsize):
    # not sorted by chrm1
    with h5py.File(testfile_path, 'w') as h5:
        bad_reads = MockReads({
            'chrms1': mock_reads['chrms2'],
            'cuts1':  mock_reads['cuts2'],
            'chrms2': mock_reads['chrms1'],
            'cuts2':  mock_reads['cuts1'],
        })
        assert_raises(ValueError, cooler.io.HDF5Aggregator,
            bad_reads, chromsizes, bintable,
            binsize, chunksize=66)
    # not triu
    with h5py.File(testfile_path, 'w') as h5:
        bad_reads = MockReads({
            'chrms1': mock_reads['chrms1'].copy(),
            'cuts1':  mock_reads['cuts1'].copy(),
            'chrms2': mock_reads['chrms2'].copy(),
            'cuts2':  mock_reads['cuts2'].copy(),
        })
        bad_reads['chrms1'][0] = 0
        bad_reads['chrms2'][0] = 0
        bad_reads['cuts1'][0] = 10
        bad_reads['cuts2'][0] = 9
        reader = cooler.io.HDF5Aggregator(
            bad_reads, chromsizes, bintable,
            binsize, chunksize=66)
        assert_raises(ValueError, cooler.io.create,
            h5, chromsizes, bintable, reader, binsize)


def should_work_with_int32_cols(bintable, binsize):
    with h5py.File(testfile_path, 'w') as h5:
        reader = cooler.io.HDF5Aggregator(
            mock_reads, chromsizes, bintable,
            binsize, chunksize=66)
        cooler.io.create(h5, chromsizes, bintable, reader, binsize)
        oc1 = h5['indexes']['chrom_offset'][:]
        ob1 = h5['indexes']['bin1_offset'][:]
        p1 = cooler.pixeltable(h5, join=False)

    with h5py.File(testfile_path, 'w') as h5:
        mock_reads32 = MockReads({
            'chrms1': mock_reads['chrms1'].astype(np.int32),
            'cuts1':  mock_reads['cuts1'].astype(np.int32),
            'chrms2': mock_reads['chrms2'].astype(np.int32),
            'cuts2':  mock_reads['cuts2'].astype(np.int32),
        })
        reader = cooler.io.HDF5Aggregator(
            mock_reads32, chromsizes, bintable,
            binsize, chunksize=66)
        cooler.io.create(h5, chromsizes, bintable, reader, binsize)
        oc2 = h5['indexes']['chrom_offset'][:]
        ob2 = h5['indexes']['bin1_offset'][:]
        p2 = cooler.pixeltable(h5, join=False)

    assert np.all(oc1 == oc2)
    assert np.all(ob1 == ob2)
    assert np.all(p1.values == p2.values)


@with_setup(teardown=teardown_func)
def test_from_readhdf5():
    binsize = 100
    bintable = cooler.make_bintable(chromsizes, binsize)
    # assuming uniform bins
    yield should_not_depend_on_chunksize, bintable, binsize
    yield should_raise_if_input_not_sorted, bintable, binsize
    yield should_work_with_int32_cols, bintable, binsize
    # not assuming uniform bins
    yield should_not_depend_on_chunksize, bintable, None
    yield should_raise_if_input_not_sorted, bintable, None
    yield should_work_with_int32_cols, bintable, None
