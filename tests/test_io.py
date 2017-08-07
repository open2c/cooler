# -*- coding: utf-8 -*-
from __future__ import division, print_function
from collections import OrderedDict
from six import iteritems
import tempfile
import os.path as op
import os

import numpy as np
import pandas
import h5py

from nose.tools import assert_raises
from nose import with_setup
import cooler.io
import cooler

testdir = op.dirname(op.realpath(__file__))
tmp = tempfile.gettempdir()
testfile_path = op.join(tmp, 'test.cool')


# Generate mock data
class MockReads(dict):
    pass

n_chroms = 2
clen = 2000
chromsizes = pandas.Series(index=['chr1', 'chr2'], data=[clen, clen])
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


def _alternating_bins(chromsizes, steps):
    def _each(chrom):
        clen = chromsizes[chrom]
        edges = [0]
        total = 0
        i = 0
        while edges[i] < clen:
            step = steps[i % len(steps)]
            total += step
            i += 1
            edges.append(total)
            print(edges)
        edges[-1] = clen
        return pandas.DataFrame(
            {'chrom': chrom, 'start': edges[:-1], 'end': edges[1:]},
            columns=['chrom', 'start', 'end'])
    return pandas.concat(map(_each, chromsizes.keys()), axis=0, ignore_index=True)


def teardown_func():
    try:
        os.remove(testfile_path)
    except OSError:
        pass


def should_not_depend_on_chunksize(bintable):
    # try different chunk sizes
    binner = cooler.io.HDF5Aggregator(
        mock_reads, chromsizes, bintable, chunksize=66)
    cooler.io.create(testfile_path, bintable, binner)
    with h5py.File(testfile_path, 'r') as h5:
        oc1 = h5['indexes']['chrom_offset'][:]
        ob1 = h5['indexes']['bin1_offset'][:]
        p1 = cooler.pixels(h5, join=False)

    binner = cooler.io.HDF5Aggregator(
        mock_reads, chromsizes, bintable, chunksize=666)
    cooler.io.create(testfile_path, bintable, binner)
    with h5py.File(testfile_path, 'r') as h5:
        oc2 = h5['indexes']['chrom_offset'][:]
        ob2 = h5['indexes']['bin1_offset'][:]
        p2 = cooler.pixels(h5, join=False)

    assert np.all(oc1 == oc2)
    assert np.all(ob1 == ob2)
    assert np.all(p1.values == p2.values)


def should_raise_if_input_not_sorted(bintable):
    # not sorted by chrm1
    #with h5py.File(testfile_path, 'w') as h5:
    bad_reads = MockReads({
        'chrms1': mock_reads['chrms2'],
        'cuts1':  mock_reads['cuts2'],
        'chrms2': mock_reads['chrms1'],
        'cuts2':  mock_reads['cuts1'],
    })
    assert_raises(ValueError, cooler.io.HDF5Aggregator,
        bad_reads, chromsizes, bintable, chunksize=66)

    # not triu
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
    binner = cooler.io.HDF5Aggregator(
        bad_reads, chromsizes, bintable, chunksize=66)
    assert_raises(ValueError, cooler.io.create,
        testfile_path, bintable, binner)


def should_work_with_int32_cols(bintable):
    # int64
    binner = cooler.io.HDF5Aggregator(
        mock_reads, chromsizes, bintable, chunksize=66)
    cooler.io.create(testfile_path, bintable, binner)
    with h5py.File(testfile_path, 'r') as h5:
        oc1 = h5['indexes']['chrom_offset'][:]
        ob1 = h5['indexes']['bin1_offset'][:]
        p1 = cooler.pixels(h5, join=False)

    # int32
    mock_reads32 = MockReads({
        'chrms1': mock_reads['chrms1'].astype(np.int32),
        'cuts1':  mock_reads['cuts1'].astype(np.int32),
        'chrms2': mock_reads['chrms2'].astype(np.int32),
        'cuts2':  mock_reads['cuts2'].astype(np.int32),
    })
    binner = cooler.io.HDF5Aggregator(
        mock_reads32, chromsizes, bintable, chunksize=66)
    cooler.io.create(testfile_path, bintable, binner)
    with h5py.File(testfile_path, 'r') as h5:
        oc2 = h5['indexes']['chrom_offset'][:]
        ob2 = h5['indexes']['bin1_offset'][:]
        p2 = cooler.pixels(h5, join=False)

    assert np.all(oc1 == oc2)
    assert np.all(ob1 == ob2)
    assert np.all(p1.values == p2.values)


@with_setup(teardown=teardown_func)
def test_from_readhdf5():
    # uniform bins
    binsize = 100
    bintable = cooler.binnify(chromsizes, binsize)
    yield should_not_depend_on_chunksize, bintable
    yield should_raise_if_input_not_sorted, bintable
    yield should_work_with_int32_cols, bintable

    # non-uniform bins
    steps = [10, 100]
    bintable = _alternating_bins(chromsizes, steps)
    yield should_not_depend_on_chunksize, bintable
    yield should_raise_if_input_not_sorted, bintable
    yield should_work_with_int32_cols, bintable


def test_create_append():
    import dask.dataframe as dd

    c = cooler.Cooler(op.join(testdir, 'data', 'GM12878-MboI-matrix.2000kb.cool'))
    chromsizes = c.chromsizes
    bins = c.bins()[:]
    pixels = c.pixels()[:]

    cooler.io.create(op.join(tmp, 'test.df.2000kb.cool'), bins, pixels)
    cooler.io.create(op.join(tmp, 'test.dict.2000kb.cool'), bins, {k:v for k,v in iteritems(pixels)})
    cooler.io.create(op.join(tmp, 'test.iter_df.2000kb.cool'), bins, [pixels])
    cooler.io.create(op.join(tmp, 'test.iter_dict.2000kb.cool'), bins, [{k:v for k,v in iteritems(pixels)}])
    ddf = dd.from_pandas(pixels, npartitions=3)
    cooler.io.create(op.join(tmp, 'test.ddf.2000kb.cool'), bins, ddf)

    # Append
    cooler.io.append(op.join(tmp, 'test.df.2000kb.cool'), 
                     'bins', 
                     {'start_1based': bins.apply(lambda x: x.start + 1, axis=1)})
    cooler.io.append(op.join(tmp, 'test.df.2000kb.cool'), 
                     'bins', 
                     {'ones': 1})
    series = (ddf['count'] / ddf['count'].sum())
    series.name = 'normed'
    cooler.io.append(op.join(tmp, 'test.df.2000kb.cool'), 
                     'pixels', 
                     series)
    cooler.io.append(op.join(tmp, 'test.df.2000kb.cool'), 
                     'pixels', 
                     series, force=True)
    cooler.io.append(op.join(tmp, 'test.df.2000kb.cool'), 
                     'bins', 
                     {'twos': [np.ones(1000, dtype=int)*2, 
                               np.ones(561, dtype=int)*2]}, 
                     chunked=True, force=True)
    c2 = cooler.Cooler(op.join(tmp, 'test.df.2000kb.cool'))
    assert (len(c2.bins().columns) == 6)
    assert (len(c2.pixels().columns) == 4)
