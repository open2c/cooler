# -*- coding: utf-8 -*-
from functools import partial
import os.path as op
import tempfile
import filecmp
import os
import numpy as np
import pandas as pd
import h5py

import cooler
from cooler.cli.cload import (
    tabix as cload_tabix,
    pairix as cload_pairix,
    pairs as cload_pairs
)
from cooler.cli.load import load
import pytest

tmp = tempfile.gettempdir()
testdir = op.realpath(op.dirname(__file__))
testcool_path = op.join(tmp, 'test.cool')


@pytest.mark.filterwarnings("ignore")
@pytest.mark.parametrize("bins_path,pairs_path,ref_path", [(
    op.join(testdir, 'data', 'hg19.bins.2000kb.bed.gz'),
    op.join(testdir, 'data', 'hg19.GM12878-MboI.pairs.subsample.sorted.txt.gz'),
    op.join(testdir, 'data', 'hg19.GM12878-MboI.matrix.2000kb.cool')
)])
def test_cload_tabix(bins_path, pairs_path, ref_path):

    cload_tabix.callback(
        bins_path,
        pairs_path,
        testcool_path,
        metadata=None,
        assembly='hg19',
        nproc=8,
        zero_based=False,
        max_split=2,
    )
    with h5py.File(testcool_path, 'r') as f1, \
         h5py.File(ref_path, 'r') as f2:
        assert np.all(f1['pixels/bin1_id'][:] == f2['pixels/bin1_id'][:])
        assert np.all(f1['pixels/bin2_id'][:] == f2['pixels/bin2_id'][:])
        assert np.all(f1['pixels/count'][:] == f2['pixels/count'][:])
    try:
        os.remove(testcool_path)
    except OSError:
        pass


@pytest.mark.parametrize("bins_path,pairs_path,ref_path", [(
    op.join(testdir, 'data', 'hg19.bins.2000kb.bed.gz'),
    op.join(testdir, 'data', 'hg19.GM12878-MboI.pairs.subsample.blksrt.txt.gz'),
    op.join(testdir, 'data', 'hg19.GM12878-MboI.matrix.2000kb.cool')
)])
def test_cload_pairix(bins_path, pairs_path, ref_path):

    cload_pairix.callback(
        bins_path,
        pairs_path,
        testcool_path,
        metadata=None,
        assembly='hg19',
        nproc=8,
        zero_based=False,
        max_split=2,
    )
    with h5py.File(testcool_path, 'r') as f1, \
         h5py.File(ref_path, 'r') as f2:
        assert np.all(f1['pixels/bin1_id'][:] == f2['pixels/bin1_id'][:])
        assert np.all(f1['pixels/bin2_id'][:] == f2['pixels/bin2_id'][:])
        assert np.all(f1['pixels/count'][:] == f2['pixels/count'][:])
    try:
        os.remove(testcool_path)
    except OSError:
        pass


@pytest.mark.parametrize("bins_path,pairs_path,ref_path", [(
    op.join(testdir, 'data', 'hg19.bins.2000kb.bed.gz'),
    op.join(testdir, 'data', 'hg19.GM12878-MboI.pairs.subsample.blksrt.txt.gz'),
    op.join(testdir, 'data', 'hg19.GM12878-MboI.matrix.2000kb.cool')
)])
def test_cload_pairs(bins_path, pairs_path, ref_path):

    cload_pairs.callback(
        bins_path,
        pairs_path,
        testcool_path,
        metadata=None,
        assembly='hg19',
        chunksize=int(15e6),
        zero_based=False,
        comment_char='#',
        symmetric_input='unique',
        no_symmetric_storage=False,
        field=(),
        temp_dir=None,
        no_delete_temp=False,
        storage_options=None,
        chrom1=1, pos1=2,
        chrom2=4, pos2=5,
    )
    with h5py.File(testcool_path, 'r') as f1, \
         h5py.File(ref_path, 'r') as f2:
        assert np.all(f1['pixels/bin1_id'][:] == f2['pixels/bin1_id'][:])
        assert np.all(f1['pixels/bin2_id'][:] == f2['pixels/bin2_id'][:])
        assert np.all(f1['pixels/count'][:] == f2['pixels/count'][:])
    try:
        os.remove(testcool_path)
    except OSError:
        pass


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
        return pd.DataFrame(
            {'chrom': chrom, 'start': edges[:-1], 'end': edges[1:]},
            columns=['chrom', 'start', 'end'])
    return pd.concat(map(_each, chromsizes.keys()), axis=0, ignore_index=True)


def test_from_hdf5_pairs():

    def should_not_depend_on_chunksize(chromsizes, bintable, mock_pairs):
        # try different chunk sizes
        binner = cooler.io.HDF5Aggregator(
            mock_pairs, chromsizes, bintable, chunksize=66)
        cooler.io.create(testcool_path, bintable, binner)
        with h5py.File(testcool_path, 'r') as h5:
            oc1 = h5['indexes']['chrom_offset'][:]
            ob1 = h5['indexes']['bin1_offset'][:]
            p1 = cooler.pixels(h5, join=False)

        binner = cooler.io.HDF5Aggregator(
            mock_pairs, chromsizes, bintable, chunksize=666)
        cooler.io.create(testcool_path, bintable, binner)
        with h5py.File(testcool_path, 'r') as h5:
            oc2 = h5['indexes']['chrom_offset'][:]
            ob2 = h5['indexes']['bin1_offset'][:]
            p2 = cooler.pixels(h5, join=False)

        assert np.all(oc1 == oc2)
        assert np.all(ob1 == ob2)
        assert np.all(p1.values == p2.values)


    def should_raise_if_input_not_sorted(chromsizes, bintable, mock_pairs):
        # not sorted by chrm1
        #with h5py.File(testcool_path, 'w') as h5:
        bad_reads = {
            'chrms1': mock_pairs['chrms2'],
            'cuts1':  mock_pairs['cuts2'],
            'chrms2': mock_pairs['chrms1'],
            'cuts2':  mock_pairs['cuts1'],
        }
        with pytest.raises(ValueError):
            cooler.io.HDF5Aggregator(bad_reads, chromsizes, bintable, chunksize=66)

        # not triu
        bad_reads = {
            'chrms1': mock_pairs['chrms1'].copy(),
            'cuts1':  mock_pairs['cuts1'].copy(),
            'chrms2': mock_pairs['chrms2'].copy(),
            'cuts2':  mock_pairs['cuts2'].copy(),
        }
        bad_reads['chrms1'][0] = 0
        bad_reads['chrms2'][0] = 0
        bad_reads['cuts1'][0] = 10
        bad_reads['cuts2'][0] = 9
        binner = cooler.io.HDF5Aggregator(
            bad_reads, chromsizes, bintable, chunksize=66)
        with pytest.raises(ValueError):
            cooler.io.create(testcool_path, bintable, binner)


    def should_work_with_int32_cols(chromsizes, bintable, mock_pairs):
        # int64
        binner = cooler.io.HDF5Aggregator(
            mock_pairs, chromsizes, bintable, chunksize=66)
        cooler.io.create(testcool_path, bintable, binner)
        with h5py.File(testcool_path, 'r') as h5:
            oc1 = h5['indexes']['chrom_offset'][:]
            ob1 = h5['indexes']['bin1_offset'][:]
            p1 = cooler.pixels(h5, join=False)

        # int32
        mock_pairs32 = {
            'chrms1': mock_pairs['chrms1'].astype(np.int32),
            'cuts1':  mock_pairs['cuts1'].astype(np.int32),
            'chrms2': mock_pairs['chrms2'].astype(np.int32),
            'cuts2':  mock_pairs['cuts2'].astype(np.int32),
        }
        binner = cooler.io.HDF5Aggregator(
            mock_pairs32, chromsizes, bintable, chunksize=66)
        cooler.io.create(testcool_path, bintable, binner)
        with h5py.File(testcool_path, 'r') as h5:
            oc2 = h5['indexes']['chrom_offset'][:]
            ob2 = h5['indexes']['bin1_offset'][:]
            p2 = cooler.pixels(h5, join=False)

        assert np.all(oc1 == oc2)
        assert np.all(ob1 == ob2)
        assert np.all(p1.values == p2.values)


    def _mock_hdf5_pairs():
        np.random.seed(1)
        chrms = np.random.randint(0, n_chroms, n_records * 2)
        cuts = np.random.randint(0, clen, n_records * 2)
        abs_cuts = np.array([clen * chrm + cut for chrm, cut in zip(chrms, cuts)])
        abs_cuts1, abs_cuts2 = abs_cuts[:n_records], abs_cuts[n_records:]
        mock_pairs = {
            'chrms1': chrms[:n_records],
            'cuts1':  cuts[:n_records],
            'chrms2': chrms[n_records:],
            'cuts2':  cuts[n_records:],
        }
        # Triu-sort
        mask = abs_cuts1 > abs_cuts2
        mock_pairs['chrms1'][mask], mock_pairs['chrms2'][mask] = mock_pairs['chrms2'][mask], mock_pairs['chrms1'][mask]
        mock_pairs['cuts1'][mask], mock_pairs['cuts2'][mask] = mock_pairs['cuts2'][mask], mock_pairs['cuts1'][mask]
        abs_cuts1[mask], abs_cuts2[mask] = abs_cuts2[mask], abs_cuts1[mask]
        idx = np.lexsort([abs_cuts2, abs_cuts1])
        for key in mock_pairs:
            mock_pairs[key] = mock_pairs[key][idx]
        return mock_pairs

    n_chroms = 2
    clen = 2000
    n_records = 3000
    chromsizes = pd.Series(index=['chr1', 'chr2'], data=[clen, clen])
    mock_pairs = _mock_hdf5_pairs()

    # uniform bins
    bintable = cooler.binnify(chromsizes, 100)
    should_not_depend_on_chunksize(chromsizes, bintable, mock_pairs)
    should_raise_if_input_not_sorted(chromsizes, bintable, mock_pairs)
    should_work_with_int32_cols(chromsizes, bintable, mock_pairs)

    # non-uniform bins
    bintable = _alternating_bins(chromsizes, [10, 100])
    should_not_depend_on_chunksize(chromsizes, bintable, mock_pairs)
    should_raise_if_input_not_sorted(chromsizes, bintable, mock_pairs)
    should_work_with_int32_cols(chromsizes, bintable, mock_pairs)


def test_load_bg2_vs_coo():

    out_path1 = op.join(tmp, 'test1.cool')
    out_path2 = op.join(tmp, 'test2.cool')

    load.callback(
        op.join(testdir, 'data', 'hg19.bins.2000kb.bed.gz'),
        op.join(testdir, 'data', 'hg19.GM12878-MboI.matrix.2000kb.bg2.gz'),
        out_path1,
        format='bg2',
        metadata=None,
        assembly='hg19',
        chunksize=int(20e6),
        field=(),
        count_as_float=False,
        one_based=False,
        comment_char='#',
        symmetric_input='unique',
        no_symmetric_storage=False,
        storage_options=None,
    )
    load.callback(
        op.join(testdir, 'data', 'hg19.bins.2000kb.bed.gz'),
        op.join(testdir, 'data', 'hg19.GM12878-MboI.matrix.2000kb.coo.txt'),
        out_path2,
        format='coo',
        metadata=None,
        assembly='hg19',
        chunksize=int(20e6),
        field=(),
        count_as_float=False,
        one_based=False,
        comment_char='#',
        symmetric_input='unique',
        no_symmetric_storage=False,
        storage_options=None,
    )

    with h5py.File(out_path1, 'r') as f1, \
         h5py.File(out_path2, 'r') as f2:

        for col in ['bin1_id', 'bin2_id', 'count']:
            assert np.all(f1['pixels'][col][:] == f2['pixels'][col][:])

    for fp in [out_path1, out_path2]:
        try:
            os.remove(fp)
        except OSError:
            pass
