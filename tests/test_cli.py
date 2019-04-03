"""

Useful things for debugging CLI failures:
pytest --pdb
print(result.output)
import traceback; traceback.print_tb(result.exception.__traceback__)

"""
from __future__ import absolute_import, division
from io import StringIO
from glob import glob
import os.path as op
import tempfile
import json
import os

from pandas.api import types
import numpy as np
import pandas as pd

from _common import cooler_cmp
from click.testing import CliRunner
import cooler
import pytest
import click
tmp = tempfile.gettempdir()
testdir = op.realpath(op.dirname(__file__))
datadir = op.join(testdir, 'data')


### INGEST AND AGGREGATION ###
from cooler.cli.cload import pairs as cload_pairs
from cooler.cli.load import load
from cooler.cli.merge import merge
from cooler.cli.coarsen import coarsen
from cooler.cli.zoomify import zoomify

def _run_cload_pairs(runner, binsize, extra_args):
    args = [
        op.join(datadir, 'toy.chrom.sizes') + ':' + str(binsize),
        op.join(datadir, 'toy.pairs'),
        'toy.{}.cool'.format(binsize),
        '-c1', '1',
        '-p1', '2',
        '-c2', '3',
        '-p2', '4',
        '--assembly', 'toy',
        '--chunksize', '10',
    ] + extra_args
    return runner.invoke(cload_pairs, args)


def _cmp_pixels_2_bg(f_out, f_ref, one_based_ref=True):
    # output, 1-based starts
    out_df = cooler.Cooler(f_out).pixels(join=True)[:]
    if one_based_ref:
        out_df['start1'] += 1
        out_df['start2'] += 1

    # reference
    ref_df = pd.read_csv(f_ref,
        sep='\t',
        names=['chrom1', 'start1', 'end1',
               'chrom2', 'start2', 'end2', 'count'])

    assert np.all(out_df == ref_df)


#'--no-symmetric-upper'
#'--input-copy-status', 'unique|duplex',
@pytest.mark.parametrize("ref,extra_args", [
    ('symm.upper', []), # reflect triu pairs
    ('symm.upper', ['--input-copy-status', 'unique']),  # reflect triu pairs
    ('asymm', ['--no-symmetric-upper']),
])
def test_cload_symm_asymm(ref, extra_args):
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = _run_cload_pairs(runner, 2, extra_args)
        assert result.exit_code == 0
        _cmp_pixels_2_bg(
            'toy.2.cool',
            op.join(datadir, 'toy.{}.2.bg2'.format(ref)))


#'--temp-dir', '',
#'--no-delete-temp',
#'--max-merge', '',
@pytest.mark.parametrize("ref,extra_args", [
    ('symm.upper', ['--temp-dir', '.', '--no-delete-temp']),
])
def test_cload_mergepass(ref, extra_args):
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = _run_cload_pairs(runner, 2, extra_args)
        assert result.exit_code == 0
        _cmp_pixels_2_bg(
            'toy.2.cool',
            op.join(datadir, 'toy.{}.2.bg2'.format(ref)))
        assert len(cooler.fileops.list_coolers(glob('*.cool')[0])) > 0


#'--field', '',
#'--no-count', '',
def test_cload_field():
    runner = CliRunner()
    with runner.isolated_filesystem():
        extra_args =  ['--field', 'score=7']
        result = _run_cload_pairs(runner, 2, extra_args)
        assert result.exit_code == 0
        pixels = cooler.Cooler('toy.2.cool').pixels()[:]
        assert 'count' in pixels.columns and types.is_integer_dtype(pixels.dtypes['count'])
        assert 'score' in pixels.columns and types.is_float_dtype(pixels.dtypes['score'])

        extra_args =  ['--field', 'count=7']
        result = _run_cload_pairs(runner, 2, extra_args)
        assert result.exit_code == 0
        pixels = cooler.Cooler('toy.2.cool').pixels()[:]
        assert 'count' in pixels.columns and types.is_integer_dtype(pixels.dtypes['count'])
        assert np.allclose(pixels['count'][:], 0)

        extra_args =  ['--field', 'count=7:dtype=float']
        result = _run_cload_pairs(runner, 2, extra_args)
        assert result.exit_code == 0
        pixels = cooler.Cooler('toy.2.cool').pixels()[:]
        assert 'count' in pixels.columns and types.is_float_dtype(pixels.dtypes['count'])
        assert np.allclose(pixels['count'][:], 0.2)

        extra_args =  ['--field', 'count=7:agg=min,dtype=float']
        result = _run_cload_pairs(runner, 2, extra_args)
        assert result.exit_code == 0
        pixels = cooler.Cooler('toy.2.cool').pixels()[:]
        assert 'count' in pixels.columns and types.is_float_dtype(pixels.dtypes['count'])
        assert np.allclose(pixels['count'][:], 0.1)

        ## don't implement the --no-count for now
        # extra_args =  ['--field', 'score=7:dtype=float', '--no-count']
        # result = _run_cload_pairs(runner, 2, extra_args)
        # assert result.exit_code == 0
        # pixels = cooler.Cooler('toy.2.cool').pixels()[:]
        # assert 'count' not in pixels.columns
        # assert 'score' in pixels.columns and types.is_float_dtype(pixels.dtypes['score'])


#'--metadata', '',
#'--zero-based',
#'--comment-char', '',
#'--storage-options', '',
def test_cload_other_options():
    runner = CliRunner()
    with runner.isolated_filesystem():
        meta = {'foo': 'bar', 'number': 42}
        with open('meta.json', 'w') as f:
            json.dump(meta, f)
        extra_args =  [
            '--metadata', 'meta.json',
            '--zero-based',
            '--storage-options', 'shuffle=True,fletcher32=True,compression=lzf']
        result = _run_cload_pairs(runner, 2, extra_args)
        assert result.exit_code == 0
        c = cooler.Cooler('toy.2.cool')
        assert c.info['metadata'] == meta
        with c.open('r') as h5:
            dset = h5['bins/start']
            assert dset.shuffle
            assert dset.fletcher32
            assert dset.compression == 'lzf'


def _run_load(runner, matrix_file, format, binsize, extra_args):
    args = [
        '-f', format,
        op.join(datadir, 'toy.chrom.sizes') + ':' + str(binsize),
        op.join(datadir, matrix_file),
        'toy.{}.cool'.format(binsize),
        '--assembly', 'toy',
        '--chunksize', '10',
    ] + extra_args
    return runner.invoke(load, args)


#'--no-symmetric-upper'
#'--input-copy-status', 'unique|duplex',
@pytest.mark.parametrize("ref,extra_args", [
    ('symm.upper', []),                                              # reflect triu pairs
    ('symm.upper', ['--one-based', '--input-copy-status', 'unique']),  # reflect triu pairs
    ('asymm', ['--one-based', '--no-symmetric-upper']),
])
def test_load_symm_asymm(ref, extra_args):
    runner = CliRunner()
    with runner.isolated_filesystem():
        ref = op.join(datadir, 'toy.{}.2.bg2'.format(ref))
        result = _run_load(runner, ref, 'bg2', 2, extra_args)
        assert result.exit_code == 0
        _cmp_pixels_2_bg('toy.2.cool', ref)

#'--field', '',
def test_load_field():
    runner = CliRunner()
    with runner.isolated_filesystem():
        extra_args =  ['--field', 'count=7:dtype=float']
        result = _run_load(runner, 'toy.symm.upper.2.bg2', 'bg2', 2, extra_args)
        assert result.exit_code == 0
        pixels1 = cooler.Cooler(op.join(datadir, 'toy.symm.upper.2.cool')).pixels()[:]
        pixels2 = cooler.Cooler('toy.2.cool').pixels()[:]
        assert 'count' in pixels2.columns and types.is_float_dtype(pixels2.dtypes['count'])
        assert np.allclose(pixels1['count'][:], pixels2['count'][:])


def test_load_field2():
    runner = CliRunner()
    with runner.isolated_filesystem():
        extra_args =  ['--count-as-float']
        result = _run_load(runner, 'toy.symm.upper.2.bg2', 'bg2', 2, extra_args)
        assert result.exit_code == 0
        pixels1 = cooler.Cooler(op.join(datadir, 'toy.symm.upper.2.cool')).pixels()[:]
        pixels2 = cooler.Cooler('toy.2.cool').pixels()[:]
        assert 'count' in pixels2.columns and types.is_float_dtype(pixels2.dtypes['count'])
        assert np.allclose(pixels1['count'][:], pixels2['count'][:])


def test_merge():
    runner = CliRunner()
    with runner.isolated_filesystem():
        f_in = op.join(datadir, 'toy.symm.upper.2.cool')
        result = runner.invoke(merge, [
            'toy.2.double.cool',
            f_in,
            f_in,
            '--field', 'count:dtype=int'
        ])
        assert result.exit_code == 0
        total1 = cooler.Cooler(f_in).pixels()['count'][:].sum()
        total2 = cooler.Cooler('toy.2.double.cool').pixels()['count'][:].sum()
        assert total2 == 2 * total1


def test_coarsen():
    runner = CliRunner()
    with runner.isolated_filesystem():
        f_in = op.join(datadir, 'toy.symm.upper.2.cool')
        f_ref = op.join(datadir, 'toy.symm.upper.4.cool')
        result = runner.invoke(coarsen, [
            f_in,
            '--factor', '2',
            '--nproc', '2',
            '-o', 'toy.2.coarsen_2.cool',
        ])
        assert result.exit_code == 0
        pix1 = cooler.Cooler(f_ref).pixels()['count'][:]
        pix2 = cooler.Cooler('toy.2.coarsen_2.cool').pixels()['count'][:]
        assert np.allclose(pix1, pix2)


def test_zoomify():
    runner = CliRunner()
    with runner.isolated_filesystem():
        f_in = op.join(datadir, 'toy.symm.upper.2.cool')
        result = runner.invoke(zoomify, [
            f_in,
            '--balance',
            '--nproc', '2',
            '-o', 'toy.2.mcool',
        ])
        assert result.exit_code == 0
        # pix1 = cooler.Cooler(f_ref).pixels()['count'][:]
        # pix2 = cooler.Cooler('toy.4.cool').pixels()['count'][:]
        # assert np.allclose(pix1, pix2)


### COMPUTE ###
from cooler.cli.balance import balance

def test_balance():
    runner = CliRunner()
    with runner.isolated_filesystem():
        f_in = op.join(datadir, 'toy.symm.upper.2.cool')
        result = runner.invoke(balance, [
            f_in,
            '--ignore-diags', '2',
            '--mad-max', '0',
            '--min-nnz', '0',
            '--tol', '0.05',
            '--nproc', '2',
            '--stdout'
        ])
        assert result.exit_code == 0
        assert len(result.output.split('\n')) == 32


### OUTPUT ###
from cooler.cli.info import info
from cooler.cli.dump import dump
from cooler.cli.show import show

def test_info():
    runner = CliRunner()
    with runner.isolated_filesystem():
        f_in = op.join(datadir, 'toy.symm.upper.2.cool')
        result = runner.invoke(info, [f_in,])
        assert result.exit_code == 0

def test_dump():
    runner = CliRunner()
    with runner.isolated_filesystem():
        f_in = op.join(datadir, 'toy.symm.upper.2.cool')
        result = runner.invoke(dump, [f_in,])
        assert result.exit_code == 0

        # roundtrip symm-upper data
        bins = pd.read_csv(StringIO(
            runner.invoke(dump, [f_in, '-H', '-t', 'bins']).output), sep='\t')
        pixels = pd.read_csv(
            StringIO(runner.invoke(dump, [f_in, '-H']).output), sep='\t')
        cooler.create_cooler('out.cool', bins, pixels, symmetric_upper=True)
        cooler_cmp(f_in, 'out.cool')

        # duplexed output
        pixels2 = pd.read_csv(StringIO(
            runner.invoke(dump, [f_in, '--matrix', '-H']).output), sep='\t')
        assert len(pixels2) > len(pixels)
        upper = pixels2[pixels2['bin1_id'] <= pixels2['bin2_id']].reset_index(drop=True)
        assert np.allclose(pixels, upper)

        # lower triangle
        trans_lower = pd.read_csv(StringIO(
            runner.invoke(dump, [f_in, '-H', '-r', 'chr2', '-r2', 'chr1']).output), sep='\t')
        assert len(trans_lower) == 0
        trans_lower = pd.read_csv(StringIO(
            runner.invoke(dump, [f_in, '-m', '-H', '-r', 'chr2', '-r2', 'chr1']).output), sep='\t')
        assert len(trans_lower) > 0

        # roundtrip square data
        f_in = op.join(datadir, 'toy.asymm.2.cool')
        bins = pd.read_csv(StringIO(
            runner.invoke(dump, [f_in, '-H', '-t', 'bins']).output), sep='\t')
        pixels = pd.read_csv(StringIO(
            runner.invoke(dump, [f_in, '-H']).output), sep='\t')
        cooler.create_cooler('out.cool', bins, pixels, symmetric_upper=False)
        cooler_cmp(f_in, 'out.cool')
        pixels2 = pd.read_csv(StringIO(
            runner.invoke(dump, [f_in, '--matrix', '-H']).output), sep='\t')
        assert np.allclose(pixels, pixels2)

        # for square data, -m is a no-op
        lower1 = pd.read_csv(StringIO(
            runner.invoke(dump, [f_in, '-H', '-r', 'chr2', '-r2', 'chr1']).output), sep='\t')
        lower2 = pd.read_csv(StringIO(
            runner.invoke(dump, [f_in, '-m', '-H', '-r', 'chr2', '-r2', 'chr1']).output), sep='\t')
        assert np.allclose(lower1, lower2)



def test_show():
    pass


### FILEOPS ###
def test_cp():
    pass

def test_mv():
    pass

def test_ln():
    pass

def test_tree():
    pass

def test_attrs():
    pass


### HELPERS ###
def test_csort():
    pass

def test_makebins():
    pass

def test_digest():
    pass
