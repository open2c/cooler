# -*- coding: utf-8 -*-
from functools import partial
import os.path as op
import subprocess
import traceback
import tempfile
import filecmp
import sys
import os

import numpy as np
import h5py

import logging
import nose
from nose.tools import with_setup, set_trace
from click.testing import CliRunner

from cooler.cli.cload import cload, tabix as cload_tabix
from cooler.cli.csort import csort


if sys.version_info[0] == 3 and sys.version_info[1] == 3:
    raise nose.SkipTest


testdir = op.realpath(op.dirname(__file__))
tmp = tempfile.gettempdir()
testcool_path = op.join(tmp, 'test.cool')
testcsort_path = op.join(tmp, 'test.sorted.txt.gz')
testtbi_path = op.join(tmp, 'test.sorted.txt.gz.tbi')


def teardown_func(*filepaths):
    for fp in filepaths:
        try:
            os.remove(fp)
        except OSError:
            pass


@with_setup(teardown=partial(teardown_func, testcsort_path, testtbi_path))
def test_csort():
    runner = CliRunner()
    result = runner.invoke(
        csort, [
            op.join(testdir, 'data', 'GM12878-MboI-contacts.subsample.shuffled.txt.gz'),
            op.join(testdir, 'data', 'hg19-chromsizes.select.txt'),
            '-i', 'tabix',
            '-c1', '1', '-p1', '2', '-s1', '3', '-c2', '4', '-p2', '5', '-s2', '6',
            '--out', testcsort_path,
        ]
    )
    assert result.exit_code == 0, ''.join(traceback.format_exception(*result.exc_info))

    ref_path = op.join(testdir, 'data', 'GM12878-MboI-contacts.subsample.sorted.txt.gz')
    retcode = subprocess.call(['zcmp', ref_path, testcsort_path])
    assert retcode == 0


@with_setup(teardown=partial(teardown_func, testcool_path))
def test_cload_tabix():
    runner = CliRunner()
    result = runner.invoke(
        cload_tabix, [
            op.join(testdir, 'data', 'hg19-bins.2000kb.bed.gz'),
            op.join(testdir, 'data', 'GM12878-MboI-contacts.subsample.sorted.txt.gz'),
            testcool_path
        ]
    )
    # set_trace()
    # import traceback
    # traceback.print_tb(result.exc_info[2])
    assert result.exit_code == 0, ''.join(traceback.format_exception(*result.exc_info))

    ref_path = op.join(testdir, 'data', 'GM12878-MboI-matrix.2000kb.cool')
    with h5py.File(testcool_path, 'r') as f1, \
         h5py.File(ref_path, 'r') as f2:
        assert np.all(f1['pixels/bin1_id'][:] == f2['pixels/bin1_id'][:])
        assert np.all(f1['pixels/bin2_id'][:] == f2['pixels/bin2_id'][:])
        assert np.all(f1['pixels/count'][:] == f2['pixels/count'][:])

    # test loading a file with long chromosome names which should raise
    # an error
    result = runner.invoke(
        cload_tabix, [
            op.join(testdir, 'data', 'UBR4_chromsize_bins.1nt.bed'),
            op.join(testdir, 'data', 'dec2_20_pluslig_1pGene_grch38_UBR4_D_1nt.pairwise.sorted.txt.gz'),
            testcool_path
        ]
    )

    # the cload command should raise a ValueError because the chromosome names are too
    # long (greater than 32 characters)
    # UPDATE: no longer a limit
    assert result.exit_code == 0, ''.join(traceback.format_exception(*result.exc_info))

