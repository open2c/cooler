# -*- coding: utf-8 -*-
from functools import partial
import os.path as op
import tempfile
import logging
import sys
import os

import numpy as np
import h5py
import nose
from nose.tools import with_setup, set_trace
from click.testing import CliRunner
import traceback

from cooler.cli.load import load


testdir = op.realpath(op.dirname(__file__))
tmp = tempfile.gettempdir()
out_path1 = op.join(tmp, 'test1.cool')
out_path2 = op.join(tmp, 'test2.cool')


def teardown_func(*filepaths):
    for fp in filepaths:
        try:
            os.remove(fp)
        except OSError:
            pass


@with_setup(teardown=partial(teardown_func, out_path1, out_path2))
def test_load_bg2_vs_coo():
    runner = CliRunner()
    result = runner.invoke(
        load, [
            op.join(testdir, 'data', 'hg19-bins.2000kb.bed.gz'),
            op.join(testdir, 'data', 'GM12878-MboI-matrix.2000kb.bg2.blksrt.gz'),
            '-f', 'bg2',
            out_path1
        ]
    )
    assert result.exit_code == 0, ''.join(traceback.format_exception(*result.exc_info))

    result = runner.invoke(
        load, [
            op.join(testdir, 'data', 'hg19-bins.2000kb.bed.gz'),
            op.join(testdir, 'data', 'GM12878-MboI-matrix.2000kb.txt'),
            '-f', 'coo',
            out_path2
        ]
    )
    assert result.exit_code == 0, ''.join(traceback.format_exception(*result.exc_info))

    with h5py.File(out_path1, 'r') as f1, \
         h5py.File(out_path2, 'r') as f2:

        for col in ['bin1_id', 'bin2_id', 'count']:
            assert np.all(f1['pixels'][col][:] == f2['pixels'][col][:])
