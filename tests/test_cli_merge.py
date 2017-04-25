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

import cooler
from cooler.cli.merge import merge


testdir = op.realpath(op.dirname(__file__))
tmp = tempfile.gettempdir()
merged_path = op.join(tmp, 'test_merged.cool')


def teardown_func(*filepaths):
    for fp in filepaths:
        try:
            os.remove(fp)
        except OSError:
            pass


@with_setup(teardown=partial(teardown_func, merged_path))
def test_merge():
    runner = CliRunner()
    result = runner.invoke(
        merge, [
            merged_path,
            op.join(testdir, 'data', 'GM12878-MboI-matrix.2000kb.cool'),
            op.join(testdir, 'data', 'GM12878-MboI-matrix.2000kb.cool'),
        ]
    )
    assert result.exit_code == 0

    single = cooler.Cooler(op.join(testdir, 'data', 'GM12878-MboI-matrix.2000kb.cool'))
    merged = cooler.Cooler(merged_path)
    assert merged.pixels()['count'][:].sum() == 2 * single.pixels()['count'][:].sum()
