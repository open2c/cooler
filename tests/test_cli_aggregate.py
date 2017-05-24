# -*- coding: utf-8 -*-
from functools import partial
import os.path as op
import subprocess
import traceback
import tempfile
import filecmp
import sys
import os
import cooler
import logging

import numpy as np
import h5py

import nose
from nose.tools import with_setup, set_trace
from click.testing import CliRunner

from cooler.cli.cload import cload, tabix as cload_tabix
from cooler.cli.aggregate import coarsen, zoomify, multires_aggregate


if sys.version_info[0] == 3 and sys.version_info[1] == 3:
    raise nose.SkipTest


testdir = op.realpath(op.dirname(__file__))
tmp = tempfile.gettempdir()

multires_path = op.join(tmp, 'test.multires.cool')


def teardown_func(*filepaths):
    for fp in filepaths:
        try:
            os.remove(fp)
        except OSError:
            pass


def test_recursive_agg():
    infile = op.join(testdir, 'data', 'GM12878-MboI-matrix.2000kb.cool')
    outfile = '/tmp/bla.cool'
    chunksize = int(10e6)
    n_zooms = 2
    n_cpus = 8
    multires_aggregate(infile, outfile, n_cpus, chunksize)
    #ccc.multires_balance(outfile, n_zooms, chunksize, n_cpus)


@with_setup(teardown=partial(teardown_func, multires_path))
def test_zoomify():
    runner = CliRunner()
    result = runner.invoke(
        zoomify, [
            op.join(testdir, 'data', 
                'dec2_20_pluslig_1pGene_grch38_UBR4_D_1nt.pairwise.sorted.cool'),
            '--out', multires_path,
            '--balance', '--balance-args', '--mad-max 3 --min-nnz 100'
        ]
    )

    #sys.stdout.write(result.output)
    assert result.exit_code == 0 , ''.join(traceback.format_exception(*result.exc_info))

    # this file should have base + 6 zoom levels
    assert(len(cooler.io.ls(multires_path)) == 7)

    # inconsistent chromosome names in chrom table (truncated) and bin table 
    # (full length) of the input file are now resolved by forcing use of the 
    # chrom table names in the bin tables of the output file
    c = cooler.Cooler(multires_path + '::' + '1')
    names = c.bins()['chrom'][:].cat.categories
    assert names[0] == 'ENSG00000127481|ENST00000375254|'

    # FIXME: with the exception of the base resolution
    c = cooler.Cooler(multires_path + '::' + '6')
    names = c.bins()['chrom'][:].cat.categories
    assert names[0] != 'ENSG00000127481|ENST00000375254|'
