# -*- coding: utf-8 -*-
from __future__ import division, print_function
from six import StringIO
import os.path as op
import tempfile
import os

import numpy as np
import pandas as pd
import h5py
from cooler.create import sanitize_pixels, sanitize_records, BadInputError
import cooler
import pytest


testdir = op.dirname(op.realpath(__file__))
tmp = tempfile.gettempdir()
testfile_path = op.join(tmp, 'test.cool')


columns = ['chrom1', 'pos1', 'strand1',
           'chrom2', 'pos2', 'strand2',
           'name', 'pair_type', 'triu']

valid_data = """chr1\t1\t+\tchr2\t100\t-\t.\tLL\t1
chr2\t99\t+\tchr1\t13\t-\t.\tLL\t0
chr2\t13\t+\tchr2\t60\t-\t.\tLL\t1
chr1\t200\t+\tchr2\t50\t-\t.\tLL\t1
chr3\t11\t+\tchr3\t40\t-\t.\tLL\t1
chr1\t234\t+\tchr3\t30\t-\t.\tLL\t1
chr3\t3\t+\tchr2\t20\t-\t.\tLL\t0
chr2\t23\t+\tchr3\t11\t-\t.\tLL\t1
chr1\t123\t+\tchr1\t200\t-\t.\tLL\t1
"""

nuisance_chroms = """chr1\t222\t+\tchr9\t200\t-\t.\tLL\t1
chr9\t222\t+\tchr9\t200\t-\t.\tLL\t1"""
oob_lower = """chr1\t-1\t+\tchr1\t10\t+\t.\tLL\t1"""
oob_upper = """chr1\t123\t+\tchr1\t301\t+\t.\tLL\t1"""

binsize = 10
chromsizes = pd.Series(
    index=['chr1', 'chr2', 'chr3'],
    data=[300, 300, 300])
bins = cooler.util.binnify(chromsizes, binsize)


def _insert_lines(d, new):
    lines = d.split('\n')
    for line in new.split('\n'):
        lines.insert(np.random.randint(len(lines)), line)
    return '\n'.join(lines)


def test_sanitize_triu_action():
    text = valid_data
    chunk = pd.read_csv(StringIO(text), sep='\t', names=columns)
    out = sanitize_records(
        bins,
        schema='pairs',
        validate=True,
        tril_action='reflect',
    )(chunk.copy())
    is_tril = ~np.array(out['triu'], dtype=bool)
    assert np.all(out.loc[is_tril, 'chrom1'] == chunk.loc[is_tril, 'chrom2'])
    assert np.all(out.loc[is_tril, 'chrom2'] == chunk.loc[is_tril, 'chrom1'])
    assert np.all(out.loc[is_tril, 'strand1'] == '+')

    text = valid_data
    chunk = pd.read_csv(StringIO(text), sep='\t', names=columns)
    out = sanitize_records(
        bins,
        schema='pairs',
        validate=True,
        tril_action='drop',
    )(chunk.copy())
    is_tril = ~np.array(out['triu'], dtype=bool)
    assert np.all(out.loc[is_tril, 'chrom1'] == chunk.loc[is_tril, 'chrom2'])
    assert np.all(out.loc[is_tril, 'chrom2'] == chunk.loc[is_tril, 'chrom1'])
    assert np.all(out.loc[is_tril, 'strand1'] == '+')
    assert len(out) == chunk['triu'].sum()

    func = sanitize_records(
        bins,
        schema='pairs',
        validate=True,
        tril_action='raise',
    )
    text = valid_data
    chunk = pd.read_csv(StringIO(text), sep='\t', names=columns)
    with pytest.raises(BadInputError):
        func(chunk)


def test_sanitize_with_strand_column():
    text = valid_data
    chunk = pd.read_csv(StringIO(text), sep='\t', names=columns)
    out = sanitize_records(
        bins,
        schema='pairs',
        validate=True,
        tril_action='reflect',
        sided_fields=('chrom', 'pos', 'strand'),
    )(chunk.copy())
    is_tril = ~np.array(out['triu'], dtype=bool)
    assert np.all(out.loc[is_tril, 'chrom1'] == chunk.loc[is_tril, 'chrom2'])
    assert np.all(out.loc[is_tril, 'chrom2'] == chunk.loc[is_tril, 'chrom1'])
    assert np.all(out.loc[is_tril, 'strand1'] == '-')


def test_sanitize_with_nuisance_records():
    text = _insert_lines(valid_data, nuisance_chroms)
    chunk = pd.read_csv(StringIO(text), sep='\t', names=columns)
    out = sanitize_records(
        bins,
        schema='pairs',
        validate=True,
        tril_action='reflect',
    )(chunk.copy())
    assert ('chr9' not in out['chrom1']) and ('chr9' not in out['chrom2'])


def test_sanitize_with_bad_records():
    func = sanitize_records(
        bins,
        schema='pairs',
        validate=True,
        tril_action='reflect',
    )

    text = _insert_lines(valid_data, oob_lower)
    chunk = pd.read_csv(StringIO(text), sep='\t', names=columns)
    with pytest.raises(BadInputError):
        func(chunk)

    text = _insert_lines(valid_data, oob_upper)
    chunk = pd.read_csv(StringIO(text), sep='\t', names=columns)
    with pytest.raises(BadInputError):
        func(chunk)
