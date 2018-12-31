# -*- coding: utf-8 -*-
import os.path as op
import tempfile
import shutil
import os

import numpy as np
import h5py

from _common import isolated_filesystem, cooler_cmp
from cooler import fileops
import cooler
import pytest


testdir = op.realpath(op.dirname(__file__))


def test_ls():
    listing = fileops.list_coolers(op.join(testdir, 'data', 'toy.symm.upper.2.mcool'))
    paths = set(listing)
    for path in (
            '/resolutions/2',
            '/resolutions/4',
            '/resolutions/8',
            '/resolutions/16',
            '/resolutions/32'):
        assert path in paths


def test_cp():
    with isolated_filesystem() as fs:
        src_file = op.join(testdir, 'data', 'toy.symm.upper.2.mcool')

        # file-to-file
        src_uri = src_file + '::resolutions/2'
        fileops.cp(src_uri, 'test.2.cool')
        cooler_cmp(src_uri, 'test.2.cool')

        # within-file
        test_file = 'test.src.mcool'
        shutil.copyfile(src_file, test_file)
        fileops.cp(test_file + '::resolutions/2', test_file + '::abc/d')
        cooler_cmp(test_file + '::resolutions/2', test_file + '::abc/d')
        with h5py.File(test_file) as f:
            assert 'resolutions/2' in f
            assert 'abc/d' in f
            assert f['resolutions/2'].id != f['abc/d'].id


def test_mv():
    with isolated_filesystem():
        ref_file = 'test.ref.mcool'
        src_file = 'test.src.mcool'
        shutil.copyfile(op.join(testdir, 'data', 'toy.symm.upper.2.mcool'), ref_file)
        shutil.copyfile(op.join(testdir, 'data', 'toy.symm.upper.2.mcool'), src_file)
        fileops.mv(src_file + '::resolutions/2', src_file + '::abc/d')
        with h5py.File(src_file) as f:
            assert 'resolutions/2' not in f
            assert 'abc/d' in f
        cooler_cmp(ref_file + '::resolutions/2', src_file + '::abc/d')


def test_ln():
    with isolated_filesystem() as fs:
        src_file = op.join(testdir, 'data', 'toy.symm.upper.2.mcool')

        # within-file hard link
        test_file = 'test.hardlink.mcool'
        shutil.copyfile(src_file, test_file)
        fileops.ln(test_file + '::resolutions/2', test_file + '::abc/d')
        with h5py.File(test_file) as f:
            assert 'resolutions/2' in f
            assert 'abc/d' in f
            assert f['resolutions/2'].id == f['abc/d'].id
        cooler_cmp(test_file + '::resolutions/2', test_file + '::abc/d')

        # within-file soft link
        test_file = 'test.softlink.mcool'
        shutil.copyfile(src_file, test_file)
        fileops.ln(test_file + '::resolutions/2', test_file + '::abc/d', soft=True)
        with h5py.File(test_file) as f:
            assert 'resolutions/2' in f
            assert 'abc/d' in f
            assert f['resolutions/2'].id == f['abc/d'].id
        cooler_cmp(test_file + '::resolutions/2', test_file + '::abc/d')

        # between-file external link
        test_file = 'test.extlink.mcool'
        dst_file = 'test.dst.cool'
        shutil.copyfile(src_file, test_file)
        fileops.ln(test_file + '::resolutions/2', dst_file + '::abc/d', soft=True)
        cooler_cmp(test_file + '::resolutions/2', dst_file + '::abc/d')


def test_read_attr_tree():
    src_file = op.join(testdir, 'data', 'toy.symm.upper.2.mcool')
    with h5py.File(src_file, 'r') as f:
        attr_dict = fileops.read_attr_tree(f, level=None)
    assert attr_dict['@attrs']['format'] == 'HDF5::MCOOL'
