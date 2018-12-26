# -*- coding: utf-8 -*-
from functools import partial
import os.path as op
import tempfile
import logging
import sys
import os
import numpy as np
import h5py

import cooler
from cooler.reduce import merge
import pytest


testdir = op.realpath(op.dirname(__file__))
merged_path = op.join(tempfile.gettempdir(), 'test_merged.cool')


@pytest.mark.parametrize("path1,path2", [(
    op.join(testdir, 'data', 'hg19.GM12878-MboI.matrix.2000kb.cool'),
    op.join(testdir, 'data', 'hg19.GM12878-MboI.matrix.2000kb.cool')
)])
def test_merge(path1, path2):

    merge(merged_path, [path1, path2], mergebuf=int(15e6))
    single = cooler.Cooler(path1)
    merged = cooler.Cooler(merged_path)
    assert merged.pixels()['count'][:].sum() == 2 * single.pixels()['count'][:].sum()

    try:
        os.remove(merged_path)
    except OSError:
        pass
