# -*- coding: utf-8 -*-
import os.path as op
import tempfile
import shutil
import os

from _common import isolated_filesystem, cooler_cmp
from cooler.reduce import merge_coolers, coarsen_cooler, zoomify_cooler, legacy_zoomify
import numpy as np
import cooler
import h5py
import pytest


testdir = op.realpath(op.dirname(__file__))
datadir = op.join(testdir, 'data')


@pytest.mark.parametrize("path1,path2", [(
    op.join(datadir, 'hg19.GM12878-MboI.matrix.2000kb.cool'),
    op.join(datadir, 'hg19.GM12878-MboI.matrix.2000kb.cool')
)])
def test_merge(path1, path2):
    with isolated_filesystem():
        merge_coolers('test.cool', [path1, path2], mergebuf=int(15e6))
        single = cooler.Cooler(path1)
        merged = cooler.Cooler('test.cool')
        assert merged.pixels()['count'][:].sum() == 2 * single.pixels()['count'][:].sum()


def test_recursive_agg():
    infile = op.join(datadir, 'hg19.GM12878-MboI.matrix.2000kb.cool')
    chunksize = int(10e6)
    n_zooms = 2
    n_cpus = 1
    with isolated_filesystem():
        legacy_zoomify(infile, 'test.multires.cool', n_cpus, chunksize)


@pytest.mark.parametrize("input_uri,factor,ref_uri", [
    (op.join(datadir, 'toy.symm.upper.2.cool'), 2,
     op.join(datadir, 'toy.symm.upper.4.cool')),
    (op.join(datadir, 'toy.asymm.2.cool'), 2,
     op.join(datadir, 'toy.asymm.4.cool')),
    (op.join(datadir, 'toy.symm.upper.var.cool'), 2,
     op.join(datadir, 'toy.symm.upper.var2x.cool'))
    ])
def test_coarsen(input_uri, factor, ref_uri):
    kwargs = dict(
        chunksize=10,
        nproc=1,
        columns=None,
        dtypes=None,
        agg=None
    )
    with isolated_filesystem():
        coarsen_cooler(
            input_uri,
            'test.cool',
            factor,
            **kwargs
        )
        cooler_cmp('test.cool', ref_uri)


def test_coarsen_partitions_correctly():
    kwargs = dict(
        nproc=1,
        columns=None,
        dtypes=None,
        agg=None
    )
    with isolated_filesystem():
        f_ref = op.join(datadir, 'odd.4.cool')
        f_in = op.join(datadir, 'odd.1.cool')
        coarsen_cooler(
            f_in,
            'odd.1.coarsen_4.cool',
            factor=4,
            chunksize=2,
            **kwargs
        )
        pix1 = cooler.Cooler(f_ref).pixels()[:]
        pix2 = cooler.Cooler('odd.1.coarsen_4.cool').pixels()[:]
        assert len(pix1) == len(pix2)
        assert sum(pix2[['bin1_id', 'bin2_id']].duplicated()) == 0
        assert np.allclose(pix1, pix2)


def test_zoomify():
    kwargs = dict(
        chunksize=10,
        nproc=1,
        columns=None,
        dtypes=None,
        agg=None,
    )
    with isolated_filesystem():
        zoomify_cooler(
            op.join(datadir, 'toy.asymm.2.cool'),
            'test.2.mcool',
            resolutions=[4, 8, 16, 32],
            **kwargs
        )
        for res in [2, 4, 8, 16, 32]:
            cooler_cmp(
                'test.2.mcool::resolutions/{}'.format(res),
                op.join(datadir, 'toy.asymm.{}.cool'.format(res)))


# def test_zoomify():
#     zoomify(
#         op.join(datadir,
#             'dec2_20_pluslig_1pGene_grch38_UBR4_D_1nt.pairwise.sorted.cool'),
#         out=multires_path,
#         balance=True,
#         balance_args='--mad-max 3 --min-nnz 100'
#     )
#     # this file should have base + 6 zoom levels
#     assert(len(cooler.io.ls(multires_path)) == 7)

#     # inconsistent chromosome names in chrom table (truncated) and bin table
#     # (full length) of the input file are now resolved by forcing use of the
#     # chrom table names in the bin tables of the output file
#     c = cooler.Cooler(multires_path + '::' + '1')
#     names = c.bins()['chrom'][:].cat.categories
#     assert names[0] == 'ENSG00000127481|ENST00000375254|'

#     # FIXME: with the exception of the base resolution
#     c = cooler.Cooler(multires_path + '::' + '6')
#     names = c.bins()['chrom'][:].cat.categories
#     assert names[0] != 'ENSG00000127481|ENST00000375254|'

#     try:
#         os.remove(multires_path)
#     except OSError:
#         pass
