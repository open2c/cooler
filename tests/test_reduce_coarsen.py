# -*- coding: utf-8 -*-
import os.path as op
import tempfile
import os
import cooler

import numpy as np
from cooler.reduce import coarsen, zoomify, multires_aggregate


testdir = op.realpath(op.dirname(__file__))
tmp = tempfile.gettempdir()
multires_path = op.join(tmp, 'test.multires.cool')


def test_recursive_agg():
    infile = op.join(testdir, 'data', 'hg19.GM12878-MboI.matrix.2000kb.cool')
    outfile = multires_path
    chunksize = int(10e6)
    n_zooms = 2
    n_cpus = 8
    multires_aggregate(infile, outfile, n_cpus, chunksize)

    try:
        os.remove(multires_path)
    except OSError:
        pass


# def test_zoomify():
#     zoomify(
#         op.join(testdir, 'data',
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
