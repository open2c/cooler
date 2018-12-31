# # -*- coding: utf-8 -*-
# from __future__ import division, print_function
# import multiprocess as mp
# import sys

# import numpy as np
# import h5py

# import click
# from . import cli
# from ..api import Cooler
# from ..tools import split
# from .. import util


# class maskdiags(object):
#     def __init__(self, n_diags):
#         self.n_diags = n_diags
#     def __call__(self, chunk):
#         pixels = chunk['pixels']
#         mask = np.abs(pixels['bin1_id'] - pixels['bin2_id']) < self.n_diags
#         pixels['count'][mask] = 0
#         return chunk


# def marginalize(chunk):
#     n = len(chunk['bins']['chrom'])
#     pixels = chunk['pixels']
#     marg = (
#           np.bincount(pixels['bin1_id'], weights=pixels['count'], minlength=n)
#         + np.bincount(pixels['bin2_id'], weights=pixels['count'], minlength=n)
#     )
#     return marg


# def mad(data, axis=None):
#     return np.median(np.abs(data - np.median(data, axis)), axis)


# @cli.command()
# @click.argument(
#     "cool_path",
#     type=click.Path(exists=True))
# @click.option(
#     "--bins",
#     type=int)
# def marg(cool_path, bins):
#     import matplotlib.pyplot as plt

#     c = Cooler(cool_path)

#     with mp.Pool(10) as pool:
#         marg = np.sum(
#             split(c, pool.map, chunksize=int(10e6))
#                 .pipe(maskdiags(2)) 
#                 .pipe(marginalize)
#                 .combine(),
#             axis=0)
#     logNzMarg = np.log10(marg[marg>0])

#     logMedMarg = np.median(logNzMarg)
#     madSigma = mad(logNzMarg)

#     plt.hist(logNzMarg, bins=bins)
#     plt.axvline(logMedMarg, c='k')
#     plt.axvline(logMedMarg - 5*madSigma, c='r', ls='--')
#     plt.axvline(logMedMarg + 5*madSigma, c='r', ls='--')
#     plt.show()
