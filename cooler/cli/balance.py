# -*- coding: utf-8 -*-
from __future__ import division, print_function
from multiprocessing import Pool
import sys

import numpy as np
import h5py

import click
from . import cli
from .. import ice


@cli.command()
@click.argument(
    "cool_path")
@click.option(
    "--nproc", "-p",
    help="Number of processes to split the work between.",
    type=int,
    default=8,
    show_default=True)
@click.option(
    "--chunksize", "-c",
    help="Control the number of pixels handled by each worker process at a time.",
    type=int,
    default=int(100e6),
    show_default=True)
@click.option(
    "--mad-max",
    help="Ignore bins from the contact matrix using the 'MAD-max' filter: "
         "bins whose log marginal sum is less than ``mad-max`` mean absolute "
         "deviations below the median log marginal sum of all the bins in the "
         "same chromosome.",
    type=int,
    default=3,
    show_default=True)
@click.option(
    "--min-nnz",
    help="Ignore bins from the contact matrix whose marginal number of "
         "nonzeros is less than this number.",
    type=int,
    default=10,
    show_default=True)
@click.option(
    "--min-count",
    help="Ignore bins from the contact matrix whose marginal count is less "
         "than this number.",
    type=int,
    default=0,
    show_default=True)
@click.option(
    "--ignore-diags",
    help="Number of diagonals of the contact matrix to ignore, including the "
         "main diagonal. Examples: 0 ignores nothing, 1 ignores the main "
         "diagonal, 2 ignores diagonals (-1, 0, 1), etc.",
    type=int,
    default=3,
    show_default=True)
@click.option(
    "--tol",
    help="Threshold value of variance of the marginals for the algorithm to "
         "converge.",
    type=float,
    default=1e-5,
    show_default=True)
@click.option(
    "--max-iters",
    help="Maximum number of iterations to perform if convergence is not achieved.",
    type=int,
    default=200,
    show_default=True)
@click.option(
    "--cis-only",
    help="Calculate weights against intra-chromosomal data only instead of "
         "genome-wide.",
    is_flag=True,
    default=False)
@click.option(
    "--force", "-f",
    help="Overwrite the target dataset, 'weight', if it already exists.",
    is_flag=True,
    default=False)
def balance(cool_path, nproc, chunksize, mad_max, min_nnz, min_count,
            ignore_diags, tol, cis_only, max_iters, force):
    """
    Out-of-core contact matrix balancing.

    Assumes uniform binning. See the help for various filtering options to
    ignore poorly mapped bins.

    COOL_PATH : Path to a COOL file.

    """
    with h5py.File(cool_path, 'r') as h5:
        if 'weight' in h5['bins']:
            if not force:
                print("'weight' column already exists. "
                      "Use --force option to overwrite.", file=sys.stderr)
                sys.exit(1)
            else:
                del h5['bins']['weight']

    try:
        pool = Pool(nproc)
        with h5py.File(cool_path, 'a') as h5:

            bias = ice.iterative_correction(
                h5,
                chunksize=chunksize,
                cis_only=cis_only,
                tol=tol,
                min_nnz=min_nnz,
                min_count=min_count,
                mad_max=mad_max,
                max_iters=max_iters,
                ignore_diags=ignore_diags,
                map=pool.map)

            # add the bias column to the file
            h5opts = dict(compression='gzip', compression_opts=6)
            h5['bins'].create_dataset('weight', data=bias, **h5opts)

    finally:
        pool.close()
