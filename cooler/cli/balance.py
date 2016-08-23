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
    "cooler_path")
@click.option(
    "--ncpu", "-p",
    type=int,
    default=8)
@click.option(
    "--chunksize", "-c",
    type=int,
    default=int(100e6))
@click.option(
    "--mad-max",
    type=int,
    default=0)
@click.option(
    "--min-nnz",
    type=int,
    default=0)
@click.option(
    "--min-count",
    type=int,
    default=0)
@click.option(
    "--ignore-diags",
    type=int,
    default=3)
@click.option(
    "--tol",
    type=float,
    default=1e-5)
@click.option(
    "--cis-only",
    is_flag=True,
    default=False)
@click.option(
    "--force", "-f",
    help="Overwrite the target dataset 'weight' if it already exists",
    is_flag=True,
    default=False)
def balance(cooler_path, ncpu, chunksize, mad_max, min_nnz, min_count, ignore_diags, tol, cis_only, force):
    """
    Compute a genome-wide balancing/bias/normalization vector. Assumes uniform binning.

    COOLER_PATH : Cooler file

    """
    with h5py.File(cooler_path, 'r') as h5:
        if 'weight' in h5['bins']:
            if not force:
                print("'weight' column already exists. Use --force option to overwrite.", file=sys.stderr)
                sys.exit(1)
            else:
                del h5['bins']['weight']     

    try:
        pool = Pool(ncpu)
        with h5py.File(cooler_path, 'a') as h5:

            bias = ice.iterative_correction(
                h5,
                chunksize=chunksize,
                cis_only=cis_only,
                tol=tol,
                min_nnz=min_nnz,
                min_count=min_count,
                mad_max=mad_max,
                ignore_diags=ignore_diags,
                map=pool.map)

            # add the bias column to the file
            h5opts = dict(compression='gzip', compression_opts=6)
            h5['bins'].create_dataset('weight', data=bias, **h5opts)

    finally:
        pool.close()
