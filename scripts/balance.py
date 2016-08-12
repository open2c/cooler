#!/usr/bin/env python
from __future__ import division, print_function
from multiprocessing import Pool
import argparse

import numpy as np
import h5py

import cooler
import cooler.ice



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Compute a genome-wide balancing/bias/normalization vector. Assumes uniform binning.")
    parser.add_argument(
        "cooler_file",
        help="Cooler file",
        metavar="COOLER_PATH")
    parser.add_argument(
        "--nproc", "-p",
        type=int,
        default=8)
    parser.add_argument(
        "--chunksize", "-c",
        type=int,
        default=int(100e6))
    parser.add_argument(
        "--mad-max",
        type=int,
        default=0)
    parser.add_argument(
        "--min-nnz",
        type=int,
        default=0)
    parser.add_argument(
        "--min-count",
        type=int,
        default=0)
    parser.add_argument(
        "--ignore-diags",
        type=int,
        default=3)
    parser.add_argument(
        "--tol",
        type=float,
        default=1e-5)
    parser.add_argument(
        "--cis-only",
        action='store_true',
        default=False)

    args = vars(parser.parse_args())
    chunksize = args['chunksize']
    N_CPUS = args['nproc']

    try:
        pool = Pool(N_CPUS)
        with h5py.File(args['cooler_file'], 'a') as h5:
            bias = cooler.ice.iterative_correction(
                h5,
                chunksize=chunksize,
                cis_only=args['cis_only'],
                tol=args['tol'],
                min_nnz=args['min_nnz'],
                min_count=args['min_count'],
                mad_max=args['mad_max'],
                ignore_diags=args['ignore_diags'],
                map=pool.map)

            # add the bias column to the file
            if 'weight' in h5['bins']:
                del h5['bins']['weight']
            h5opts = dict(compression='gzip', compression_opts=6)
            h5['bins'].create_dataset('weight', data=bias, **h5opts)

    finally:
        pool.close()
