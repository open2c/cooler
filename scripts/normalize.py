from __future__ import division, print_function
from collections import OrderedDict
from multiprocessing import Pool

import numpy as np
import pandas
import h5py

import pyfaidx
import cooler


N_CPUS = 4

# TODO: bundle a clean h5dict wrapper?

if __name__ == '__main__':
    # Compute a genome-wide balancing/bias/normalization vector
    # *** assumes uniform binning ***
    from cooler import balancing
    chunksize = int(100e6)
    try:
        pool = Pool(N_CPUS)
        with h5py.File(COOLER_PATH, 'a') as h5:
            bias = balancing.iterative_correction(
                h5, chunksize=chunksize, tol=1e-05, min_nnz=100,
                cis_only=False, ignore_diags=3, map=pool.map)

            # add the bias column to the file
            # TODO: if already exists, provide dialog to overwrite
            h5['bins'].create_dataset('weight', data=bias, **h5opts)
    finally:
        pool.close()

