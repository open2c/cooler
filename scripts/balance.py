#!/usr/bin/env python
from __future__ import division, print_function
from multiprocessing import Pool

import numpy as np
import h5py

import cooler
import cooler.ice

N_CPUS = 5

if __name__ == '__main__':
    # Compute a genome-wide balancing/bias/normalization vector
    # *** assumes uniform binning ***
    chunksize = int(100e6)
    try:
        pool = Pool(N_CPUS)
        with h5py.File(COOLER_PATH, 'a') as h5:
            bias = cooler.ice.iterative_correction(
                h5, chunksize=chunksize, tol=1e-05, min_nnz=100,
                cis_only=False, ignore_diags=3, map=pool.map)

            # add the bias column to the file
            if 'weight' in h5['bins']:
                del h5['bins']['weight']
            h5['bins'].create_dataset('weight', data=bias, **h5opts)

    finally:
        pool.close()
