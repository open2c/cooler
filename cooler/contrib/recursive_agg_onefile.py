#!/usr/bin/env python
from __future__ import division, print_function
import multiprocess as mp
from six.moves import map
import tempfile
import argparse
import sys

import numpy as np
import h5py

from cooler.io import CoolerAggregator
import cooler


FACTOR = 2
TILESIZE = 256


def aggregate(infile, outfile, n_zooms, chunksize, n_cpus):
    """
    Generate a multires cooler in 2X bin size increments from a base-level 
    resolution.

    """
    with h5py.File(infile, 'r') as f:
        binsize = cooler.info(f)['bin-size']
        chromtable = cooler.chroms(f)
    chromsizes = chromtable.set_index('name')['length']
    chroms, lengths = chromtable['name'].values, chromtable['length'].values

    print(
        "Copying base matrix to level {0} and producing {0} new zoom levels counting down to 0...".format(n_zooms),
        file=sys.stderr
    )

    # copy base matrix
    with h5py.File(infile, 'r') as src, \
         h5py.File(outfile, 'w') as dest:

        zoomLevel = str(n_zooms)
        src.copy('/', dest, zoomLevel)

        binsize = src.attrs['bin-size']
        dest.attrs[str(n_zooms)] = binsize
        dest.attrs['max-zoom'] = n_zooms

        print("ZoomLevel:", zoomLevel, binsize, file=sys.stderr)
        new_binsize = binsize

    # aggregate
    for i in range(n_zooms - 1, -1, -1):

        new_binsize *= FACTOR
        new_bins = cooler.util.binnify(chromsizes, new_binsize)

        prevLevel = str(i+1)
        zoomLevel = str(i)
        print("ZoomLevel:", zoomLevel, new_binsize, file=sys.stderr)

        # Note: If using HDF5 file in a process pool, fork before opening
        try:
            if n_cpus > 1:
                pool = mp.Pool(n_cpus)
            with h5py.File(outfile, 'r+') as fw:
                reader = CoolerAggregator(
                    outfile, 
                    new_bins, 
                    chunksize, 
                    cooler_root=prevLevel, 
                    map=pool.imap if n_cpus > 1 else map)
                
                cooler.io.create(
                    fw.create_group(zoomLevel), 
                    chroms, lengths, new_bins, reader)

                fw.attrs[zoomLevel] = new_binsize
                fw.flush()
        finally:
            if n_cpus > 1:
                pool.close()


def balance(outfile, n_zooms, chunksize, n_cpus, too_close=10000, include_base=False):
    """
    Balance a multires file.

    Bin-level filters applied
    -------------------------
    min_nnz = 0
    mad_max = 3
    too_close : use ~ 10000 for 6-cutter, ~1000 for 4-cutter
        (determines number of diagonals to ignore)
    
    """
    print("Performing matrix balancing...", file=sys.stderr)
    if include_base:
        n = n_zooms
    else:
        n = n_zooms - 1

    # balance
    for i in range(n, -1, -1):
        zoomLevel = str(i)

        with h5py.File(outfile, 'r') as fr:
            binsize = fr.attrs[zoomLevel]
            ignore_diags = 1 + int(np.ceil(too_close / binsize))
            print("ZoomLevel:", zoomLevel, binsize, file=sys.stderr)
            try:
                if n_cpus > 1:
                    pool = mp.Pool(n_cpus)
                bias, stats = cooler.ice.iterative_correction(
                    fr, zoomLevel,
                    chunksize=chunksize,
                    min_nnz=10,
                    mad_max=3,
                    ignore_diags=ignore_diags,
                    normalize_marginals=True,
                    map=pool.map if n_cpus > 1 else map)
            finally:
                if n_cpus > 1:
                    pool.close()       

        with h5py.File(outfile, 'r+') as fw:
            h5opts = dict(compression='gzip', compression_opts=6)
            grp = fw[zoomLevel]
            dset = grp['bins'].require_dataset(
                'weight', bias.shape, bias.dtype, **h5opts)
            dset[:] = bias
            dset.attrs.update(stats)



def check_ncpus(arg_value):
    arg_value = int(arg_value)

    if arg_value <= 0:
        raise argparse.ArgumentTypeError("n_cpus must be >= 1")
    else:
        return min(arg_value, mp.cpu_count())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Recursively aggregate a single resolution cooler file into a multi-resolution file.")
    parser.add_argument(
        "cooler_file",
        help="Cooler file",
        metavar="COOLER_PATH")
    parser.add_argument(
        "--out", "-o",
        help="Output multires file")
    parser.add_argument(
        "--n_cpus", "-n",
        help="Number of cpus to use in process pool (Default=1, i.e. no pool)",
        default=1,
        type=check_ncpus)
    parser.add_argument(
        "--chunk-size", "-c",
        help="Chunk size",
        default=int(10e6),
        type=int)
    args = vars(parser.parse_args())


    infile = args['cooler_file']
    if args['out'] is None:
        outfile = infile.replace('.cool', '.multires.cool')
    else:
        outfile = args['out']
    chunksize = args['chunk_size']
    n_cpus = args['n_cpus']

    with h5py.File(infile, 'r') as f:
        binsize = cooler.info(f)['bin-size']
        chromsizes= cooler.chroms(f).set_index('name')['length']
    total_length = np.sum(chromsizes.values)
    n_tiles = total_length / binsize / TILESIZE
    n_zooms = int(np.ceil(np.log2(n_tiles)))

    print("binsize:", binsize, file=sys.stderr)
    print("total_length (bp):", total_length, file=sys.stderr)
    print('n_tiles:', n_tiles, file=sys.stderr)
    print('n_zooms:', n_zooms, file=sys.stderr)
    aggregate(infile, outfile, n_zooms, chunksize, n_cpus)
    balance(outfile, n_zooms, chunksize, n_cpus)
