#!/usr/bin/env python
from __future__ import division, print_function
from multiprocess import Pool

import argparse
import sys

import numpy as np
from cooler.io import CoolerAggregator

import cooler
import h5py


FACTOR = 2
TILESIZE = 256


def main(infile, outfile, chunksize, n_cpus):
    c = cooler.Cooler(infile)
    binsize = c.info['bin-size']
    chromtable = c.chroms()[:]
    chromsizes = chromtable.set_index('name')['length']
    chroms = chromtable['name'].values
    lengths = chromtable['length'].values
    total_length = np.sum(chromsizes.values)
    n_tiles = total_length / binsize / TILESIZE
    n_zooms = int(np.ceil(np.log2(n_tiles)))

    print("binsize:", binsize)
    print("total_length (bp):", total_length)
    print(
        "Copying base matrix to level {0} and producing {0} new zoom levels counting down to 0...".format(n_zooms),
        file=sys.stderr
    )

    try:
        # If using HDF5 file in a process pool, fork before opening
        pool = Pool(n_cpus)

        print('n_zooms:', n_zooms, file=sys.stderr)

        # transfer base matrix
        with h5py.File(outfile, 'w') as dest, \
             h5py.File(infile, 'r') as src:

            zoomLevel = str(n_zooms)
            src.copy('/', dest, zoomLevel)

            binsize = src.attrs['bin-size']
            dest.attrs[str(n_zooms)] = binsize
            dest.attrs['max-zoom'] = n_zooms

            print("ZoomLevel:", zoomLevel, binsize, file=sys.stderr)
            new_binsize = binsize

        with h5py.File(outfile, 'r+') as f:
            # aggregate
            for i in range(n_zooms - 1, -1, -1):

                new_binsize *= FACTOR
                new_bins = cooler.util.binnify(chromsizes, new_binsize)

                prevLevel = str(i+1)
                zoomLevel = str(i)
                print("ZoomLevel:", zoomLevel, new_binsize, file=sys.stderr)

                c = cooler.Cooler(f[prevLevel])

                reader = CoolerAggregator(c, new_bins, chunksize, map=pool.imap)

                cooler.io.create(f.create_group(zoomLevel), chroms, lengths, new_bins, reader)
                f.attrs[zoomLevel] = new_binsize
                f.flush()

        with h5py.File(outfile, 'r+') as f:
            # balance
            for i in range(n_zooms - 1, -1, -1):
                zoomLevel = str(i)
                grp = f[zoomLevel]
                binsize = f.attrs[zoomLevel]

                # balance
                too_close = 10000  # for HindIII
                # too_close = 1000  # for DpnII
                ignore_diags = 1 + int(np.ceil(too_close / binsize))

                bias, stats = cooler.ice.iterative_correction(
                    f, zoomLevel,
                    chunksize=chunksize,
                    min_nnz=10,
                    mad_max=0,
                    ignore_diags=ignore_diags,
                    normalize_marginals=True,
                    map=pool.map)
                h5opts = dict(compression='gzip', compression_opts=6)
                grp['bins'].create_dataset('weight', data=bias, **h5opts)
                grp['bins']['weight'].attrs.update(stats)

    except Exception as e:
        print(e)

    finally:
            pool.close()


def check_ncpus(arg_value):
    arg_value = int(arg_value)

    if arg_value <= 0:
        raise argparse.ArgumentTypeError("n_cpus must be >= 1")
    else:
        return arg_value


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Recursively aggregate a single resolution cooler file into a multi-resolution file.")
    parser.add_argument(
        "cooler_file",
        help="Cooler file",
        metavar="COOLER_PATH")
    parser.add_argument(
        "--out", "-o",
        help="Output file")
    parser.add_argument(
        "--n_cpus", "-n",
        help="Number of cpus to use in process pool",
        default=1,
        type=check_ncpus)
    args = vars(parser.parse_args())

    infile = args['cooler_file']
    if args['out'] is None:
        outfile = infile.replace('.cool', '.multires.cool')
    else:
        outfile = args['out']

    chunksize = int(10e6)
    main(infile, outfile, chunksize, args["n_cpus"])
