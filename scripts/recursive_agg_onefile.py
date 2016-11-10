from __future__ import division, print_function
from multiprocessing import Pool
import argparse
import sys

import numpy as np
from cooler.io import CoolerAggregator
import cooler.ice
import cooler
import h5py


FACTOR = 2
TILESIZE = 256
N_CPU = 8


def set_postmortem_hook():
    import sys, traceback, ipdb
    def _excepthook(exc_type, value, tb):
        traceback.print_exception(exc_type, value, tb)
        print()
        ipdb.pm()
    sys.excepthook = _excepthook

set_postmortem_hook()


def main(infile, outfile, chunksize):
    c = cooler.Cooler(infile)
    binsize = c.info['bin-size']
    chromtable = c.chroms()[:]
    chromsizes = chromtable.set_index('name')['length']
    chroms = chromtable['name'].values
    lengths = chromtable['length'].values
    total_length = np.sum(chromsizes.values)
    n_tiles = total_length / binsize / TILESIZE
    n_zooms = int(np.ceil(np.log2(n_tiles)))

    print(
        "Copying base matrix to level {0} and producing {0} zoom levels starting from 0...".format(n_zooms),
        file=sys.stderr
    )

    # transfer base matrix
    with h5py.File(outfile, 'w') as dest, \
         h5py.File(infile, 'r') as src:

        zoomLevel = str(n_zooms)
        src.copy('/', dest, zoomLevel)

        print(zoomLevel, file=sys.stderr)


    # produce aggregations
    with h5py.File(outfile, 'r+') as f:
        grp = f[str(n_zooms)]
        c = cooler.Cooler(grp)
        binsize = cooler.info(grp)['bin-size']

        for i in range(n_zooms - 1, -1, -1):
            zoomLevel = str(i)

            # aggregate
            new_binsize = binsize * FACTOR
            new_bins = cooler.util.binnify(chromsizes, new_binsize)
     
            reader = CoolerAggregator(c, new_bins, chunksize)
            
            grp = f.create_group(zoomLevel)
            f.attrs[zoomLevel] = new_binsize
            cooler.io.create(grp, chroms, lengths, new_bins, reader)

            # balance
            #with Pool(N_CPU) as pool:
            too_close = 20000  # for HindIII
            ignore_diags = max(int(np.ceil(too_close / new_binsize)), 3)

            bias = cooler.ice.iterative_correction(
                f, zoomLevel,
                chunksize=chunksize,
                min_nnz=10,
                mad_max=3,
                ignore_diags=ignore_diags,
                map=map)
            h5opts = dict(compression='gzip', compression_opts=6)
            grp['bins'].create_dataset('weight', data=bias, **h5opts)

            print(zoomLevel, file=sys.stderr)

            c = cooler.Cooler(grp)
            binsize = new_binsize


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Recursively aggregate a single resolution cooler file into a multi-resolution file.")
    parser.add_argument(
        "cooler_file",
        help="Cooler file",
        metavar="COOLER_PATH")
    parser.add_argument(
        "--out", "-o",
        help="Output text file")
    args = vars(parser.parse_args())


    infile = args['cooler_file']
    if args['out'] is None:
        outfile = infile.replace('.cool', '.multires.cool')
    else:
        outfile = args['out']

    chunksize = int(1e6)
    main(infile, outfile, chunksize)

