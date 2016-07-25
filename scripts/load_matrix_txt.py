from __future__ import division, print_function
import argparse
import sys

import numpy as np
import pandas as pd
import h5py

import cooler


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Build a contact matrix by aggregating contacts.")
    parser.add_argument(
        "bins",
        help="BED-like file containing genomic bin segmentation",
        metavar="BINS_PATH")
    parser.add_argument(
        "pixels",
        help="Non-zero aggregated contact counts",
        metavar="PIXELS_PATH")
    parser.add_argument(
        "out",
        help="Output cooler file"
        metavar="COOLER_PATH")
    args = vars(parser.parse_args())

    # Bin table
    bins = pd.read_csv(args['bins'], sep='\t')
    binsizes = (bins.end - bins.start).unique()
    if len(binsizes) == 1:
        binsize = int(binsize)
    else:
        binsize = None

    # Chrom table
    chroms = bins.drop_duplicates(['chrom'], keep='last')[['chrom', 'end']]
    chroms = chroms.reset_index(drop=True)
    chroms = chroms.rename(columns={'chrom': 'name', 'end': 'length'})

    pixels = pd.read_csv(args['pixels'], sep='\t', compression='gzip',
        names=['bin1_id', 'bin2_id', 'count'])

    # Aggregate the contacts
    h5opts = {'compression': 'gzip', 'compression_opts': 6, 'chunks': True}
    chunksize = int(100e6)
    with h5py.File(args['out'], 'w') as h5:
        cooler.io.from_sparse(
            h5,
            chroms,
            bins,
            pixels,
            binsize=binsize,
            info={'genome-assembly': args.get('assembly', 'unknown')},
            h5opts=h5opts,
            chunksize=chunksize)

