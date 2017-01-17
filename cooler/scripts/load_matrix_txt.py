#!/usr/bin/env python
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
        help="Output cooler file",
        metavar="COOLER_PATH")
    args = vars(parser.parse_args())

    # Bin table
    bins = pd.read_csv(args['bins'], sep='\t')

    # Chrom table
    chromtable = (
        bins.drop_duplicates(['chrom'], keep='last')[['chrom', 'end']]
            .reset_index(drop=True)
            .rename(columns={'chrom': 'name', 'end': 'length'})
    )
    chroms, lengths = list(chromtable['name']), list(chromtable['length'])

    # Load the binned contacts
    chunksize = int(100e6)
    reader = cooler.io.SparseLoader(args['pixels'], chunksize)
    with h5py.File(args['out'], 'w') as h5:
        cooler.io.create(h5, chroms, lengths, bins, reader) # metadata, assembly)
