#!/usr/bin/env python
from __future__ import division, print_function
import argparse
import sys

import numpy as np
import pandas as pd
import h5py

import cooler


# TODO: accept optional metadata as JSON or YAML file
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Build a contact matrix by aggregating contacts.")
    parser.add_argument(
        "pairs",
        help="Sorted hiclib HDF5 read pairs file",
        metavar="PAIRS_PATH")
    parser.add_argument(
        "bins",
        help="BED-like file containing genomic bin segmentation",
        metavar="BINS_PATH")
    parser.add_argument(
        "out",
        help="Output cooler file",
        metavar="COOLER_PATH")
    args = vars(parser.parse_args())

    # Bin table
    bins = pd.read_csv(
        args['bins'],
        sep='\t',
        names=['chrom', 'start', 'end'],
        dtype={'chrom': str})

    # Chrom sizes from bin table
    chromtable = (
        bins.drop_duplicates(['chrom'], keep='last')[['chrom', 'end']]
            .reset_index(drop=True)
            .rename(columns={'chrom': 'name', 'end': 'length'})
    )
    chroms, lengths = list(chromtable['name']), list(chromtable['length'])
    chromsizes = pd.Series(index=chroms, data=lengths)

    # Aggregate the contacts
    chunksize = int(10e6)
    with h5py.File(args['pairs'], 'r') as h5pairs, \
         h5py.File(args['out'], 'w') as h5:
        reader = cooler.io.HDF5Aggregator(h5pairs, chromsizes, bins, chunksize)
        cooler.io.create(h5, chroms, lengths, bins, reader) # metadata, assembly)
