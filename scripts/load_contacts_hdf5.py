from __future__ import division, print_function
from collections import OrderedDict
import argparse
import sys

import numpy as np
import pandas as pd
import h5py

import pyfaidx
import cooler


# pass metadata as json or yaml file
CHROMINFO_PATH = 'UCSC chromInfo-like file'
PAIRS_PATH = ('filtered, merged and **SORTED** input HDF5 file containing'
              'datasets: chrms1 cuts1 chrms2 cuts2')
COOLER_PATH = 'output binned sparse contact map file path'

# FIXME: Chromosomes currently need to be listed in the same order as in the pairs file :(
# TODO: accept optional metadata as JSON or YAML file


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Build a contact matrix by aggregating contacts.")
    parser.add_argument(
        "pairs",
        help="HDF5 read pairs file",
        metavar="PAIRS_PATH")
    parser.add_argument(
        "chromsizes",
        help="UCSC-like chromsizes file, with chromosomes in desired order",
        metavar="CHROMSIZES_PATH")
    parser.add_argument(
        "bins",
        help="BED-like file containing genomic bin segmentation",
        metavar="BINS_PATH")
    parser.add_argument(
        "out",
        help="Output cooler file"
        metavar="COOLER_PATH")

    args = vars(parser.parse_args())

    # Chrom table
    chroms = pd.read_csv(
        args['chromsizes'], sep='\t', usecols=[0, 1], names=['name', 'length'])
    chroms.index = chroms['name']
    chromsizes = chroms['length']

    # Bin table
    bins = pd.read_csv(args['bins'], sep='\t')
    binsizes = (bins.end - bins.start).unique()
    if len(binsizes) == 1:
        binsize = int(binsize)
    else:
        binsize = None

    assert np.all(chroms['name'] == bins['chrom'].unique())

    # Aggregate the contacts
    h5opts = {'compression': 'gzip', 'compression_opts': 9, 'chunks': True}
    chunksize = int(100e6)
    with h5py.File(args['pairs'], 'r') as h5pairs:
        with h5py.File(args['out'], 'w') as h5binned:
            cooler.io.from_readhdf5(
                h5binned,
                chroms,
                bins,
                h5pairs,
                binsize=binsize,
                info={'genome-assembly': args.get('assembly', 'unknown')},
                h5opts=h5opts,
                chunksize=chunksize)

