#!/usr/bin/env python
from __future__ import division, print_function
import argparse
import sys

import pandas as pd
import cooler


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Output a genome segmentation of restriction fragments as a BED file.")
    parser.add_argument(
        "chromsizes",
        help="UCSC-like chromsizes file, with chromosomes in desired order",
        metavar="CHROMSIZES_PATH")
    parser.add_argument(
        "--binsize",
        help="Resolution (bin size) in base pairs <int>")
    parser.add_argument(
        "--out", "-o",
        help="Output file (defaults to stdout)")
    args = vars(parser.parse_args())

    binsize = int(args['binsize'])

    chroms = pd.read_csv(
        args['chromsizes'], sep='\t', usecols=[0, 1], names=['name', 'length'])
    chroms.index = chroms['name']
    chromsizes = chroms['length']

    bins = cooler.make_bintable(chromsizes, binsize)

    # Write output
    out = args['out']
    try:
        if out is None:
            f = sys.stdout
        else:
            f = open(out, 'wt')
        bins.to_csv(f, sep='\t', index=False)
    finally:
        if f is not sys.stdout:
            f.close()

