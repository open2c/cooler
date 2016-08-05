#!/usr/bin/env python
from __future__ import division, print_function
import os.path as op
import subprocess
import argparse
import sys

from mirnylib.genome import Genome
from hiclib.fragmentHiC import HiCdataset


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Sort contacts by position and order the reads of each pair so that all "
                    "contacts are upper triangular with respect to the chromosome ordering "
                    "given by the chromsizes file.")
    parser.add_argument(
        "genome",
        help="hiclib genome path",
        metavar="GENOME_PATH")
    parser.add_argument(
        "pairs",
        help="HDF5 hiclib read pairs file",
        metavar="PAIRS_PATH")

    args = vars(parser.parse_args())

    genome_db = Genome(args['genome'])
    infile = args['pairs']
    if args['out'] is not None:
        outfile = args['out']
        ds = HiCdataset(outfile, genome_db, 'HindIII')
        ds.load(infile)
        ds._sortData()
    else:
        outfile = args['out']
        ds = HiCdataset(infile, genome_db, 'HindIII')
        ds._sortData()
