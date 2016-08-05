#!/usr/bin/env python
from __future__ import division, print_function
import os.path as op
import argparse
import glob
import sys

import cooler


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Output a genome segmentation of restriction fragments as a BED file.")
    parser.add_argument("chromsizes",
        help="UCSC-like chromsizes file, with chromosomes in desired order",
        metavar="CHROMSIZES_PATH")
    parser.add_argument("fasta",
        help="Genome assembly FASTA file or folder containing FASTA files (uncompressed)",
        metavar="FASTA_PATH")
    parser.add_argument("enzyme",
        help="Name of restriction enzyme",
        metavar="ENZYME")
    parser.add_argument("--out", "-o",
        help="Output file (defaults to stdout)")
    args = vars(parser.parse_args())


    chromsizes = cooler.read_chromsizes(args['chromsizes'])
    chroms = list(chromsizes.keys())

    # Load sequences
    fasta = args['fasta']
    if op.isdir(fasta):
        filepaths = glob.glob(op.join(fasta, '*.fa'))
        filepaths.extend(glob.glob(op.join(fasta, '*.fasta')))
    else:
        filepaths = [fasta]
    fasta_records = cooler.util.load_fasta(chroms, *filepaths)

    # Digest sequences
    enzyme = args['enzyme']
    frags = cooler.util.digest(fasta_records, enzyme)

    # Write output
    out = args['out']
    try:
        if out is None:
            f = sys.stdout
        else:
            f = open(f, 'wt')
        frags.to_csv(f, sep='\t', index=False, header=False)
    except OSError:
        pass
    finally:
        f.close()

