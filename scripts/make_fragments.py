#!/usr/bin/env python
from __future__ import division, print_function
from collections import OrderedDict
import argparse
import sys

import numpy as np
import pandas as pd

import Bio.Restriction as biorst
import Bio.Seq as bioseq
import pyfaidx


def load_fasta(chromosomes, *filepaths):
    if len(filepaths) == 1:
        fa = pyfaidx.Fasta(filepaths[0], as_raw=True)

    else:
        fa = {}
        for filepath in filepaths:
            fa.update(pyfaidx.Fasta(filepath, as_raw=True).records)

    records = OrderedDict((chrom, fa[chrom]) for chrom in chromosomes)
    return records


def digest(fasta_records, enzyme):
    # http://biopython.org/DIST/docs/cookbook/Restriction.html#mozTocId447698
    chroms = fasta_records.keys()
    try:
        cut_finder = getattr(biorst, enzyme).search
    except AttributeError:
        raise ValueError('Unknown enzyme name: {}'.format(enzyme))

    def _each(chrom):
        seq = bioseq.Seq(str(fasta_records[chrom]))
        cuts = np.r_[0, np.array(cut_finder(seq)) + 1, len(seq)].astype(int)
        n_frags = len(cuts) - 1

        frags = pd.DataFrame({
            'chrom': [chrom] * n_frags,
            'start': cuts[:-1],
            'end': cuts[1:]},
            columns=['chrom', 'start', 'end'])
        return frags

    return pd.concat(map(_each, chroms), axis=0, ignore_index=True)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Output a genome segmentation of restriction fragments as a BED file.")
    parser.add_argument("chromsizes",
                        help="UCSC-like chromsizes file, with chromosomes in desired order",
                        metavar="CHROMSIZES_PATH")
    parser.add_argument("--fasta", "-f",
                        help="Genome assembly FASTA path",
                        metavar="FASTA_PATH")
    parser.add_argument("--enzyme", "-e",
                        help="Name of restriction enzyme")
    parser.add_argument("--out", "-o",
                        help="Output file (defaults to stdout)")
    args = vars(parser.parse_args())


    # Need a chromInfo.txt style tab-separated file
    # Two columns: 1) chromosome label and 2) length in bp.
    chroms = pd.read_csv(
        args['chromsizes'], sep='\t', usecols=[0, 1], names=['name', 'length'])
    chroms.index = chroms['name']

    # Ordered mapping of chromosome sequences
    filepaths = args['fasta'].split(',')
    fasta_records = load_fasta(list(chroms['name']), *filepaths)

    # Restriction enzyme
    enzyme = args['enzyme']

    # Digest the genome
    frags = digest(fasta_records, enzyme)

    # Write output
    out = args['out']
    try:
        if out is None:
            out = sys.stdout
        else:
            out = open(out, 'wt')
        frags.to_csv(out, sep='\t', index=False)
    finally:
        if out is not sys.stdout:
            out.close()

