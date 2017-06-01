#!/usr/bin/env python
from __future__ import division, print_function
import argparse
import gzip
import sys

import numpy as np
import cooler


def main():
    parser = argparse.ArgumentParser(
        description="Convert a cooler file to a tsv file.")
    parser.add_argument(
        "cooler_file",
        help="Cooler file",
        metavar="COOLER_PATH")
    parser.add_argument(
        "--bins", "-b",
        help="Dump the bin table.",
        action='store_true',
        default=False)
    parser.add_argument(
        "--join",
        help="Print chromosome bin coordinates instead of bin IDs",
        action='store_true',
        default=False)
    parser.add_argument(
        "--balanced",
        help="Apply balancing weights to data",
        action='store_true',
        default=False)
    parser.add_argument(
        '--chunksize',
        help="Sets the amount of data loaded from disk at one time",
        type=int,
        default=int(1e6))
    parser.add_argument(
        "--out", "-o",
        help="Output text file")
    args = vars(parser.parse_args())


    c = cooler.Cooler(args['cooler_file'])
    chunksize = args['chunksize']
    
    # output stream
    out = args['out']
    if out is None or out == '-':
        f = sys.stdout
    elif out.endswith('.gz'):
        f = gzip.open(out, 'wt')
    else:
        f = open(out, 'wt')

    # choose the source
    if args['bins']:
        selector = c.bins()
        n = c.info['nbins']
    else:
        selector = c.pixels()
        n = c.info['nnz']
        bins = c.bins()[:]  # load all the bins
        if args['balanced'] and 'weight' not in bins.columns:
            print('Balancing weights not found', file=sys.stderr)
            sys.exit(1)

    # write in chunks
    edges = np.arange(0, n+chunksize, chunksize)
    edges[-1] = n

    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = selector[lo:hi]

        if not args['bins']:

            # apply any relational joins
            cols_to_drop = []

            if args['join']:
                ncols = len(sel.columns)
                sel = (
                    sel.merge(bins[['chrom', 'start', 'end']],
                              how='left',
                              left_on='bin1_id', 
                              right_index=True)
                       .merge(bins[['chrom', 'start', 'end']],
                              how='left',
                              left_on='bin2_id', 
                              right_index=True, 
                              suffixes=('1', '2'))
                )
                sel = sel[list(sel.columns[ncols:]) + list(sel.columns[:ncols])]
                cols_to_drop.extend(['bin1_id', 'bin2_id'])

            if args['balanced']:
                sel = (
                    sel.merge(bins[['weight']],
                              how='left', 
                              left_on='bin1_id',
                              right_index=True)
                       .merge(bins[['weight']],
                              how='left',
                              left_on='bin2_id',
                              right_index=True,
                              suffixes=('1', '2'))
                    )
                sel['balanced'] = sel['weight1'] * sel['weight2'] * sel['count']
                cols_to_drop.extend(['weight1', 'weight2', 'count'])

            if cols_to_drop:
                sel = sel.drop(cols_to_drop, axis=1)

        try:
            sel.to_csv(f, sep='\t', index=False, header=True)

        except OSError as e:
            if e.errno == 32:  # broken pipe
                try:
                    f.close()
                except OSError:
                    pass
                break
            else:
                raise
    else:
        f.close()


if __name__ == '__main__':
    main()
