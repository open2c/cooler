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
        "--join",
        help="Print chromosome bin coordinates instead of bin IDs",
        action='store_true',
        default=False)
    parser.add_argument('--chunksize',
            help="Reduce the amount of data loaded by cooler at one time",
            type=int,
            default=100e6)
    parser.add_argument(
        "--out", "-o",
        help="Output text file")
    args = vars(parser.parse_args())

    c = cooler.Cooler(args['cooler_file'])
    table = c.pixels(join=args['join'])

    out = args['out']
    if out is None or out == '-':
        f = sys.stdout
    elif out.endswith('.gz'):
        f = gzip.open(out, 'wt')
    else:
        f = open(out, 'wt')

    chunksize = int(args['chunksize'])
    print("chunksize:", chunksize, file=sys.stderr)
    spans = np.arange(0, c.info['nnz']+chunksize, chunksize)

    for lo, hi in zip(spans[:-1], spans[1:]):
        print(lo, hi, file=sys.stderr)

        pix = table[lo:hi]

        try:
            pix.to_csv(f, sep='\t', index=False, header=False)
        except OSError:
            print("File probably closed:", lo, hi, file=sys.stderr)
            f.close()
            break;


if __name__ == '__main__':
    main()
