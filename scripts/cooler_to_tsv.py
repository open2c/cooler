#!/usr/bin/env python
from __future__ import division, print_function
import argparse
import gzip
import sys


def main():
    parser = argparse.ArgumentParser(description="""
    
    python cooler_to_tsv.py cooler_file output_file
​
    Convert a cooler file to a tsv file with chromosome positions and bin counts like this:
​
    chr1     5000    10000    chr1    5000    10000    1.0

""")

    parser.add_argument('cooler_file', nargs=1)
    parser.add_argument('output_file', nargs=1)
    #parser.add_argument('-o', '--options', default='yo',
    #                    help="Some option", type='str')
    #parser.add_argument('-u', '--useless', action='store_true', 
    #                    help='Another useless option')

    args = parser.parse_args()

    import cooler
    import numpy as np

    cooler_file = args.cooler_file[0]
    print("Loading cooler file... " + cooler_file, file=sys.stderr)

    c = cooler.Cooler(args.cooler_file[0])
    print("done", file=sys.stderr)

    chunksize = 10000000
    spans = np.arange(0, c.info['nnz']+chunksize, chunksize)
    table = c.pixeltable(join=True)

    if args.output_file[0] == '-':
        f = sys.stdout
    elif args.output_file[0].endswith('.gz'):
        f = gzip.open(args.output_file[0], 'wt')
    else:
        f = open(args.output_file[0], 'w')

    for lo, hi in zip(spans[:-1], spans[1:]):
        print('loading chunk {}..{}'.format(lo, hi), file=sys.stderr)
        pix = table[lo:hi]
        print('writing', file=sys.stderr)
        pix.to_csv(f, sep='\t', index=False, header=False)    


if __name__ == '__main__':
    main()
