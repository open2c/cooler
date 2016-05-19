#!/usr/bin/python

import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description="""
    
    python cooler_to_tsv.py cooler_file output_file

    Convert a cooler file to a tsv file with chromosome positions and bin counts like this:

    chr1     5000    10000    chr1    5000    10000    1.0
""")

    parser.add_argument('cooler_file', nargs=1)
    parser.add_argument('output_file', nargs=1)
    #parser.add_argument('-o', '--options', default='yo',
    #					 help="Some option", type='str')
    #parser.add_argument('-u', '--useless', action='store_true', 
    #					 help='Another useless option')

    args = parser.parse_args()


    import cooler
    import numpy as np

    cooler_file = args.cooler_file[0]
    print >>sys.stderr, "Loading cooler file... " + cooler_file
    c = cooler.Cooler(args.cooler_file[0])
    print >>sys.stderr, "done"

    chunksize = 1000000
    spans = np.arange(0, c.info['nnz']+chunksize, chunksize)
    table = c.pixeltable(join=True)

    if args.output_file == '-':
        f = sys.stdout
    else:
        f = open(args.output_file[0], 'w')

    for lo, hi in zip(spans[:-1], spans[1:]):
        print(lo, hi)
        pix = table[lo:hi]
        pix.to_csv(f, sep='\t', index=False, header=False)    

if __name__ == '__main__':
    main()


