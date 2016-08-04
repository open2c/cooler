#!/usr/bin/env python
import sys
import h5py

infiles = sys.argv[1:]

for infile in infiles:
    with h5py.File(infile, 'a') as h5:
        print(infile)
        if 'format-version' in h5.attrs and h5.attrs['format-version'] < 1:

            if 'matrix' in h5 and not 'pixels' in h5:
                print('renaming matrix --> pixels')
                h5['pixels'] = h5['matrix']

            if 'scaffolds' in h5 and not 'chroms' in h5:
                print('renaming scaffolds --> chroms')
                h5['chroms'] = h5['scaffolds']

            h5.attrs['format-version'] = 1

