#! /usr/bin/env python
import sys
import h5py

infiles = sys.argv[1:]

for infile in infiles:
    with h5py.File(infile, 'a') as h5:
        print(infile)
        if 'matrix' in h5:
            print('renaming matrix')
            h5['pixels'] = h5['matrix']
            #del h5['matrix']

        if 'scaffolds' in h5:
            print('renaming scaffolds')
            h5['chroms'] = h5['scaffolds']
            #del h5['scaffolds']


