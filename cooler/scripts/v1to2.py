#!/usr/bin/env python
import sys
import numpy as np
import h5py

CHROMID_DTYPE = np.int32

infiles = sys.argv[1:]

for infile in infiles:
    with h5py.File(infile, 'a') as h5:
        print(infile)

        if 'format-version' in h5.attrs and h5.attrs['format-version'] == 1:
            chroms = h5['chroms']['name'][:].astype('U')
            idmap = dict(zip(chroms, range(len(chroms))))
            enum_dtype = h5py.special_dtype(enum=(CHROMID_DTYPE, idmap))

            if 'bins' in h5 and 'chrom_id' in h5['bins']:
                h5opts = dict(compression='gzip', compression_opts=6)
                bin_chrom_ids = h5['bins']['chrom_id'][:]
                n_bins = len(bin_chrom_ids)

                print('renaming chrom_id [int] --> chrom [enum]')

                h5['bins'].create_dataset('chrom',
                                   shape=(n_bins,),
                                   dtype=enum_dtype,
                                   data=bin_chrom_ids,
                                   **h5opts)
                del h5['bins']['chrom_id']

            h5.attrs['format-version'] = 2

