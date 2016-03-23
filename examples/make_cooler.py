from __future__ import division, print_function
from multiprocessing import Pool
from collections import OrderedDict

import numpy as np
import pandas
import h5py

import Bio.Restriction as biorst
import Bio.Seq as bioseq
import pyfaidx
import cooler


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

        frags = pandas.DataFrame({
            'chrom': [chrom] * n_frags,
            'start': cuts[:-1],
            'end': cuts[1:]},
            columns=['chrom', 'start', 'end'])
        return frags

    return pandas.concat(map(_each, chroms), axis=0, ignore_index=True)


if __name__ == '__main__':

    ENZYME = 'HindIII'
    FASTA_PATH = 'fasta file'
    CHROMINFO_PATH = 'UCSC chromInfo-like file'
    HIC_PATH = ('filtered, merged and **SORTED** input HDF5 file containing'
                'datasets: chrms1 cuts1 chrms2 cuts2')
    COOLER_PATH = 'output binned sparse contact map file path'
    BINSIZE = 'integer or "frag"'
    N_CPUS = 4


    # Index and read a single FASTA file
    # If using multiple files, read them separately and put the records into one
    # ordered dictionary.
    # Pyfaidx will autogenerate fai index files.
    fasta_records = OrderedDict(pyfaidx.Fasta(FASTA_PATH))


    # Need a chromInfo.txt style tab-separated file
    # Two columns: 1) chromosome label and 2) length in bp. 
    # (the fai file usually works)
    # Chromosomes should be listed in the same order as in the fasta records.
    chromtable = pandas.read_csv(
        CHROMINFO_PATH, sep='\t', usecols=[0, 1], names=['name', 'length'])
    chromtable.index = chromtable['name']


    if BINSIZE == 'frag':
        # Make a fragment-level "bin table" 
        fragtable = digest(fasta_records, ENZYME)

        # Bin the data (non-uniform fragment-level binning), 
        # i.e. bintable == fragtable
        # Note that matrix balancing does not yet support non-uniformly binned
        # data
        h5opts = {'compression': 'gzip', 'compression_opts': 6}
        chunksize = int(100e6)
        with h5py.File(HIC_PATH, 'r') as h5read:
            with h5py.File(COOLER_PATH, 'w') as h5binned:
                cooler.io.from_readhdf5(
                    h5binned, 
                    chromtable, 
                    fragtable,
                    h5read, 
                    info={'genome-assembly': 'myAssembly'},
                    h5opts=h5opts,
                    chunksize=chunksize)

    else:
        # For uniform bins, no need to assign fragments, just use:
        bintable = cooler.make_bintable(chromtable['length'], BINSIZE)

        # Bin the data
        h5opts = {'compression': 'gzip', 'compression_opts': 6}
        chunksize = int(100e6)
        with h5py.File(HIC_PATH, 'r') as h5frag:
            with h5py.File(COOLER_PATH, 'w') as h5binned:
                cooler.io.from_readhdf5(
                    h5binned, 
                    chromtable, 
                    bintable,
                    h5read,
                    binsize=BINSIZE,
                    info={'genome-assembly': 'myAssembly'},
                    h5opts=h5opts,
                    chunksize=chunksize)

        # Compute a genome-wide balancing/bias/normalization vector
        # *** assumes uniform binning ***
        from cooler import balancing
        chunksize = int(100e6)
        try:
            pool = Pool(N_CPUS)
            with h5py.File(COOLER_PATH, 'a') as h5:
                bias = balancing.iterative_correction(
                    h5, chunksize=chunksize, tol=1e-05, min_nnz=100,
                    cis_only=False, ignore_diags=3, map=pool.map)
                # add the bias column to the file (optional)
                #h5['bins'].create_dataset('weight', data=bias, **h5opts)
        finally:
            pool.close()

        """
            # example range query + applying balancing weights
            c = cooler.Cooler(COOLER_PATH)
            
            # fetch a scipy sparse matrix
            mat = c.matrix().fetch('chr1:20,000,000-40,000,000')

            # apply the balancing weights
            i0, i1 = c.extent('chr1')
            b = bias[i0:i1]
            mat.data = b[mat.row] * b[mat.col] * mat.data
            
            # convert to dense numpy array
            A = mat.toarray()
            np.save('chr1.npy', A)

        """
