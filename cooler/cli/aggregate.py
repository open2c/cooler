from __future__ import print_function

#!/usr/bin/env python
from __future__ import division, print_function
import multiprocess as mp
from six.moves import map
import tempfile
import argparse
import sys

import numpy as np
import h5py

from ..tools import lock
from ..io import CoolerAggregator, create
from .. import chroms, info, get_logger
from ..ice import iterative_correction
from ..util import binnify

logger = get_logger()

import click
from . import cli

def check_ncpus(arg_value):
    arg_value = int(arg_value)

    if arg_value <= 0:
        raise argparse.ArgumentTypeError("n_cpus must be >= 1")
    else:
        return min(arg_value, mp.cpu_count())

FACTOR = 2
TILESIZE = 256

def multires_aggregate(infile, outfile, n_zooms, chunksize, n_cpus):
    """
    Generate a multires cooler in 2X bin size increments from a base-level 
    resolution.

    """
    with h5py.File(infile, 'r') as f:
        binsize = info(f)['bin-size']
        chromtable = chroms(f)
    chromsizes = chromtable.set_index('name')['length']
    _, lengths = chromtable['name'].values, chromtable['length'].values

    logger.info(
        "Copying base matrix to level {0} and producing {0} new zoom levels counting down to 0...".format(n_zooms)
    )

    # copy base matrix
    with h5py.File(infile, 'r') as src, \
         h5py.File(outfile, 'w') as dest:

        zoomLevel = str(n_zooms)
        src.copy('/', dest, zoomLevel)

        binsize = src.attrs['bin-size']
        dest.attrs[str(n_zooms)] = binsize
        dest.attrs['max-zoom'] = n_zooms

        logger.info("Aggregating at zoom level: " + str(zoomLevel) + " bin size: " + str(binsize))
        new_binsize = binsize

    # Aggregate
    for i in range(n_zooms - 1, -1, -1):

        new_binsize *= FACTOR
        new_bins = binnify(chromsizes, new_binsize)

        prevLevel = str(i+1)
        zoomLevel = str(i)
        logger.info("Aggregating at zoom level: " + str(zoomLevel) + " bin size: " + str(new_binsize))

        # Note: If using HDF5 file in a process pool, fork before opening
        try:
            if n_cpus > 1:
                pool = mp.Pool(n_cpus)
            with h5py.File(outfile, 'r+') as fw:
                reader = CoolerAggregator(
                    outfile, 
                    new_bins, 
                    chunksize, 
                    cooler_root=prevLevel, 
                    map=pool.imap if n_cpus > 1 else map)
                
                create(
                    outfile, 
                    chromsizes,
                    new_bins, 
                    reader,
                    group=zoomLevel,
                    lock=lock)

                fw.attrs[zoomLevel] = new_binsize
                fw.flush()
        finally:
            if n_cpus > 1:
                pool.close()


def multires_balance(outfile, n_zooms, chunksize, n_cpus, too_close=10000, include_base=False):
    """
    Balance a multires file.

    Bin-level filters applied
    -------------------------
    min_nnz = 0
    mad_max = 3
    too_close : use ~ 10000 for 6-cutter, ~1000 for 4-cutter
        (determines number of diagonals to ignore)
    
    """
    logger.info("Performing matrix balancing...")
    if include_base:
        n = n_zooms
    else:
        n = n_zooms - 1

    # balance
    for i in range(n, -1, -1):
        zoomLevel = str(i)

        with h5py.File(outfile, 'r') as fr:
            binsize = fr.attrs[zoomLevel]
            ignore_diags = 1 + int(np.ceil(too_close / binsize))
            logger.info("balancing at zoom level: " + str(zoomLevel) + " bin size: " + str(binsize))
            try:
                if n_cpus > 1:
                    pool = mp.Pool(n_cpus)
                bias, stats = iterative_correction(
                    fr, zoomLevel,
                    chunksize=chunksize,
                    min_nnz=10,
                    mad_max=3,
                    ignore_diags=ignore_diags,
                    rescale_marginals=True,
                    map=pool.map if n_cpus > 1 else map)
            finally:
                if n_cpus > 1:
                    pool.close()       

        with h5py.File(outfile, 'r+') as fw:
            h5opts = dict(compression='gzip', compression_opts=6)
            grp = fw[zoomLevel]
            dset = grp['bins'].require_dataset(
                'weight', bias.shape, bias.dtype, **h5opts)
            dset[:] = bias
            dset.attrs.update(stats)

@cli.command()
@click.argument(
        'cooler_file',
        metavar="COOLER_PATH")
@click.option(
        '--output-file',
        '-o',
        help="Output multires file")
@click.option(
        '--n_cpus', '-n',
        help="Number of cpus to use in process pool (Default=1, i.e. no pool)",
        default=1,
        type=check_ncpus)
@click.option(
        "--chunk-size", "-c",
        help="Chunk size",
        default=int(10e6),
        type=int)
@click.option(
        '--balance/--no-balance',
        default=True,
        help="Don't balance each level while recursing")
def aggregate(cooler_file, output_file, n_cpus, chunk_size, balance):
    """
    Aggregation to multi-res cooler file.

    Converts a single resolution cooler file to a multi-resolution representation
    by recursively aggregating (summing) adjacent bins.

    COOL_PATH : Path to a COOL file
    """
    infile = cooler_file
    if output_file is None:
        outfile = infile.replace('.cool', '.multires.cool')
    else:
        outfile = output_file

    chunksize = chunk_size
    n_cpus = n_cpus

    with h5py.File(infile, 'r') as f:
        binsize = info(f)['bin-size']
        chromsizes= chroms(f).set_index('name')['length']
    total_length = np.sum(chromsizes.values)
    n_tiles = total_length / binsize / TILESIZE
    n_zooms = int(np.ceil(np.log2(n_tiles)))

    print("binsize:", binsize, file=sys.stderr)
    print("total_length (bp):", total_length, file=sys.stderr)
    print('n_tiles:', n_tiles, file=sys.stderr)
    print('n_zooms:', n_zooms, file=sys.stderr)
    print("balance:", balance)
    multires_aggregate(infile, outfile, n_zooms, chunk_size, n_cpus)

    if balance:
        multires_balance(outfile, n_zooms, chunk_size, n_cpus)

    pass
