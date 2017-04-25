# -*- coding: utf-8 -*-
from __future__ import division, print_function
from collections import OrderedDict
from six.moves import map
import multiprocess as mp
import os.path as op
import math
import sys

import numpy as np
import h5py

from ..io import CoolerAggregator, create, parse_cooler_uri
from ..ice import iterative_correction
from ..util import binnify
from ..tools import lock
from .. import api

import click
from click.testing import CliRunner
from . import cli, logger
from .balance import balance as balance_cmd


QUAD_TILE_SIZE_PIXELS = 256


def aggregate(input_uri, output_uri, factor, nproc, chunksize, lock):
    c = api.Cooler(input_uri)
    chromsizes = c.chromsizes
    new_binsize = c.binsize * factor
    new_bins = binnify(chromsizes, new_binsize)

    try:
        # Note: fork before opening to prevent inconsistent global HDF5 state
        if nproc > 1:
            pool = mp.Pool(nproc)

        iterator = CoolerAggregator(
            input_uri,
            new_bins,
            chunksize,
            batchsize=nproc,
            map=pool.map if nproc > 1 else map)

        create(
            output_uri,
            new_bins,
            iterator,
            lock=lock,
            append=True)

    finally:
        if nproc > 1:
            pool.close()


def get_quadtree_depth(chromsizes, binsize):
    """
    Depth of quad tree necessary to tesselate the concatenated genome with quad
    tiles such that linear dimension of the tiles is a preset multiple of the 
    genomic resolution.

    """
    tile_size_bp = QUAD_TILE_SIZE_PIXELS * binsize
    min_tile_cover = math.ceil(sum(chromsizes) / tile_size_bp)
    return int(math.ceil(np.log2(min_tile_cover)))


def multires_aggregate(input_uri, outfile, nproc, chunksize, lock=None):
    """
    Quad-tree tiling for HiGlass

    """
    infile, ingroup = parse_cooler_uri(input_uri)

    clr = api.Cooler(infile, ingroup)
    n_zooms = get_quadtree_depth(clr.chromsizes, clr.binsize)
    factor = 2

    logger.info("total_length (bp): {}".format(np.sum(clr.chromsizes)))
    logger.info("binsize: {}".format(clr.binsize))
    logger.info("n_zooms: {}".format(n_zooms))
    logger.info("quad tile cover: {}".format(2**n_zooms))
    logger.info(
        "Copying base matrix to level " +
        "{0} and producing {0} new zoom levels ".format(n_zooms) +
        "counting down to 0..."
    )

    zoom_levels = OrderedDict()
    zoomLevel = str(n_zooms)
    binsize = clr.binsize
    logger.info(
        "Zoom level: "
        + str(zoomLevel)
        + " bin size: "
        + str(binsize))
    
    # Copy base matrix
    with h5py.File(infile, 'r') as src, \
         h5py.File(outfile, 'w') as dest:

        src.copy(ingroup, dest, str(zoomLevel))
        zoom_levels[zoomLevel] = binsize
    
    # Aggregate
    # Use lock to sync read/write ops on same file
    for i in range(n_zooms - 1, -1, -1):
        prev_binsize = binsize
        binsize *= factor
        prevLevel = str(i+1)
        zoomLevel = str(i)
        logger.info(
            "Aggregating at zoom level: "
            + str(zoomLevel)
            + " bin size: "
            + str(binsize))

        aggregate(
            outfile + '::' + str(prevLevel), 
            outfile + '::' + str(zoomLevel),
            factor, 
            nproc, 
            chunksize,
            lock
        )
        zoom_levels[zoomLevel] = binsize

    with h5py.File(outfile, 'r+') as fw:
        fw.attrs.update({'max-zoom': n_zooms})
        # grp = fw.require_group('.zooms')
        # grp.attrs.update({'max-zoom': n_zooms})
        # # hard links
        # for level, res in zoom_levels.items():
        #     fw['.zooms'][level] = fw[str(res)]

    return n_zooms, zoom_levels


@cli.command()
@click.argument(
    'cool_uri',
    metavar="COOL_PATH")
@click.option(
    '--factor', '-k',
    help="Gridding factor. The contact matrix is coarsegrained by grouping "
         "each chromosomal contact block into k-by-k element tiles",
    type=int,
    default=2,
    show_default=True)
@click.option(
    '--nproc', '-n', '-p',
    help="Number of processes to use for batch processing chunks of pixels "
         "[default: 1, i.e. no process pool]",
    default=1,
    type=int)
@click.option(
    '--chunksize', '-c',
    help="Number of pixels allocated to each process",
    type=int,
    default=int(10e6),
    show_default=True)
@click.option(
    '--out', '-o',
    required=True,
    help="Output file or URI")
def coarsen(cool_uri, factor, nproc, chunksize, out):
    """
    Coarsen a contact matrix by uniformly gridding the elements of each 
    chromosomal block and summing the elements inside the grid tiles, i.e. a
    2-D histogram.

    \b\bArguments:

    COOL_PATH : Path to a COOL file or Cooler URI.

    """
    infile, _ = parse_cooler_uri(cool_uri)
    outfile, _ = parse_cooler_uri(out)
    same_file = op.realpath(infile) == op.realpath(outfile)
    aggregate(cool_uri, out, factor, nproc, chunksize, 
              lock=lock if same_file else None)


@cli.command()
@click.argument(
    'cool_uri',
    metavar="COOL_PATH")
@click.option(
    '--nproc', '-n', '-p',
    help="Number of processes to use for batch processing chunks of pixels "
         "[default: 1, i.e. no process pool]",
    default=1,
    type=int)
@click.option(
    '--chunksize', '-c',
    help="Number of pixels allocated to each process",
    type=int,
    default=int(10e6),
    show_default=True)
@click.option(
    '--balance/--no-balance',
    help="Apply balancing to each zoom level",
    default=False,
    show_default=True)
@click.option(
    '--balance-args',
    help="Additional arguments to pass to cooler balance",
    type=str)
@click.option(
    '--out', '-o',
    help="Output file or URI")
@click.pass_context
def zoomify(ctx, cool_uri, nproc, chunksize, balance, balance_args, out):
    """
    Generate zoom levels for HiGlass by recursively generating 2-by-2 element 
    tiled aggregations of the contact matrix until reaching a minimum
    dimension. The aggregations are stored in a multi-resolution file.

    \b\bArguments:

    COOL_PATH : Path to a COOL file or Cooler URI.

    """
    infile, _ = parse_cooler_uri(cool_uri)

    if out is None:
        outfile = infile.replace('.cool', '.multi.cool')
    else:
        outfile, _ = parse_cooler_uri(out)

    logger.info('Recursively aggregating "{}"'.format(cool_uri))
    logger.info('Writing to "{}"'.format(outfile))
    n_zooms, zoom_levels = multires_aggregate(cool_uri, outfile, nproc, chunksize, lock=lock)

    if balance:
        runner = CliRunner()

        if balance_args is None: 
            balance_args = []

        for level, res in reversed(zoom_levels.items()):
            uri = outfile + '::' + str(res)
            if level == str(n_zooms):
                if 'weight' in api.Cooler(uri).bins():
                    continue
            logger.info('Balancing zoom level {}, bin size {}'.format(level, res))
            result = runner.invoke(balance_cmd, args=[uri] + balance_args)
            if result.exit_code != 0:
                raise result.exception


@cli.command(context_settings={
    'ignore_unknown_options': True,
    'allow_extra_args': True})
@click.pass_context
def coarsegrain(ctx):
    """
    Deprecated in favor of separate "coarsen" and "zoomify" commands. Do not use.

    """
    click.echo(
        '"cooler coarsegrain" is deprecated.\nUse "cooler coarsen" for '
        'single aggregations.\nUse "cooler zoomify" for multiresolution '
        'aggregation.')
