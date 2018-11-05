# -*- coding: utf-8 -*-
from __future__ import division, print_function
from collections import OrderedDict
from six.moves import map
import multiprocess as mp
import os.path as op
import shlex
import math
import sys

import numpy as np
import h5py

from ..reduce import multires_aggregate, zoomify as _zoomify
from ..io import parse_cooler_uri, create
from ..tools import lock
from .. import api

import click
from click.testing import CliRunner
from . import cli, get_logger


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
    '--balance',
    help="Apply balancing to each zoom level. Off by default.",
    is_flag=True,
    default=False)
@click.option(
    '--balance-args',
    help="Additional arguments to pass to cooler balance",
    type=str)
@click.option(
    '--out', '-o',
    help="Output file or URI")
@click.option(
    '--resolutions', '-r',
    help="Comma-separated list of target resolutions")
def zoomify(cool_uri, nproc, chunksize, balance, balance_args, out,
            resolutions):
    """
    Generate zoom levels for HiGlass by recursively generating 2-by-2 element
    tiled aggregations of the contact matrix until reaching a minimum
    dimension. The aggregations are stored in a multi-resolution file.

    \b\bArguments:

    COOL_PATH : Path to a COOL file or Cooler URI.

    """
    from .balance import balance as balance_cmd
    logger = get_logger(__name__)
    infile, _ = parse_cooler_uri(cool_uri)

    if out is None:
        outfile = infile.replace('.cool', '.mcool')
    else:
        outfile, _ = parse_cooler_uri(out)

    logger.info('Recursively aggregating "{}"'.format(cool_uri))
    logger.info('Writing to "{}"'.format(outfile))

    if resolutions is not None:
        resolutions = [int(s.strip()) for s in resolutions.split(',')]
        _zoomify([cool_uri], outfile, resolutions, nproc,
            chunksize, lock=lock)

        if balance:
            runner = CliRunner()

            if balance_args is None:
                balance_args = []
            else:
                balance_args = shlex.split(balance_args)
            logger.debug('Balancing args: {}'.format(balance_args))

            for res in resolutions:
                uri = outfile + '::resolutions/' + str(res)
                if 'weight' in api.Cooler(uri).bins():
                    continue
                logger.info('Balancing zoom level with bin size {}'.format(res))
                result = runner.invoke(balance_cmd, args=[uri] + balance_args)
                if result.exit_code != 0:
                    raise result.exception

    else:
        n_zooms, zoom_levels = multires_aggregate(cool_uri, outfile, nproc,
            chunksize, lock=lock)

        if balance:
            runner = CliRunner()

            if balance_args is None:
                balance_args = []
            else:
                balance_args = shlex.split(balance_args)
            logger.debug('Balancing args: {}'.format(balance_args))

            for level, res in reversed(list(zoom_levels.items())):
                uri = outfile + '::' + str(level)
                if level == str(n_zooms):
                    if 'weight' in api.Cooler(uri).bins():
                        continue
                logger.info(
                    'Balancing zoom level {}, bin size {}'.format(level, res))
                result = runner.invoke(balance_cmd, args=[uri] + balance_args)
                if result.exit_code != 0:
                    raise result.exception
