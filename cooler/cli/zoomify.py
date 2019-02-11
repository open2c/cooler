# -*- coding: utf-8 -*-
from __future__ import division, print_function
import os.path as op
import shlex

from ._util import parse_field_param
from . import cli, get_logger
from click.testing import CliRunner
import click

from ..reduce import (
    legacy_zoomify,
    zoomify_cooler,
    get_quadtree_depth,
    HIGLASS_TILE_DIM
)
from ..util import parse_cooler_uri
from ..create import create
from ..tools import lock
from .. import api


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
    '--resolutions', '-r',
    help="Comma-separated list of target resolutions.")
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
    '--base-uri', '-i',
    help="One or more additional base coolers to aggregate from, if needed.",
    multiple=True)
@click.option(
    '--out', '-o',
    help="Output file or URI")
@click.option(
    "--field",
    help="Specify the names of value columns to merge as '<name>'. "
         "Repeat the `--field` option for each one. "
         "Use '<name>:dtype=<dtype>' to specify the dtype. Include "
         "',agg=<agg>' to specify an aggregation function different from 'sum'.",
    type=str,
    multiple=True)
@click.option(
    '--legacy',
    help="Use the legacy layout of integer-labeled zoom levels.",
    is_flag=True,
    default=False)
def zoomify(cool_uri, nproc, chunksize, resolutions, balance, balance_args,
            field, legacy, base_uri, out):
    """
    Generate a multi-resolution cooler file by coarsening.

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


    if legacy:
        n_zooms, zoom_levels = legacy_zoomify(
            cool_uri, outfile, nproc, chunksize, lock=lock)

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
    else:
        if resolutions is not None:
            resolutions = [int(s.strip()) for s in resolutions.split(',')]
        else:
            clr = api.Cooler(cool_uri)
            n_zooms = get_quadtree_depth(clr.chromsizes, clr.binsize, HIGLASS_TILE_DIM)
            resolutions = [clr.binsize * 2**i for i in range(n_zooms)]

        if len(field):
            field_specifiers = [
                parse_field_param(arg, includes_colnum=False) for arg in field
            ]
            columns, _, dtypes, agg = zip(*field_specifiers)
            columns = list(columns)
            dtypes = {col: dt for col, dt in zip(columns, dtypes) if dt is not None}
            agg = {col: f for col, f in zip(columns, agg) if f is not None}
        else:
            # If no other fields are given, 'count' is implicitly chosen.
            # Default aggregation. Dtype will be inferred.
            columns, dtypes, agg = ['count'], None, None

        zoomify_cooler(
            [cool_uri] + list(base_uri),
            outfile,
            resolutions,
            chunksize,
            nproc=nproc,
            lock=lock,
            columns=columns,
            dtypes=dtypes,
            agg=agg)

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
