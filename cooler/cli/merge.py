# -*- coding: utf-8 -*-
from __future__ import division, print_function

from ._util import parse_field_param
from . import cli, get_logger
import click

from ..reduce import merge_coolers

@cli.command()
@click.argument(
    "out_path",
    type=click.Path(exists=False))
@click.argument(
    "in_paths",
    nargs=-1,
    type=click.Path(exists=False))
@click.option(
    "--chunksize", "-c",
    help="Size of the merge buffer in number of pixel table rows.",
    type=int,
    default=int(20e6),
    show_default=True)
@click.option(
    "--field",
    help="Specify the names of value columns to merge as '<name>'. "
         "Repeat the `--field` option for each one. "
         "Use '<name>,dtype=<dtype>' to specify the dtype. Include "
         "',agg=<agg>' to specify an aggregation function different from 'sum'.",
    type=str,
    multiple=True)
def merge(out_path, in_paths, chunksize, field):
    """
    Merge multiple coolers with identical axes.

    OUT_PATH : Output file path or URI.

    IN_PATHS : Input file paths or URIs of coolers to merge.

    **Notes**

    Data columns merged:

        pixels/bin1_id, pixels/bin2_id, pixels/<value columns>

    Data columns preserved:

        chroms/name, chroms/length
        bins/chrom, bins/start, bins/end

    Additional columns in the the input files are not transferred to the output.

    """
    logger = get_logger(__name__)

    if len(field):
        field_specifiers = [
            parse_field_param(arg, includes_colnum=False) for arg in field
        ]
        columns, _, dtypes, agg = zip(*field_specifiers)
        dtypes = {col: dt for col, dt in zip(columns, dtypes) if dt is not None}
        agg = {col: f for col, f in zip(columns, agg) if f is not None}
    else:
        # If no other fields are given, 'count' is implicitly chosen.
        # Default aggregation. Dtype will be inferred.
        columns, dtypes, agg = ['count'], None, None

    merge_coolers(out_path, in_paths, mergebuf=chunksize, columns=columns,
           dtypes=dtypes, agg=agg)

