# -*- coding: utf-8 -*-
from __future__ import division, print_function
import gzip
import sys

import numpy as np
import pandas as pd

from ._util import DelimitedTuple, exit_on_broken_pipe
from . import cli
import click

from .. import api


@cli.command()
@click.argument(
    "cool_uri",
    metavar="COOL_PATH")
@click.option(
    "--table", "-t",
    help="Which table to dump. Choosing 'chroms' or 'bins' will cause all "
         "pixel-related options to be ignored. "
         "Note that dumping 'pixels' will only provide data for the upper "
         "triangle of the contact matrix. ",
    type=click.Choice(['chroms', 'bins', 'pixels']),
    default='pixels',
    show_default=True)
@click.option(
    "--columns", "-c",
    help="Restrict output to a subset of columns, provided as a "
         "comma-separated list.",
    type=DelimitedTuple(sep=','))
@click.option(
    "--header", "-H",
    help="Print the header of column names as the first row.",
    is_flag=True,
    default=False,
    show_default=True)
@click.option(
    "--range", "-r",
    help="The coordinates of a genomic region shown along the row dimension, "
         "in UCSC notation. (Example: chr1:10,000,000-11,000,000). "
         "If omitted, the entire contact matrix is printed.",
    type=str)
@click.option(
    "--range2", "-r2",
    type=str,
    help="The coordinates of a genomic region shown along the column dimension. "
         "If omitted, the column range is the same as the row range.")
@click.option(
    "--balanced/--no-balance", "-b",
    help="Apply balancing weights to data. This will print an extra column "
         "called `balanced`",
    is_flag=True,
    default=False,
    show_default=True)
@click.option(
    "--join",
    help="Print the full chromosome bin coordinates instead of bin IDs. "
         "This will replace the `bin1_id` column with `chrom1`, `start1`, and "
         "`end1`, and the `bin2_id` column with `chrom2`, `start2` and `end2`.",
    is_flag=True,
    default=False,
    show_default=True)
@click.option(
    "--annotate",
    help="Join additional columns from the bin table against the pixels. "
         "Provide a comma separated list of column names (no spaces). "
         "The merged columns will be suffixed by '1' and '2' accordingly.",
    type=DelimitedTuple(sep=','))
@click.option(
    "--na-rep",
    help="Missing data representation. Default is empty ''.",
    default='')
@click.option(
    "--float-format",
    help="Format string for floating point numbers (e.g. '.12g', '03.2f').",
    default='g',
    show_default=True)
@click.option(
    '--chunksize', "-k",
    help="Sets the amount of pixel data loaded from disk at one time. "
         "Can affect the performance of joins on high resolution datasets. "
         "Default is to load as many rows as there are bins.",
    type=int)
@click.option(
    "--out", "-o",
    help="Output text file If .gz extension is detected, file is written "
         "using zlib. Default behavior is to stream to stdout.")
# duplex (not use unique values for symmetric matrices)
@exit_on_broken_pipe(1)
def dump(cool_uri, table, columns, header, range, range2, balanced, join,
         annotate, na_rep, float_format, chunksize, out):
    """
    Dump a cooler's data to a text stream.

    COOL_PATH : Path to COOL file or cooler URI.

    """
    c = api.Cooler(cool_uri)

    # output stream
    if out is None or out == '-':
        f = sys.stdout
    elif out.endswith('.gz'):
        f = gzip.open(out, 'wt')
    else:
        f = open(out, 'wt')

    # choose the source
    if table == 'chroms':
        selector = c.chroms()
        n = c.info['nchroms']
        chunksize = n
    elif table == 'bins':
        selector = c.bins()
        n = c.info['nbins']
        chunksize = n
    else:
        # load all the bins
        bins = c.bins()[:]
        if balanced and 'weight' not in bins.columns:
            print('Balancing weights not found', file=sys.stderr)
            sys.exit(1)

        if range:
            selector = (c.matrix(as_pixels=True, balance=balanced)
                         .fetch(range, range2))
            n = len(selector)
        else:
            selector = c.pixels()
            n = c.info['nnz']

        if chunksize is None:
            chunksize = len(bins)

    if columns is not None:
        selector = selector[list(columns)]

    # write in chunks
    edges = np.arange(0, n+chunksize, chunksize)
    edges[-1] = n

    first = True
    if float_format is not None:
        float_format = '%' + float_format

    for lo, hi in zip(edges[:-1], edges[1:]):

        sel = selector[lo:hi]

        if range:
            sel = sel.copy()  # suppress pandas warning

        if table == 'pixels':

            if annotate is not None:
                extra_fields = list(annotate)
                try:
                    extra_cols = bins[extra_fields]
                except KeyError as e:
                    print('Column not found:\n {}'.format(e))
                    sys.exit(1)
                extra = api.annotate(sel[['bin1_id', 'bin2_id']], extra_cols)

            if balanced:
                df = api.annotate(sel, bins[['weight']])
                sel['balanced'] = df['weight1'] * df['weight2'] * sel['count']

            if join:
                sel = api.annotate(sel, bins[['chrom', 'start', 'end']])

            if annotate is not None:
                sel = pd.concat([sel, extra], axis=1)

        if first:
            if header:
                sel[0:0].to_csv(
                    f,
                    sep='\t',
                    index=False,
                    header=True,
                    float_format=float_format)
            first = False

        sel.to_csv(
            f,
            sep='\t',
            index=False,
            header=False,
            float_format=float_format,
            na_rep=na_rep)

    else:
        f.flush()
