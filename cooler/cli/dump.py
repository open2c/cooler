# -*- coding: utf-8 -*-
from __future__ import division, print_function
import gzip
import sys

import numpy as np
import pandas as pd

import click
from . import cli
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
    "--header",
    help="Print the header of column names as the first row.",
    is_flag=True,
    default=False,
    show_default=True)
@click.option(
    '--chunksize', "-k",
    help="Sets the amount of pixel data loaded from disk at one time. "
         "Can affect the performance of joins on high resolution datasets. "
         "Default is to load as many rows as there are bins.",
    type=int)
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
    default='')
@click.option(
    "--out", "-o",
    help="Output text file If .gz extension is detected, file is written "
         "using zlib. Default behavior is to stream to stdout.")
def dump(cool_uri, table, chunksize, range, range2, join,
         annotate, balanced, header, out):
    """
    Dump a contact matrix.
    Print the contents of a COOL file to tab-delimited text.

    COOL_PATH : Path to COOL file or Cooler URI.

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

    # write in chunks
    edges = np.arange(0, n+chunksize, chunksize)
    edges[-1] = n

    first = True

    for lo, hi in zip(edges[:-1], edges[1:]):

        sel = selector[lo:hi]

        if range:
            sel = sel.copy()  # suppress pandas warning

        if table == 'pixels':

            if annotate:
                extra_fields = annotate.split(',')
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

            if annotate:
                sel = pd.concat([sel, extra], axis=1)

        if first:
            if header:
                sel[0:0].to_csv(
                    f, sep='\t', index=False, header=True, float_format='%g')
            first = False

        try:
            sel.to_csv(
                f, sep='\t', index=False, header=False, float_format='%g')

        except (IOError, OSError) as e:
            if e.errno == 32:  # broken pipe
                try:
                    f.close()
                except OSError:
                    pass
                break
            else:
                raise
    else:
        f.close()
