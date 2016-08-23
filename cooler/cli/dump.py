# -*- coding: utf-8 -*-
from __future__ import division, print_function
import gzip
import sys

import numpy as np

import click
from . import cli
from ..api import Cooler, annotate


@cli.command()
@click.argument(
    "cooler_file",
    metavar="COOLER_PATH")
@click.option(
    "--bins", "-b",
    help="Dump the bin table.",
    is_flag=True,
    default=False)
@click.option(
    "--join",
    help="Print chromosome bin coordinates instead of bin IDs",
    is_flag=True,
    default=False)
@click.option(
    "--balanced",
    help="Apply balancing weights to data",
    is_flag=True,
    default=False)
@click.option(
    '--chunksize',
    help="Sets the amount of data loaded from disk at one time",
    type=int,
    default=int(1e6))
@click.option(
    "--out", "-o",
    help="Output text file")
def dump(cooler_file, output_bins, join, balanced, chunksize, out):
    """
    Dump a cooler file to a tsv file.

    COOLER_PATH : Cooler file.

    """
    c = Cooler(cooler_file)
    
    # output stream
    if out is None or out == '-':
        f = sys.stdout
    elif out.endswith('.gz'):
        f = gzip.open(out, 'wt')
    else:
        f = open(out, 'wt')

    # choose the source
    if output_bins:
        selector = c.bins()
        n = c.info['nbins']
    else:
        selector = c.pixels()
        n = c.info['nnz']
        bins = c.bins()[:]  # load all the bins
        if balanced and 'weight' not in bins.columns:
            print('Balancing weights not found', file=sys.stderr)
            sys.exit(1)

    # write in chunks
    edges = np.arange(0, n+chunksize, chunksize)
    edges[-1] = n

    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = selector[lo:hi]

        if not output_bins:

            if join:
                sel = annotate(sel, bins[['chrom', 'start', 'end']])

            if balanced:
                df = annotate(sel, bins[['weight']])
                sel['balanced'] = df['weight1'] * df['weight2'] * sel['count']

        try:
            sel.to_csv(f, sep='\t', index=False, header=True)

        except OSError as e:
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
