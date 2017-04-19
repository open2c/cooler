# -*- coding: utf-8 -*-
from __future__ import division, print_function
import sys

import click
from . import cli
from .. import util


@cli.command()
@click.option(
    "--out", "-o",
    help="Output file (defaults to stdout)")
@click.argument(
    "chromsizes",
    type=str,
    metavar="CHROMSIZES_PATH")
@click.argument(
    "binsize",
    type=int,
    metavar="BINSIZE")
def makebins(chromsizes, binsize, out):
    """
    Generate fixed-width genomic bins.
    Output a genome segmentation at a fixed resolution as a BED file.

    CHROMSIZES_PATH : UCSC-like chromsizes file, with chromosomes in desired
    order.

    BINSIZE : Resolution (bin size) in base pairs <int>.

    """
    chromsizes = util.read_chromsizes(chromsizes, all_names=True)
    bins = util.binnify(chromsizes, binsize)

    # Write output
    try:
        if out is None:
            f = sys.stdout
        else:
            f = open(out, 'wt')
        bins.to_csv(f, sep='\t', index=False, header=False)
    except (IOError, OSError) as e:
        if e.errno == 32:  # broken pipe
            try:
                f.close()
            except OSError:
                pass
        else:
            raise
    else:
        f.close()
