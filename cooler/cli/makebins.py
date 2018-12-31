# -*- coding: utf-8 -*-
from __future__ import division, print_function
import sys

from ._util import exit_on_broken_pipe
from . import cli
import click

from .. import util


@cli.command()
@click.argument(
    "chromsizes",
    type=str,
    metavar="CHROMSIZES_PATH")
@click.argument(
    "binsize",
    type=int,
    metavar="BINSIZE")
@click.option(
    "--out", "-o",
    help="Output file (defaults to stdout)")
@click.option(
    "--header", "-H",
    help="Print the header of column names as the first row.",
    is_flag=True,
    default=False,
    show_default=True)
@click.option(
    "--rel-ids", "-i",
    type=click.Choice(['0', '1']),
    help="Include a column of relative bin IDs for each chromosome. "
         "Choose whether to report them as 0- or 1-based.")
@exit_on_broken_pipe(1)
def makebins(chromsizes, binsize, out, header, rel_ids):
    """
    Generate fixed-width genomic bins.

    Output a genome segmentation at a fixed resolution as a BED file.

    CHROMSIZES_PATH : UCSC-like chromsizes file, with chromosomes in desired
    order.

    BINSIZE : Resolution (bin size) in base pairs <int>.

    """
    chromsizes = util.read_chromsizes(chromsizes, all_names=True)
    bins = util.binnify(chromsizes, binsize)

    if rel_ids is not None:
        bins['id'] = bins.groupby('chrom').cumcount()
        if int(rel_ids) == 1:
            bins['id'] += 1

    # Write output
    if out is None:
        f = sys.stdout
    else:
        f = open(out, 'wt')

    if header:
        bins[0:0].to_csv(f, sep='\t', index=False, header=True)
    bins.to_csv(f, sep='\t', index=False, header=False)
    f.flush()
