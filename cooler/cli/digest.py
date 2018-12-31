# -*- coding: utf-8 -*-
from __future__ import division, print_function
import os.path as op
import glob
import sys

from ._util import exit_on_broken_pipe
from . import cli
import click

from .. import util


@cli.command()
@click.argument(
    "chromsizes",
    metavar="CHROMSIZES_PATH")
@click.argument(
    "fasta",
    metavar="FASTA_PATH")
@click.argument(
    "enzyme",
    metavar="ENZYME")
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
def digest(chromsizes, fasta, enzyme, out, header, rel_ids):
    """
    Generate fragment-delimited genomic bins.

    Output a genome segmentation of restriction fragments as a BED file.

    CHROMSIZES_PATH : UCSC-like chromsizes file, with chromosomes in desired
    order.

    FASTA_PATH : Genome assembly FASTA file or folder containing FASTA files
    (uncompressed).

    ENZYME : Name of restriction enzyme

    """
    chromsizes = util.read_chromsizes(chromsizes, all_names=True)
    chroms = list(chromsizes.keys())

    # Load sequences
    if op.isdir(fasta):
        filepaths = glob.glob(op.join(fasta, '*.fa'))
        filepaths.extend(glob.glob(op.join(fasta, '*.fasta')))
    else:
        filepaths = [fasta]
    fasta_records = util.load_fasta(chroms, *filepaths)

    # Digest sequences
    frags = util.digest(fasta_records, enzyme)

    if rel_ids is not None:
        frags['id'] = frags.groupby('chrom').cumcount()
        if int(rel_ids) == 1:
            frags['id'] += 1

    # Write output
    if out is None:
        f = sys.stdout
    else:
        f = open(out, 'wt')

    if header:
        frags[0:0].to_csv(f, sep='\t', index=False, header=True)
    frags.to_csv(f, sep='\t', index=False, header=False)
    f.flush()
