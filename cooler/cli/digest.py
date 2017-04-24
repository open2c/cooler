# -*- coding: utf-8 -*-
from __future__ import division, print_function
import os.path as op
import glob
import sys

import click
from . import cli
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
def digest(chromsizes, fasta, enzyme, out):
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

    # Write output
    try:
        if out is None:
            f = sys.stdout
        else:
            f = open(out, 'wt')
        frags.to_csv(f, sep='\t', index=False, header=False)
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
