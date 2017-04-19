# -*- coding: utf-8 -*-
from __future__ import division, print_function
import json
import sys

import numpy as np
import pandas as pd
import h5py

import click
from . import cli
from ..io import create, SparseLoader, BedGraph2DLoader


# TODO: support positional and block sorted pixel data


def _parse_bins(arg):
    # Provided chromsizes and binsize
    if ":" in arg:
        chromsizes_file, binsize = arg.split(":")
        if not op.exists(chromsizes_file):
            raise ValueError('File "{}" not found'.format(chromsizes_file))
        try:
            binsize = int(binsize)
        except ValueError:
            raise ValueError(
                'Expected integer binsize argument (bp), got "{}"'.format(binsize))
        chromsizes = util.read_chromsizes(chromsizes_file, all_names=True)
        bins = util.binnify(chromsizes, binsize)

    # Provided bins
    elif op.exists(arg):
        try:
            bins = pd.read_csv(
                arg,
                sep='\t',
                names=['chrom', 'start', 'end'],
                usecols=[0, 1, 2],
                dtype={'chrom': str})
        except pd.parser.CParserError as e:
            raise ValueError(
                'Failed to parse bins file "{}": {}'.format(arg, str(e)))

        chromtable = (
            bins.drop_duplicates(['chrom'], keep='last')[['chrom', 'end']]
                .reset_index(drop=True)
                .rename(columns={'chrom': 'name', 'end': 'length'})
        )
        chroms, lengths = list(chromtable['name']), list(chromtable['length'])
        chromsizes = pd.Series(index=chroms, data=lengths)
        
    else:
        raise ValueError(
            'Expected BINS to be either <Path to bins file> or '
            '<Path to chromsizes file>:<binsize in bp>.')

    return chromsizes, bins


@cli.command()
@click.argument(
    "bins_path",
    metavar="BINS_PATH")
@click.argument(
    "pixels_path",
    metavar="PIXELS_PATH")
@click.argument(
    "cool_path",
    metavar="COOL_PATH")
@click.option(
    "--format", "-f",
    help="COO refers to a tab-delimited sparse triple file (bin1, bin2, count). "
         "BG2 refers to a 2D bedGraph-like file (chrom1, start1, end1, chrom2, start2, end2, count).",
    type=click.Choice(['BG2', 'COO']),
    default='COO',
    show_default=True)
@click.option(
    "--chunksize", "-c",
    type=int,
    default=int(10e6))
@click.option(
    "--metadata",
    help="Path to JSON file containing user metadata.")
@click.option(
    "--assembly",
    help="Name of genome assembly (e.g. hg19, mm10)")
def load(bins_path, pixels_path, cool_path, format, metadata, assembly, chunksize):
    """
    Load a contact matrix.
    Load a sparse-formatted text dump of a contact matrix into a COOL file.

    \b
    Two input format options (tab-delimited):

    * COO: COO-rdinate matrix format (i.e. ijv triple). 3 columns.

    \b
    - columns: "bin1_id, bin2_id, count",
    - lexicographically sorted by bin1_id, bin2_id

    * BG2: 2D version of the bedGraph format. 7 columns.

    \b
    - columns: "chrom1, start1, end1, chrom2, start2, end2, count"
    - each record has chrom1 <= chrom2, with respect to the chromosome order in BINS
    - sorted by chrom1, start1, chrom2, start2,  with respect to the chromosome order in BINS

    \b\bArguments:

    BINS_PATH : BED-like file containing genomic bin segmentation

    PIXELS_PATH : Text file containing nonzero pixel values. May be gzipped.

    COOL_PATH : Output COOL file path

    """
    chromsizes, bins = _parse_bins(bins)

    # User-supplied JSON file
    if metadata is not None:
        with open(metadata, 'r') as f:
            metadata = json.load(f)

    # Load the binned contacts
    if format == 'BG2':
        iterator = BedGraph2DLoader(bins, pixels_path, chunksize)
    else:
        iterator = SparseLoader(pixels_path, chunksize)
    create(cool_path, chromsizes, bins, iterator, metadata, assembly)
