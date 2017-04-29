# -*- coding: utf-8 -*-
from __future__ import division, print_function
import os.path as op
import json
import sys

import numpy as np
import pandas as pd
import h5py

import click
from . import cli
from .. import util
from ..io import create, parse_cooler_uri, SparseLoader, BedGraph2DLoader


# TODO: support dense text or memmapped binary/npy?

def _validate_fieldnum(ctx, param, value):
    if value <= 0:
        raise click.BadParameter('Field numbers are one-based')
    return value


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
    help="'coo' refers to a tab-delimited sparse triple file (bin1, bin2, count). "
         "'bg2' refers to a 2D bedGraph-like file (chrom1, start1, end1, chrom2, start2, end2, count).",
    type=click.Choice(['coo', 'bg2']),
    required=True)
@click.option(
    "--metadata",
    help="Path to JSON file containing user metadata.")
@click.option(
    "--assembly",
    help="Name of genome assembly (e.g. hg19, mm10)")
@click.option(
    "--field",
    help="Add supplemental pixel fields or override default pixel field numbers. "
         "Specify as '<name> <number>'. Field numbers are 1-based. "
         "Supplemental data columns are assumed to be floating point. "
         "Repeat for each additional field.",
    nargs=2,
    type=(str, int),
    multiple=True)
@click.option(
    "--chunksize", "-c",
    type=int,
    default=int(10e6))
def load(bins_path, pixels_path, cool_path, format, metadata, assembly, chunksize, field):
    """
    Load a contact matrix.
    Load a sparse-formatted text dump of a contact matrix into a COOL file.

    \b
    Two input format options (tab-delimited):

    * COO: COO-rdinate matrix format (i.e. ijv triple). 3 columns.

    \b
    - columns: "bin1_id, bin2_id, count",
    - lexicographically sorted by bin1_id, bin2_id
    - optionally compressed

    * BG2: 2D version of the bedGraph format. 7 columns.

    \b
    - columns: "chrom1, start1, end1, chrom2, start2, end2, count"
    - sorted by chrom1, chrom2, start1, start2
    - bgzip compressed and indexed with Pairix (see cooler csort)

    Example:

    \b
    cooler csort -c1 1 -p1 2 -c2 4 -p2 5 <in.bg2> <chrom.sizes>
    cooler load -f bg2 <chrom.sizes>:<binsize> <in.bg2.srt.gz> <out.cool>

    \b\bArguments:

    BINS_PATH : One of the following

        <TEXT:INTEGER> : 1. Path to a chromsizes file, 2. Bin size in bp
        <TEXT> : Path to BED file defining the genomic bin segmentation.

    PIXELS_PATH : Text file containing nonzero pixel values. May be gzipped.

    COOL_PATH : Output COOL file path

    """
    chromsizes, bins = _parse_bins(bins_path)

    # User-supplied JSON file
    if metadata is not None:
        with open(metadata, 'r') as f:
            metadata = json.load(f)

    # Set up the appropriate binned contacts loader
    if field is not None:
        if not all(v > 0 for k, v in field):
            raise click.BadParameter("Field numbers are assumed to be 1-based.")
        field_numbers = {k: v-1 for k, v in field}
        field_dtypes = {k: float for k, v in field}
    else:
        field_numbers = None
        field_dtypes = None

    if format == 'bg2':
        binner = BedGraph2DLoader(pixels_path, chromsizes, bins, 
                                  field_numbers, field_dtypes)
    elif format == 'coo':
        binner = SparseLoader(pixels_path, bins, chunksize, 
                              field_numbers, field_dtypes)

    create(cool_path, bins, binner, metadata, assembly, dtypes=field_dtypes)
