# -*- coding: utf-8 -*-
from __future__ import division, print_function
import os.path as op
import json
import sys

import numpy as np
from pandas.api.types import is_float_dtype
import pandas as pd
import h5py

import click
from . import cli, logger
from .. import util
from ..io import (
    parse_cooler_uri, create_from_unsorted, sanitize_records, sanitize_pixels
)


# TODO: support dense text or memmapped binary/npy?
def _parse_field_params(args):
    extra_fields = []
    bad_param = False
    for arg in args:

        parts = arg.split(',')
        if len(parts) == 1 or len(parts) > 3:
            bad_param = True
        else:
            name = parts[0]
            try:
                number = int(parts[1]) - 1
            except ValueError:
                bad_param = True

            if number < 0:
                raise click.BadParameter(
                    "Field numbers are assumed to be 1-based.")

            if len(parts) == 3:
                dtype = np.dtype(parts[2])
            else:
                dtype = None

        if bad_param:
            raise click.BadParameter(
                "Expected '--field {{name}},{{number}}' "
                "or '--field {{name}},{{number}},{{dtype}}'; "
                "got '{}'".format(arg))
        extra_fields.append((name, number, dtype))

    return extra_fields


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
    type=click.Path(exists=True, allow_dash=True),
    metavar="PIXELS_PATH")
@click.argument(
    "cool_path",
    metavar="COOL_PATH")
@click.option(
    "--format", "-f",
    help="'coo' refers to a tab-delimited sparse triplet file "
         "(bin1, bin2, count). "
         "'bg2' refers to a 2D bedGraph-like file "
         "(chrom1, start1, end1, chrom2, start2, end2, count).",
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
    help="Add supplemental value fields or override default field numbers for "
         "the specified format. Specify as '<name>,<number>' or as "
         "'<name>,<number>,<dtype>' to enforce a dtype other than `float` or "
         "the default for a standard column. Field numbers are 1-based. "
         "Repeat the `--field` option for each additional field. "
         "[Changed in v0.7.7: use a comma separator, rather than a space.]",
    type=str,
    multiple=True)
@click.option(
    "--chunksize", "-c",
    help="Size (in number of lines/records) of data chunks to read and process "
         "from the input file at a time. These chunks will be saved as "
         "temporary partial Coolers and merged at the end. Also specifies the "
         "size of the buffer during the merge step.",
    type=int,
    default=int(20e6))
@click.option(
    "--count-as-float",
    is_flag=True,
    default=False,
    help="Store the 'count' column as floating point values instead of as "
         "integers. Can also be specified using the `--field` option.")
@click.option(
    "--one-based",
    is_flag=True,
    default=False,
    help="Pass this flag if the bin IDs listed in a COO file are one-based " 
         "instead of zero-based.")
@click.option(
    "--comment-char",
    type=str,
    default='#',
    show_default=True,
    help="Comment character that indicates lines to ignore.")
@click.option(
    "--tril-action",
    type=click.Choice(['reflect', 'drop']),
    default='reflect',
    show_default=True,
    help="How to handle lower triangle pixels. " 
         "'reflect': make lower triangle pixels upper triangular. "
         "Use this if your input data comes only from a unique half of a "
         "symmetric matrix (but may not respect the specified chromosome order)."
         "'drop': discard all lower triangle pixels. Use this if your input "
         "data is derived from a complete symmetric matrix.")
def load(bins_path, pixels_path, cool_path, format, metadata, assembly,
         chunksize, field, count_as_float, one_based, comment_char, tril_action):
    """
    Load a pre-binned contact matrix into a COOL file.

    \b
    Two input format options (tab-delimited):

    * COO: COO-rdinate sparse matrix format (a.k.a. ijv triple). 3 columns.

    \b
    - columns: "bin1_id, bin2_id, count",

    * BG2: 2D version of the bedGraph format. 7 columns.

    \b
    - columns: "chrom1, start1, end1, chrom2, start2, end2, count"

    Input pixel file may be compressed.

    **New in v0.7.7: Input files no longer need to be sorted or indexed!**

    Example:

    \b
    cooler load -f bg2 <chrom.sizes>:<binsize> in.bg2.gz out.cool

    \b\bArguments:

    BINS_PATH : One of the following

        <TEXT:INTEGER> : 1. Path to a chromsizes file, 2. Bin size in bp
        <TEXT> : Path to BED file defining the genomic bin segmentation.

    PIXELS_PATH : Text file containing nonzero pixel values. May be gzipped.
                  Pass '-' to use stdin.

    COOL_PATH : Output COOL file path

    """
    chromsizes, bins = _parse_bins(bins_path)

    # User-supplied JSON file
    if metadata is not None:
        with open(metadata, 'r') as f:
            metadata = json.load(f)

    output_field_names = ['bin1_id', 'bin2_id', 'count']
    output_field_dtypes = {
        'bin1_id': int,
        'bin2_id': int,
        'count': float if count_as_float else int,
    }

    if format == 'bg2':
        input_field_names = [
            'chrom1', 'start1', 'end1', 
            'chrom2', 'start2', 'end2', 
            'count'
        ]
        input_field_dtypes = {
            'chrom1': str, 'start1': int, 'end1': int,
            'chrom2': str, 'start2': int, 'end2': int,
            'count': float if count_as_float else int,
        }
        input_field_numbers = {
            'chrom1': 0, 'start1': 1, 'end1': 2,
            'chrom2': 3, 'start2': 4, 'end2': 5,
            'count': 6,
        }
        pipeline = sanitize_records(bins, 
            schema='bg2', 
            is_one_based=one_based,
            tril_action=tril_action)

    elif format == 'coo':
        input_field_names = [
            'bin1_id', 'bin2_id', 'count'
        ]
        input_field_dtypes = {
            'bin1_id': int, 
            'bin2_id': int,
            'count': float if count_as_float else int,
        }
        input_field_numbers = {
            'bin1_id': 0, 
            'bin2_id': 1, 
            'count': 2,
        }
        pipeline = sanitize_pixels(bins, 
            is_one_based=one_based,
            tril_action=tril_action)

    # include any additional value columns
    if len(field):
        extra_fields = _parse_field_params(field)
        for name, number, dtype in extra_fields:
            if name == 'count' and count_as_float and not is_float_dtype(dtype):
                raise ValueError(
                    "Mismatch between --count-as-float and 'count' dtype "
                    "'{}' provided via the --field option".format(dtype))

            if name not in input_field_names:
                input_field_names.append(name)
                output_field_names.append(name)
            
            input_field_numbers[name] = number

            if dtype is not None:
                input_field_dtypes[name] = dtype
                output_field_dtypes[name] = dtype

    if pixels_path == '-':
        f_in = sys.stdin
    else:
        f_in = pixels_path

    reader = pd.read_table(
        f_in, 
        usecols=[input_field_numbers[name] for name in input_field_names],
        names=input_field_names,
        dtype=input_field_dtypes,
        comment=comment_char,
        iterator=True,
        chunksize=chunksize)

    logger.info('fields: {}'.format(input_field_numbers))
    logger.info('dtypes: {}'.format(input_field_dtypes))

    create_from_unsorted(
        cool_path, 
        bins, 
        map(pipeline, reader), 
        columns=output_field_names,
        dtypes=output_field_dtypes,
        metadata=metadata, 
        assembly=assembly,
        mergebuf=chunksize
    )


