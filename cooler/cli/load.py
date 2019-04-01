# -*- coding: utf-8 -*-
from __future__ import division, print_function
import os.path as op
import json
import sys

import numpy as np
from pandas.api.types import is_float_dtype
import pandas as pd
import h5py

from ._util import parse_bins, parse_field_param
from . import cli, get_logger
import click

from .. import util
from ..util import parse_cooler_uri
from ..create import (
    create_from_unordered, sanitize_records, sanitize_pixels,
    BIN_DTYPE, COUNT_DTYPE
)


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
         "the specified format. "
         "Specify quantitative input fields to aggregate into value columns "
         "using the syntax ``--field <field-name>=<field-number>``. "
         "Optionally, append ``:`` followed by ``dtype=<dtype>`` to specify "
         "the data type (e.g. float). Field numbers are 1-based. "
         "Repeat the ``--field`` option for each additional field.",
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
    "--no-symmetric-upper", "-N",
    help="Create a complete square matrix without implicit symmetry. "
         "This allows for distinct upper- and lower-triangle values",
    is_flag=True,
    default=False)
@click.option(
    "--input-copy-status",
    type=click.Choice(['unique', 'duplex']),
    default='unique',
    help="Copy status of input data when using symmetric-upper storage. | "
         "`unique`: Incoming data comes from a unique half of a symmetric "
         "matrix, regardless of how element coordinates are ordered. "
         "Execution will be aborted if duplicates are detected. "
         "`duplex`: Incoming data contains upper- and lower-triangle duplicates. "
         "All lower-triangle input elements will be discarded! | "
         "If you wish to treat lower- and upper-triangle input data as "
         "distinct, use the ``--no-symmetric-upper`` option instead. ",
    show_default=True)
@click.option(
    "--storage-options",
    help="Options to modify the data filter pipeline. Provide as a "
         "comma-separated list of key-value pairs of the form 'k1=v1,k2=v2,...'. "
         "See http://docs.h5py.org/en/stable/high/dataset.html#filter-pipeline "
         "for more details.")
def load(bins_path, pixels_path, cool_path, format, metadata, assembly,
         chunksize, field, count_as_float, one_based, comment_char,
         input_copy_status, no_symmetric_upper, storage_options, **kwargs):
    """
    Create a cooler from a pre-binned matrix.

    BINS_PATH : One of the following

        <TEXT:INTEGER> : 1. Path to a chromsizes file, 2. Bin size in bp

        <TEXT> : Path to BED file defining the genomic bin segmentation.

    PIXELS_PATH : Text file containing nonzero pixel values. May be gzipped.
    Pass '-' to use stdin.

    COOL_PATH : Output COOL file path or URI.

    **Notes**

    Two input format options (tab-delimited).
    Input pixel file may be compressed.

    COO: COO-rdinate sparse matrix format (a.k.a. ijv triple).
    3 columns: "bin1_id, bin2_id, count",

    BG2: 2D version of the bedGraph format.
    7 columns: "chrom1, start1, end1, chrom2, start2, end2, count"

    **Examples**

    cooler load -f bg2 <chrom.sizes>:<binsize> in.bg2.gz out.cool

    """
    logger = get_logger(__name__)
    chromsizes, bins = parse_bins(bins_path)

    symmetric_upper = not no_symmetric_upper
    tril_action = None
    if symmetric_upper:
        if input_copy_status == 'unique':
            tril_action = 'reflect'
        elif input_copy_status == 'duplex':
            tril_action = 'drop'

    # User-supplied JSON file
    if metadata is not None:
        with open(metadata, 'r') as f:
            metadata = json.load(f)

    # Initialize the output schema. We don't include 'count' yet.
    output_field_names = ['bin1_id', 'bin2_id']
    output_field_dtypes = {
        'bin1_id': BIN_DTYPE,
        'bin2_id': BIN_DTYPE,
        'count': COUNT_DTYPE
    }

    # Initialize the input schema and create the input santizer.
    if format == 'bg2':
        input_field_names = [
            'chrom1',
            'start1',
            'end1',
            'chrom2',
            'start2',
            'end2',
            # We don't include 'count' yet.
        ]
        input_field_dtypes = {
            'chrom1': str,
            'start1': int,
            'end1': int,
            'chrom2': str,
            'start2': int,
            'end2': int,
            'count': output_field_dtypes['count'],
        }
        input_field_numbers = {
            'chrom1': kwargs.get('chrom1', 0),
            'start1': kwargs.get('start1', 1),
            'end1':  kwargs.get('end1', 2),
            'chrom2': kwargs.get('chrom2', 3),
            'start2': kwargs.get('start2', 4),
            'end2':  kwargs.get('end2', 5),
            'count': 6,
        }
        pipeline = sanitize_records(bins,
            schema='bg2',
            is_one_based=one_based,
            tril_action=tril_action,
            sort=True)

    elif format == 'coo':
        input_field_names = [
            'bin1_id',
            'bin2_id',
            # We don't include 'count' yet.
        ]
        input_field_dtypes = {
            'bin1_id': int,
            'bin2_id': int,
            'count': output_field_dtypes['count'],
        }
        input_field_numbers = {
            'bin1_id': 0,
            'bin2_id': 1,
            'count': 2,
        }
        pipeline = sanitize_pixels(bins,
            is_one_based=one_based,
            tril_action=tril_action,
            sort=True)

    # Include input value columns
    if len(field):
        for arg in field:
            name, colnum, dtype, _ = parse_field_param(arg, includes_agg=False)

            # Special cases: omit field number to change standard dtypes.
            if colnum is None:
                if name in {'bin1_id', 'bin2_id'} and dtype is not None:
                    # No input field
                    output_field_dtypes[name] = dtype
                    continue
                elif name == 'count' and dtype is not None:
                    input_field_names.append('count')
                    output_field_names.append('count')
                    input_field_dtypes[name] = dtype
                    output_field_dtypes[name] = dtype
                    continue
                else:
                    raise click.BadParameter(
                        "A field number is required.", param_hint=arg)

            if name not in input_field_names:
                input_field_names.append(name)

            if name not in output_field_names:
                output_field_names.append(name)

            input_field_numbers[name] = colnum

            if dtype is not None:
                input_field_dtypes[name] = dtype
                output_field_dtypes[name] = dtype
    else:
        # If no other fields are given, 'count' is implicitly included.
        # Default dtype and field number are assumed.
        input_field_names.append('count')
        output_field_names.append('count')

    if 'count' in input_field_names and count_as_float:
        input_field_dtypes['count'] = np.float64
        output_field_dtypes['count'] = np.float64

    # Customize the HDF5 filters
    if storage_options is not None:
        h5opts = _parse_kv_list_param(storage_options)
        for key in h5opts:
            if isinstance(h5opts[key], list):
                h5opts[key] = tuple(h5opts[key])
    else:
        h5opts = None

    # Initialize the input stream
    if pixels_path == '-':
        f_in = sys.stdin
    else:
        f_in = pixels_path

    reader = pd.read_csv(
        f_in,
        sep='\t',
        usecols=[input_field_numbers[name] for name in input_field_names],
        names=input_field_names,
        dtype=input_field_dtypes,
        comment=comment_char,
        iterator=True,
        chunksize=chunksize)

    logger.info('fields: {}'.format(input_field_numbers))
    logger.info('dtypes: {}'.format(input_field_dtypes))
    logger.info('symmetric-upper: {}'.format(symmetric_upper))

    create_from_unordered(
        cool_path,
        bins,
        map(pipeline, reader),
        columns=output_field_names,
        dtypes=output_field_dtypes,
        metadata=metadata,
        assembly=assembly,
        mergebuf=chunksize,
        ensure_sorted=False,
        #boundscheck=True,
        #dupcheck=True,
        triucheck=True if symmetric_upper else False,
        symmetric_upper=symmetric_upper,
        h5opts=h5opts
    )
