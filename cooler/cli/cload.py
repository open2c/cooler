# -*- coding: utf-8 -*-
from __future__ import division, print_function
from multiprocess import Pool
from functools import partial
import os.path as op
import json
import six
import sys

from cytoolz import compose
import numpy as np
import pandas as pd
import h5py

from ._util import parse_bins, parse_kv_list_param, parse_field_param
from . import cli, get_logger
import click

from .. import util
from ..create import (
    create_cooler,
    sanitize_records, aggregate_records,
    TabixAggregator, HDF5Aggregator, PairixAggregator,
)


@cli.group()
def cload():
    """
    Create a cooler from genomic pairs and bins.

    Choose a subcommand based on the format of the input contact list.

    """
    pass


def register_subcommand(func):
    return (
        cload.command()(
        click.argument(
            "bins",
            type=str,
            metavar="BINS")(
        click.argument(
            "pairs_path",
            type=click.Path(exists=True, allow_dash=True),
            metavar="PAIRS_PATH")(
        click.argument(
            "cool_path",
            metavar="COOL_PATH")(
        click.option(
            "--metadata",
            help="Path to JSON file containing user metadata.")(
        click.option(
            "--assembly",
            help="Name of genome assembly (e.g. hg19, mm10)")(
        func))))))
    )


def add_arg_help(func):
    func.__doc__ = func.__doc__.format(
    """BINS : One of the following

        <TEXT:INTEGER> : 1. Path to a chromsizes file, 2. Bin size in bp

        <TEXT> : Path to BED file defining the genomic bin segmentation.

    PAIRS_PATH : Path to contacts (i.e. read pairs) file.

    COOL_PATH : Output COOL file path or URI.""")
    return func


@register_subcommand
@add_arg_help
@click.option(
    "--chunksize", "-c",
    help="Control the number of pixels handled by each worker process at a time.",
    type=int,
    default=int(100e6),
    show_default=True)
def hiclib(bins, pairs_path, cool_path, metadata, assembly, chunksize):
    """
    Bin a hiclib HDF5 contact list (frag) file.

    {}

    hiclib on BitBucket: <https://bitbucket.org/mirnylab/hiclib>.

    """

    chromsizes, bins = parse_bins(bins)

    if metadata is not None:
        with open(metadata, 'r') as f:
            metadata = json.load(f)

    with h5py.File(pairs_path, 'r') as h5pairs:
        iterator = HDF5Aggregator(h5pairs, chromsizes, bins, chunksize)
        create_cooler(
            cool_path, bins, iterator,
            metadata=metadata,
            assembly=assembly,
            ordered=True)


@register_subcommand
@add_arg_help
@click.option(
    "--nproc", "-p",
    help="Number of processes to split the work between.",
    type=int,
    default=8,
    show_default=True)
@click.option(
    "--chrom2", "-c2",
    help="chrom2 field number (one-based)",
    type=int,
    default=4)
@click.option(
    "--pos2", "-p2",
    help="pos2 field number (one-based)",
    type=int,
    default=5)
@click.option(
    "--zero-based", "-0",
    help="Positions are zero-based",
    is_flag=True,
    default=False,
    show_default=True)
@click.option(
    "--max-split", "-s",
    help="Divide the pairs from each chromosome into at most this many chunks. "
         "Smaller chromosomes will be split less frequently or not at all. "
         "Increase ths value if large chromosomes dominate the workload on "
         "multiple processors.",
    type=int,
    default=2,
    show_default=True)
def tabix(bins, pairs_path, cool_path, metadata, assembly, nproc, zero_based, max_split, **kwargs):
    """
    Bin a tabix-indexed contact list file.

    {}

    See also: 'cooler csort' to sort and index a contact list file

    Tabix manpage: <http://www.htslib.org/doc/tabix.html>.

    """
    logger = get_logger(__name__)
    chromsizes, bins = parse_bins(bins)

    if metadata is not None:
        with open(metadata, 'r') as f:
            metadata = json.load(f)

    try:
        if nproc > 1:
            pool = Pool(nproc)
            logger.info("Using {} cores".format(nproc))
            map = pool.imap
        else:
            map = six.moves.map
        opts = {}
        if 'chrom2' in kwargs:
            opts['C2'] = kwargs['chrom2'] - 1
        if 'pos2' in kwargs:
            opts['P2'] = kwargs['pos2'] - 1
        iterator = TabixAggregator(pairs_path, chromsizes, bins, map=map,
            is_one_based=(not zero_based), n_chunks=max_split, **opts)
        create_cooler(
            cool_path, bins, iterator,
            metadata=metadata,
            assembly=assembly,
            ordered=True)
    finally:
        if nproc > 1:
            pool.close()


@register_subcommand
@add_arg_help
@click.option(
    "--nproc", "-p",
    help="Number of processes to split the work between.",
    type=int,
    default=8,
    show_default=True)
@click.option(
    "--zero-based", "-0",
    help="Positions are zero-based",
    is_flag=True,
    default=False,
    show_default=True)
@click.option(
    "--max-split", "-s",
    help="Divide the pairs from each chromosome into at most this many chunks. "
         "Smaller chromosomes will be split less frequently or not at all. "
         "Increase ths value if large chromosomes dominate the workload on "
         "multiple processors.",
    type=int,
    default=2,
    show_default=True)
def pairix(bins, pairs_path, cool_path, metadata, assembly, nproc, zero_based, max_split):
    """
    Bin a pairix-indexed contact list file.

    {}

    See also: 'cooler csort' to sort and index a contact list file

    Pairix on GitHub: <https://github.com/4dn-dcic/pairix>.

    """
    logger = get_logger(__name__)
    chromsizes, bins = parse_bins(bins)

    if metadata is not None:
        with open(metadata, 'r') as f:
            metadata = json.load(f)

    try:
        if nproc > 1:
            pool = Pool(nproc)
            logger.info("Using {} cores".format(nproc))
            map = pool.imap
        else:
            map = six.moves.map
        iterator = PairixAggregator(pairs_path, chromsizes, bins, map=map,
            is_one_based=(not zero_based), n_chunks=max_split)
        create_cooler(
            cool_path, bins, iterator,
            metadata=metadata,
            assembly=assembly,
            ordered=True)
    finally:
        if nproc > 1:
            pool.close()


@register_subcommand
@add_arg_help
@click.option(
    "--chrom1", "-c1",
    help="chrom1 field number (one-based)",
    type=int,
    required=True)
@click.option(
    "--pos1", "-p1",
    help="pos1 field number (one-based)",
    type=int,
    required=True)
@click.option(
    "--chrom2", "-c2",
    help="chrom2 field number (one-based)",
    type=int,
    required=True)
@click.option(
    "--pos2", "-p2",
    help="pos2 field number (one-based)",
    type=int,
    required=True)
@click.option(
    "--chunksize",
    help="Number of input lines to load at a time",
    type=int,
    default=int(15e6))
@click.option(
    "--zero-based", "-0",
    help="Positions are zero-based",
    is_flag=True,
    default=False,
    show_default=True)
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
         "map, regardless of how the coordinates of a pair are ordered. "
         "`duplex`: Incoming data contains upper- and lower-triangle duplicates. "
         "All input records that map to the lower triangle will be discarded! | "
         "If you wish to treat lower- and upper-triangle input data as "
         "distinct, use the ``--no-symmetric-upper`` option. ",
    show_default=True)
@click.option(
    "--field",
    help="Specify quantitative input fields to aggregate into value columns "
         "using the syntax ``--field <field-name>=<field-number>``. "
         "Optionally, append ``:`` followed by ``dtype=<dtype>`` to specify "
         "the data type (e.g. float), and/or ``agg=<agg>`` to "
         "specify an aggregation function different from sum (e.g. mean). "
         "Field numbers are 1-based. Passing 'count' as the target name will "
         "override the default behavior of storing pair counts. "
         "Repeat the ``--field`` option for each additional field.",
    type=str,
    multiple=True)
# @click.option(
#     "--no-count",
#     help="Do not store the pair counts. Use this only if you use `--field` to "
#          "specify at least one input field for aggregation as an alternative.",
#     is_flag=True,
#     default=False)
@click.option(
    "--temp-dir",
    help="Create temporary files in a specified directory. Pass ``-`` to use "
         "the platform default temp dir.",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, allow_dash=True))
@click.option(
    "--no-delete-temp",
    help="Do not delete temporary files when finished.",
    is_flag=True,
    default=False)
@click.option(
    "--max-merge",
    help="Maximum number of chunks to merge before invoking recursive merging",
    type=int,
    default=200,
    show_default=True)
@click.option(
    "--storage-options",
    help="Options to modify the data filter pipeline. Provide as a "
         "comma-separated list of key-value pairs of the form 'k1=v1,k2=v2,...'. "
         "See http://docs.h5py.org/en/stable/high/dataset.html#filter-pipeline "
         "for more details.")
# @click.option(
#     "--format", "-f",
#     help="Preset data format.",
#     type=click.Choice(['4DN', 'BEDPE']))
# --sep
def pairs(bins, pairs_path, cool_path, metadata, assembly, chunksize,
          zero_based, comment_char, input_copy_status, no_symmetric_upper,
          field, temp_dir, no_delete_temp, max_merge, storage_options, **kwargs):
    """
    Bin any text file or stream of pairs.

    Pairs data need not be sorted. Accepts compressed files.
    To pipe input from stdin, set PAIRS_PATH to '-'.

    {}

    """
    chromsizes, bins = parse_bins(bins)

    symmetric_upper = not no_symmetric_upper
    tril_action = None
    if symmetric_upper:
        if input_copy_status == 'unique':
            tril_action = 'reflect'
        elif input_copy_status == 'duplex':
            tril_action = 'drop'

    if metadata is not None:
        with open(metadata, 'r') as f:
            metadata = json.load(f)

    input_field_names = [
        'chrom1', 'pos1', 'chrom2', 'pos2',
    ]
    input_field_dtypes = {
        'chrom1': str,
        'pos1': np.int64,
        'chrom2': str,
        'pos2': np.int64,
    }
    input_field_numbers = {}
    for name in ['chrom1', 'pos1', 'chrom2', 'pos2']:
        if kwargs[name] == 0:
            raise click.BadParameter(
                "Field numbers start at 1",
                param_hint=name)
        input_field_numbers[name] = kwargs[name] - 1

    # Include input value columns
    output_field_names = []
    output_field_dtypes = {}
    aggregations = {}
    if len(field):
        for arg in field:
            name, colnum, dtype, agg = parse_field_param(arg)

            # Special cases: these do not have input fields.
            # Omit field number and agg to change standard dtypes.
            if colnum is None:
                if (agg is None and dtype is not None
                        and name in {'bin1_id', 'bin2_id', 'count'}):
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

            if agg is not None:
                aggregations[name] = agg
            else:
                aggregations[name] = 'sum'

    # # Pairs counts are always produced, unless supressed explicitly
    # do_count = not no_count
    # if do_count:
    #     if 'count' not in output_field_names:
    #         output_field_names.append('count')  # default dtype and agg
    # else:
    #     if not len(output_field_names):
    #         click.BadParameter(
    #             "To pass `--no-count`, specify at least one input "
    #             "value-column using `--field`.")
    if 'count' not in output_field_names:
        output_field_names.append('count')

    # Customize the HDF5 filters
    if storage_options is not None:
        h5opts = parse_kv_list_param(storage_options)
        for key in h5opts:
            if isinstance(h5opts[key], list):
                h5opts[key] = tuple(h5opts[key])
    else:
        h5opts = None

    # Initialize the input stream
    if pairs_path == '-':
        f_in = sys.stdin
    else:
        f_in = pairs_path

    reader = pd.read_csv(
        f_in,
        sep='\t',
        usecols=[input_field_numbers[name] for name in input_field_names],
        names=input_field_names,
        dtype=input_field_dtypes,
        comment=comment_char,
        iterator=True,
        chunksize=chunksize)

    sanitize = sanitize_records(
        bins,
        schema='pairs',
        decode_chroms=True,
        is_one_based=not zero_based,
        tril_action=tril_action,
        sort=True,
        validate=True)
    aggregate = aggregate_records(agg=aggregations, count=True, sort=False)
    pipeline = compose(aggregate, sanitize)

    create_cooler(
        cool_path,
        bins,
        map(pipeline, reader),
        columns=output_field_names,
        dtypes=output_field_dtypes,
        metadata=metadata,
        assembly=assembly,
        mergebuf=chunksize,
        max_merge=max_merge,
        temp_dir=temp_dir,
        delete_temp=not no_delete_temp,
        boundscheck=False,
        triucheck=False,
        dupcheck=False,
        ensure_sorted=False,
        symmetric_upper=symmetric_upper,
        h5opts=h5opts,
        ordered=False
    )
