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

import click
from . import cli, logger
from .. import util
from ..io import (
    create, create_from_unordered,
    sanitize_records, aggregate_records,
    TabixAggregator, HDF5Aggregator, PairixAggregator
)


@cli.group()
def cload():
    """
    Create a Cooler from a sorted list of contacts and a list of genomic bins.
    Choose a subcommand based on the format of the input contact list.

    """
    pass


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

    chromsizes, bins = _parse_bins(bins)

    if metadata is not None:
        with open(metadata, 'r') as f:
            metadata = json.load(f)

    with h5py.File(pairs_path, 'r') as h5pairs:
        iterator = HDF5Aggregator(h5pairs, chromsizes, bins, chunksize)
        create(cool_path, bins, iterator, metadata, assembly)


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
    "--max-split", "-s",
    help="Divide the pairs from each chromosome into at most this many chunks. "
         "Smaller chromosomes will be split less frequently or not at all. "
         "Increase ths value if large chromosomes dominate the workload on "
         "multiple processors.",
    type=int,
    default=2,
    show_default=True)
def tabix(bins, pairs_path, cool_path, metadata, assembly, nproc, max_split, **kwargs):
    """
    Bin a tabix-indexed contact list file.

    {}

    See also: 'cooler csort' to sort and index a contact list file

    Tabix manpage: <http://www.htslib.org/doc/tabix.html>.

    """
    chromsizes, bins = _parse_bins(bins)

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
        iterator = TabixAggregator(pairs_path, chromsizes, bins, map=map, n_chunks=max_split, **opts)
        create(cool_path, bins, iterator, metadata, assembly)
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
    "--max-split", "-s",
    help="Divide the pairs from each chromosome into at most this many chunks. "
         "Smaller chromosomes will be split less frequently or not at all. "
         "Increase ths value if large chromosomes dominate the workload on "
         "multiple processors.",
    type=int,
    default=2,
    show_default=True)
def pairix(bins, pairs_path, cool_path, metadata, assembly, nproc, max_split):
    """
    Bin a pairix-indexed contact list file.

    {}

    See also: 'cooler csort' to sort and index a contact list file

    Pairix on GitHub: <https://github.com/4dn-dcic/pairix>.

    """
    chromsizes, bins = _parse_bins(bins)

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
        iterator = PairixAggregator(pairs_path, chromsizes, bins, map=map, n_chunks=max_split)
        create(cool_path, bins, iterator, metadata, assembly)
    finally:
        if nproc > 1:
            pool.close() 


@register_subcommand
@add_arg_help
@click.option(
    "--chrom1", "-c1",
    help="chrom1 field number (one-based)",
    type=int,
    required=True)  #default=1)
@click.option(
    "--pos1", "-p1",
    help="pos1 field number (one-based)",
    type=int,
    required=True)  #default=2)
@click.option(
    "--chrom2", "-c2",
    help="chrom2 field number (one-based)",
    type=int,
    required=True)  #default=4)
@click.option(
    "--pos2", "-p2",
    help="pos2 field number (one-based)",
    type=int,
    required=True)  #default=5)
# @click.option(
#     "--format", "-f",
#     help="Preset data format.",
#     type=click.Choice(['4DN', 'BEDPE']))
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
    "--tril-action",
    type=click.Choice(['reflect', 'drop']),
    default='reflect',
    show_default=True,
    help="How to handle lower triangle records. " 
         "'reflect': make lower triangle records upper triangular. "
         "Use this if your input data comes only from a unique half of a "
         "symmetric matrix (but may not respect the specified chromosome order). "
         "'drop': discard all lower triangle records. Use this if your input "
         "data has mirror duplicates, i.e. is derived from a complete symmetric "
         "matrix.")
# @click.option(
#     "--field",
#     help="Add supplemental value fields or override default field numbers for "
#          "the specified format. Specify as '<name>,<number>' or as "
#          "'<name>,<number>,<dtype>' or '<name>,<number>,<dtype>,<agg>' to enforce a dtype other than `float` or "
#          "the default for a standard column. Field numbers are 1-based. "
#          "Repeat the `--field` option for each additional field. ",
#     type=str,
#     multiple=True)
# --sep
# --count-as-float
def pairs(bins, pairs_path, cool_path, metadata, assembly, chunksize, zero_based, comment_char, tril_action, **kwargs):
    """
    Bin any text file or stream of pairs.
    
    Pairs data need not be sorted. Accepts compressed files.
    To pipe input from stdin, set PAIRS_PATH to '-'.

    {}

    """
    chromsizes, bins = _parse_bins(bins)

    if metadata is not None:
        with open(metadata, 'r') as f:
            metadata = json.load(f)

    input_field_names = [
        'chrom1', 'pos1', 'chrom2', 'pos2' 
    ]
    input_field_dtypes = {
        'chrom1': str, 'pos1': int,
        'chrom2': str, 'pos2': int,
    }
    # input_field_numbers = {
    #     'chrom1': 0, 'pos1': 1, 
    #     'chrom2': 3, 'pos2': 4,
    # }
    input_field_numbers = {}
    for name in ['chrom1', 'pos1', 'chrom2', 'pos2']:
        if kwargs[name] == 0:
            raise click.BadParameter("Field numbers start at 1", 
                param_hint=name)
        input_field_numbers[name] = kwargs[name] - 1

    output_field_names = None
    output_field_dtypes = None

    if pairs_path == '-':
        f_in = sys.stdin
    else:
        f_in = pairs_path

    reader = pd.read_table(
        f_in, 
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
    aggregate = aggregate_records(agg=None, sort=False)
    pipeline = compose(aggregate, sanitize)

    create_from_unordered(
        cool_path, 
        bins, 
        map(pipeline, reader), 
        columns=output_field_names,
        dtypes=output_field_dtypes,
        metadata=metadata, 
        assembly=assembly,
        mergebuf=chunksize,
        boundscheck=False,
        triucheck=False,
        dupcheck=False,
        ensure_sorted=False
    )



