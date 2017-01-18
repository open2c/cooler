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
from ..io import create, TabixAggregator, HDF5Aggregator, PairixAggregator


@cli.group()
def cload():
    """
    Create a COOL file from a sorted list of contacts and a list of genomic bins.
    Choose a subcommand based on the format of the input contact list.

    """
    pass


def _parse_bins(arg):
    # Enforce mutual exclusion
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
    elif op.exists(arg):
        # Bin table
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

        # Chrom table
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
            type=click.Path(exists=True),
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

    COOL_PATH : Output COOL file path.""")
    return func


@register_subcommand
@add_arg_help
def hiclib(bins, pairs_path, cool_path, metadata, assembly):
    """
    Bin a hiclib HDF5 contact list (frag) file.

    {}

    hiclib on BitBucket: <https://bitbucket.org/mirnylab/hiclib>.

    """

    chromsizes, bins = _parse_bins(bins)
    chroms, lengths = list(chromsizes.index), list(chromsizes.data)

    if metadata is not None:
        with open(metadata, 'r') as f:
            metadata = json.load(f)

    chunksize = int(100e6)
    with h5py.File(pairs_path, 'r') as h5pairs, \
         h5py.File(cool_path, 'w') as h5:
        iterator = HDF5Aggregator(h5pairs, chromsizes, bins, chunksize)
        create(h5, chroms, lengths, bins, iterator, metadata, assembly)


@register_subcommand
@add_arg_help
def tabix(bins, pairs_path, cool_path, metadata, assembly):
    """
    Bin a tabix-indexed contact list file.

    {}

    See also: 'cooler csort' to sort and index a contact list file

    Tabix manpage: <http://www.htslib.org/doc/tabix.html>.

    """
    chromsizes, bins = _parse_bins(bins)
    chroms, lengths = list(chromsizes.index), list(chromsizes.data)

    if metadata is not None:
        with open(metadata, 'r') as f:
            metadata = json.load(f)

    iterator = TabixAggregator(pairs_path, chromsizes, bins)

    with h5py.File(cool_path, 'w') as h5:
        create(h5, chroms, lengths, bins, iterator, metadata, assembly)   


@register_subcommand
@add_arg_help
def pairix(bins, pairs_path, cool_path, metadata, assembly):
    """
    Bin a pairix-indexed contact list file.

    {}

    See also: 'cooler csort' to sort and index a contact list file

    Pairix on GitHub: <https://github.com/4dn-dcic/pairix>.

    """
    chromsizes, bins = _parse_bins(bins)
    chroms, lengths = list(chromsizes.index), list(chromsizes.data)

    if metadata is not None:
        with open(metadata, 'r') as f:
            metadata = json.load(f)

    iterator = PairixAggregator(pairs_path, chromsizes, bins)

    with h5py.File(cool_path, 'w') as h5:
        create(h5, chroms, lengths, bins, iterator, metadata, assembly)   
