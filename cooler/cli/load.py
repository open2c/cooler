# -*- coding: utf-8 -*-
from __future__ import division, print_function
import json
import sys

import numpy as np
import pandas as pd
import h5py

import click
from . import cli
from ..io import create, SparseLoader #, BedGraph2DLoader

# TODO:
# cooler load bedpe
# cooler load coo
# cooler load block_coo


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


def iter_sparse(filepath, bins, chunksize):
    """
    Contact iterator for a sparse tsv Hi-C matrix with fields:
        "chrom1, start1, end1, chrom2, start2, end2, count"
    
    The fields are assumed to be defined and records assumed to 
    be sorted consistently with the bin table provided.
    
    Parameters
    ----------
    filepath : str
        Path to tsv file
    bins : DataFrame
        A bin table dataframe
    chunksize : number of rows of the matrix file to read at a time
    
    """
    iterator = pandas.read_csv(
        filepath, 
        sep='\t', 
        iterator=True,
        chunksize=chunksize,
        names=['chrom1', 'start1', 'end1', 
               'chrom2', 'start2', 'end2', 'count'])
    bins['bin'] = bins.index
    
    for chunk in iterator:
        # assign bin IDs from bin table
        df = (chunk.merge(bins, 
                          left_on=['chrom1', 'start1', 'end1'], 
                          right_on=['chrom', 'start', 'end'])
                   .merge(bins, 
                          left_on=['chrom2', 'start2', 'end2'], 
                          right_on=['chrom', 'start', 'end'], 
                          suffixes=('1', '2')))
        df = (df[['bin1', 'bin2', 'count']]
                .rename(columns={'bin1': 'bin1_id', 
                                 'bin2': 'bin2_id'})
                .sort_values(['bin1_id', 'bin2_id']))
        yield {k: v.values for k,v in six.iteritems(df)}



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
    help="",
    type=click.Choice(['BEDPE', 'COO']),
    default='BEDPE',
    show_default=True)
@click.option(
    "--metadata",
    help="Path to JSON file containing user metadata.")
@click.option(
    "--assembly",
    help="Name of genome assembly (e.g. hg19, mm10)")
def load(bins_path, pixels_path, cool_path, metadata, assembly):
    """
    Load a contact matrix.
    Load a sparse-formatted text dump of a contact matrix into a COOL file.

    Three-column, sorted sparse matrix text file in ijv-triple
    format, a.k.a. COO format.

    BINS_PATH : BED-like file containing genomic bin segmentation

    PIXELS_PATH : Text file containing nonzero pixel values. May be gzipped. 
    Must be lexically sorted by bin1_id and bin2_id.

    COOL_PATH : Output COOL file path

    """
    chromsizes, bins = _parse_bins(bins)

    # User-supplied JSON file
    if metadata is not None:
        with open(metadata, 'r') as f:
            metadata = json.load(f)

    # Load the binned contacts
    chunksize = int(100e6)
    iterator = SparseLoader(pixels_path, chunksize)
    create(cool_path, chromsizes, bins, iterator, metadata, assembly)
