# -*- coding: utf-8 -*-
from __future__ import division, print_function
import json
import sys

import numpy as np
import pandas as pd
import h5py

import click
from . import cli
from ..io import create, TabixAggregator


@cli.command()
@click.argument(
    "bins_path",
    metavar="BINS_PATH")
@click.argument(
    "pairs_path",
    metavar="PAIRS_PATH")
@click.argument(
    "out",
    metavar="COOLER_PATH")
@click.option(
    "--metadata",
    help="Path to JSON file containing user metadata.")
@click.option(
    "--assembly",
    help="Name of genome assembly (e.g. hg19, mm10)")
def cload(bins_path, pairs_path, out, metadata, assembly):
    """
    Create a COOL file from a list of read pairs (i.e. contacts) and a list of bins.

    BINS_PATH : Path to BED file describing the genomic bin segmentation.
    PAIRS_PATH : Path to Tabix-indexed text file containing columns `chrom1`, `pos1`, `chrom2`, `pos2`.
    COOLER_PATH : Output COOL file path.

    """
    # Bin table
    bins = pd.read_csv(
        bins_path,
        sep='\t',
        names=['chrom', 'start', 'end'],
        dtype={'chrom': str})

    # Chrom table
    chromtable = (
        bins.drop_duplicates(['chrom'], keep='last')[['chrom', 'end']]
            .reset_index(drop=True)
            .rename(columns={'chrom': 'name', 'end': 'length'})
    )
    chroms, lengths = list(chromtable['name']), list(chromtable['length'])
    chromsizes = pd.Series(index=chroms, data=lengths)

    if metadata is not None:
        with open(metadata, 'r') as f:
            metadata = json.load(f)

    # Aggregate contacts
    chunksize = int(100e6)
    reader = TabixAggregator(pairs_path, chromsizes, bins)
    with h5py.File(out, 'w') as h5:
        create(h5, chroms, lengths, bins, reader, metadata, assembly)
