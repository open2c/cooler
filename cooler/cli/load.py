# -*- coding: utf-8 -*-
from __future__ import division, print_function
import json
import sys

import numpy as np
import pandas as pd
import h5py

import click
from . import cli
from ..io import create, SparseLoader



@cli.command()
@click.argument(
    "bins_path",
    metavar="BINS_PATH")
@click.argument(
    "pixels_path",
    metavar="PIXELS_PATH")
@click.argument(
    "out",
    metavar="COOL_PATH")
@click.option(
    "--metadata",
    help="Path to JSON file containing user metadata.")
@click.option(
    "--assembly",
    help="Name of genome assembly (e.g. hg19, mm10)")
def load(bins_path, pixels_path, out, metadata, assembly):
    """
    Load a contact matrix.
    Load a sparse-formatted text dump of a contact matrix into a COOL file.

    BINS_PATH : BED-like file containing genomic bin segmentation

    PIXELS_PATH : Three-column, sorted sparse matrix text file in ijv-triple
    format, a.k.a. COO format. May be gzipped. Must be lexically sorted by 
    bin1_id and bin2_id.

    COOL_PATH : Output COOL file path

    """
    # Bin table
    bins = pd.read_csv(bins_path, sep='\t')

    # Chrom table
    chromtable = (
        bins.drop_duplicates(['chrom'], keep='last')[['chrom', 'end']]
            .reset_index(drop=True)
            .rename(columns={'chrom': 'name', 'end': 'length'})
    )
    chroms, lengths = list(chromtable['name']), list(chromtable['length'])

    # User-supplied JSON file
    if metadata is not None:
        with open(metadata, 'r') as f:
            metadata = json.load(f)

    # Load the binned contacts
    chunksize = int(100e6)
    reader = SparseLoader(pixels_path, chunksize)
    with h5py.File(out, 'w') as h5:
        create(h5, chroms, lengths, bins, reader, metadata, assembly)
