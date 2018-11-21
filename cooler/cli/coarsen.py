# -*- coding: utf-8 -*-
from __future__ import division, print_function
from collections import OrderedDict
from six.moves import map
import multiprocess as mp
import os.path as op
import shlex
import math
import sys

import numpy as np
import h5py

from ..io import parse_cooler_uri, create
from ..reduce import coarsen as _coarsen
from ..tools import lock
from . import cli
import click


@cli.command()
@click.argument(
    'cool_uri',
    metavar="COOL_PATH")
@click.option(
    '--factor', '-k',
    help="Gridding factor. The contact matrix is coarsegrained by grouping "
         "each chromosomal contact block into k-by-k element tiles",
    type=int,
    default=2,
    show_default=True)
@click.option(
    '--nproc', '-n', '-p',
    help="Number of processes to use for batch processing chunks of pixels "
         "[default: 1, i.e. no process pool]",
    default=1,
    type=int)
@click.option(
    '--chunksize', '-c',
    help="Number of pixels allocated to each process",
    type=int,
    default=int(10e6),
    show_default=True)
@click.option(
    '--out', '-o',
    required=True,
    help="Output file or URI")
def coarsen(cool_uri, factor, nproc, chunksize, out):
    """
    Coarsen a contact matrix by uniformly gridding the elements of each
    chromosomal block and summing the elements inside the grid tiles, i.e. a
    2-D histogram.

    \b\bArguments:

    COOL_PATH : Path to a COOL file or Cooler URI.

    """
    infile, _ = parse_cooler_uri(cool_uri)
    outfile, _ = parse_cooler_uri(out)
    same_file = op.realpath(infile) == op.realpath(outfile)
    _coarsen(cool_uri, out, factor, nproc, chunksize,
              lock=lock if same_file else None)
