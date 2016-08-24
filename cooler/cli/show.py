# -*- coding: utf-8 -*-
from __future__ import division, print_function
import sys

import numpy as np
import h5py

import click
from . import cli
from ..api import Cooler
from .. import util

def _str_to_num(string):
    if string is None:
        return string
    else:
        return float(string)


@cli.command()
@click.argument(
    "cooler_file",
    metavar="COOLER_PATH")
@click.argument(
    "range",
    type=str)
@click.option(
    "--range2", "-r2",
    type=str,
    help="The coordinates of a genomic region shown along the column dimension. "
         "If omitted, the column range is the same as the row range. "
         "Use to display asymmetric matrices or trans interactions.")
@click.option(
    "--out", "-o",
    help="Save the image of the contact matrix to a file. "
         "If not specified, the matrix is displayed in an interactive window. "
         "The figure format is deduced from the extension of the file, "
         "the supported formats are png, jpg, svg, pdf, ps and eps.")
@click.option(
    "--dpi",
    type=int,
    help="The DPI of the figure, if saving to a file")
@click.option(
    "--raw",
    is_flag=True,
    default=False,
    help="Show the raw (i.e. unbalanced) contact matrix. "
         "If not provided, display the balanced contact matrix.")
@click.option('--scale', '-s',
    type=click.Choice(['linear', 'log2', 'log10']),
    help="Scale transformation of the colormap: linear, log2 or log10. "
         "Default is log10.",
    default='log10')
@click.option(
    "--force", "-f",
    is_flag=True,
    default=False,
    help="Force display very large matrices (>=10^8 pixels). "
         "Use at your own risk as it may cause performance issues.")
@click.option(
    "--zmin",
    type=float,
    help="The minimal value of the color scale. Units must match those of the colormap scale. "
         "To provide a negative value use a equal sign and quotes, e.g. -zmin='-0.5'")
@click.option(
    "--zmax",
    type=float,
    help="The maximal value of the color scale. Units must match those of the colormap scale. "
         "To provide a negative value use a equal sign and quotes, e.g. -zmax='-0.5'")
@click.option(
    "--cmap",
    default="YlOrRd",
    help="The colormap used to display the contact matrix. "
         "See the full list at http://matplotlib.org/examples/color/colormaps_reference.html")
def show(cooler_file, range, range2, out, dpi, raw, scale, force, zmin, zmax, cmap):
    """
    Display a contact matrix.
    Display a region of a contact matrix stored in a COOL file.

    COOLER_PATH : Path to a COOL file

    RANGE : The coordinates of the genomic region to display, in UCSC notation.
    Example: chr1:10,000,000-11,000,000

    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Install matplotlib to use cooler show", file=sys.stderr)
        sys.exit(1)

    c = Cooler(cooler_file)

    range_rows = range
    range_cols = range_rows if range2 is None else range2
    nrows = c.extent(range_rows)[1] - c.extent(range_rows)[0]
    ncols = c.extent(range_cols)[1] - c.extent(range_cols)[0]
    if (ncols * nrows >= int(1e8)) and not force:
        print(
            "The matrix of the selected region is too large. "
            "Try using lower resolution, selecting a smaller region, or use "
            "the '--force' flag to override this safety limit.",
            file=sys.stderr)
        sys.exit(1)

    mat = (c.matrix(balance=(not raw))
            .fetch(range_rows, range_cols)
            .toarray())

    if scale == 'log2':
        mat = np.log2(mat)
    elif scale == 'log10':
        mat = np.log10(mat)

    chrm_row, lo_row, hi_row = util.parse_region_string(range_rows)
    chrm_col, lo_col, hi_col = util.parse_region_string(range_cols)
    vmin = _str_to_num(zmin)
    vmax = _str_to_num(zmax)

    plt.figure(figsize=(11,10))
    plt.gcf().canvas.set_window_title('Contact matrix'.format())
    plt.title('')
    plt.imshow(
        mat,
        interpolation='none',
        extent=[lo_col, hi_col, hi_row, lo_row],
        vmin=vmin,
        vmax=vmax,
        cmap=cmap)
    plt.ylabel('{} coordinate'.format(chrm_row))
    plt.xlabel('{} coordinate'.format(chrm_col))
    plt.colorbar(label=(
        {'linear' : 'relative contact frequency',
         'log2' : 'log 2 ( relative contact frequency )',
         'log10' : 'log 10 ( relative contact frequency '}[scale]))
    if out:
        plt.savefig(out, dpi=_str_to_num(dpi))
    else:
        plt.show()
