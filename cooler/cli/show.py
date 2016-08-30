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
    "--balanced", "-b",
    is_flag=True,
    default=False,
    help="Show the balanced contact matrix. "
         "If not provided, display the unbalanced counts.")
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
def show(cooler_file, range, range2, balanced, out, dpi, scale, force, zmin, zmax, cmap):
    """
    Display a contact matrix.
    Display a region of a contact matrix stored in a COOL file.

    COOLER_PATH : Path to a COOL file

    RANGE : The coordinates of the genomic region to display, in UCSC notation.
    Example: chr1:10,000,000-11,000,000

    """
    try:
        import matplotlib as mpl
        if out is not None:
            mpl.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.ticker
    except ImportError:
        print("Install matplotlib to use cooler show", file=sys.stderr)
        sys.exit(1)

    MAX_MATRIX_SIZE_FILE = int(1e8)
    MAX_MATRIX_SIZE_INTERACTIVE = int(1e7)
    c = Cooler(cooler_file)

    range_rows = range
    range_cols = range_rows if range2 is None else range2

    chrm_row, lo_row, hi_row = util.parse_region_string(range_rows)
    chrm_col, lo_col, hi_col = util.parse_region_string(range_cols)

    chrm_row_len = c.chroms()[:].set_index('name')['length'][chrm_row]
    chrm_col_len = c.chroms()[:].set_index('name')['length'][chrm_col]

    def get_matrix_size(range_rows, range_cols):
        nrows = c.extent(range_rows)[1] - c.extent(range_rows)[0]
        ncols = c.extent(range_cols)[1] - c.extent(range_cols)[0]
        return ncols * nrows

    if ((get_matrix_size(range_rows, range_cols) >= MAX_MATRIX_SIZE_FILE) 
        and not force):
        print(
            "The matrix of the selected region is too large. "
            "Try using lower resolution, selecting a smaller region, or use "
            "the '--force' flag to override this safety limit.",
            file=sys.stderr)
        sys.exit(1)

    def load_matrix(c, balanced, range_rows, range_cols, scale):
        mat = (c.matrix(balance=balanced)
                .fetch(range_rows, range_cols)
                .toarray())

        if scale == 'log2':
            mat = np.log2(mat)
        elif scale == 'log10':
            mat = np.log10(mat)

        return mat

    mat = load_matrix(c, balanced, range_rows, range_cols, scale)

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

    # If plotting into a file, plot and quit
    plt.ylabel('{} coordinate'.format(chrm_row))
    plt.xlabel('{} coordinate'.format(chrm_col))
    cb = plt.colorbar()
    cb.set_label(
        {'linear' : 'relative contact frequency',
         'log2' : 'log 2 ( relative contact frequency )',
         'log10' : 'log 10 ( relative contact frequency '}[scale])
    
    if out:
        plt.savefig(out, dpi=_str_to_num(dpi))
        sys.exit(0)

    # Otherwise, enter the interactive mode
    # The code is heavily insired by 
    # https://gist.github.com/mdboom/048aa35df685fe694330764894f0e40a
    
    def get_extent(ax):
        xstart, ystart, xdelta, ydelta = ax.viewLim.bounds
        xend = xstart + xdelta
        yend = ystart + ydelta
        return xstart, xend, ystart, yend

    def round_trim_extent(extent, chrm_row_len, chrom_col_len):
        binsize = c.info['bin-size']
        xstart = int(np.floor(extent[0] / binsize) * binsize)
        xend = int(np.ceil(extent[1] / binsize) * binsize)
        ystart = int(np.floor(extent[3] / binsize) * binsize)
        yend = int(np.ceil(extent[2] / binsize) * binsize)
        xstart = max(0, xstart)
        ystart = max(0, ystart)
        # For now, don't let users to request the last bin, b/c its end
        # lies outside of the genome
        xend = min(xend, int(np.floor(chrm_col_len / binsize) * binsize))
        yend = min(yend, int(np.floor(chrm_row_len / binsize) * binsize))

        return (xstart, xend, yend, ystart)

    plotstate = {
        'ax' : plt.gca(),
        'prev_extent' : get_extent(plt.gca())}

    def update_heatmap(event):
        ax = plotstate['ax']
        ax.set_autoscale_on(False)  # Otherwise, infinite loop

        extent = get_extent(ax)
        extent = round_trim_extent(extent, chrm_row_len, chrm_col_len)
        if extent == plotstate['prev_extent']:
            return
        plotstate['prev_extent'] = extent

        #print('move to ', extent)

        new_range_cols = '{}:{}-{}'.format(
            chrm_col, int(extent[0]), int(extent[1]))
        new_range_rows = '{}:{}-{}'.format(
            chrm_row, int(extent[3]), int(extent[2]))

        im = ax.images[-1]

        if get_matrix_size(new_range_rows, new_range_cols) >= MAX_MATRIX_SIZE_INTERACTIVE:
            im.set_data(np.ones(1)[:,None] * np.nan)
            if 'placeholders' not in plotstate:
                box, = plt.plot(
                    [0, chrm_col_len, chrm_col_len, 0, 0, chrm_col_len],
                    [0, chrm_row_len, 0, 0, chrm_row_len, chrm_row_len],
                    c='k',
                    lw=0.5)
                txt = plt.text(0.5, 0.5,
                    'The requested region is too large\n'
                    'to display at this resolution.',
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform = ax.transAxes)
                plotstate['placeholders'] = [box, txt]


        else:
            if 'placeholders' in plotstate:
                while plotstate['placeholders']:
                    plotstate['placeholders'].pop().remove()
                del(plotstate['placeholders'])

            mat = load_matrix(c, balanced, new_range_rows, new_range_cols, scale)
            im.set_data(mat)

        im.set_extent(extent)
        ax.figure.canvas.draw_idle()

    plt.gcf().canvas.mpl_connect('button_release_event', update_heatmap)
    plt.show()
