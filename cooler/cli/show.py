# -*- coding: utf-8 -*-
from __future__ import division, print_function
import sys

import numpy as np
import h5py

from . import cli
import click

from ..api import Cooler
from .. import util


MAX_MATRIX_SIZE_FILE = int(1e8)
MAX_MATRIX_SIZE_INTERACTIVE = int(1e7)


def get_matrix_size(c, row_region, col_region):
    nrows = c.extent(row_region)[1] - c.extent(row_region)[0]
    ncols = c.extent(col_region)[1] - c.extent(col_region)[0]
    return ncols * nrows


def load_matrix(c, row_region, col_region, field, balanced, scale):
    mat = (c.matrix(balance=balanced, field=field)
            .fetch(row_region, col_region))

    if scale == 'log2':
        mat = np.log2(mat)
    elif scale == 'log10':
        mat = np.log10(mat)

    return mat


def interactive(ax, c, row_chrom, col_chrom, field, balanced, scale):
    import matplotlib.pyplot as plt
    # The code is heavily insired by
    # https://gist.github.com/mdboom/048aa35df685fe694330764894f0e40a

    def get_extent(ax):
        xstart, ystart, xdelta, ydelta = ax.viewLim.bounds
        xend = xstart + xdelta
        yend = ystart + ydelta
        return xstart, xend, ystart, yend

    def round_trim_extent(extent, binsize, row_chrom_len, col_chrom_len):
        xstart = int(np.floor(extent[0] / binsize) * binsize)
        xend = int(np.ceil(extent[1] / binsize) * binsize)
        ystart = int(np.floor(extent[3] / binsize) * binsize)
        yend = int(np.ceil(extent[2] / binsize) * binsize)
        xstart = max(0, xstart)
        ystart = max(0, ystart)
        # For now, don't let users to request the last bin, b/c its end
        # lies outside of the genome
        xend = min(xend, int(np.floor(col_chrom_len / binsize) * binsize))
        yend = min(yend, int(np.floor(row_chrom_len / binsize) * binsize))
        return xstart, xend, yend, ystart

    def update_heatmap(event):
        ax.set_autoscale_on(False)  # Otherwise, infinite loop

        extent = get_extent(ax)
        extent = round_trim_extent(extent, binsize, row_chrom_len, col_chrom_len)
        if extent == plotstate['prev_extent']:
            return

        plotstate['prev_extent'] = extent
        new_col_region = col_chrom, int(extent[0]), int(extent[1])
        new_row_region = row_chrom, int(extent[3]), int(extent[2])

        im = ax.images[-1]
        nelem = get_matrix_size(c, new_row_region, new_col_region)
        if nelem  >= MAX_MATRIX_SIZE_INTERACTIVE:
            # requested area too large
            im.set_data(np.ones(1)[:, None] * np.nan)

            if not plotstate['placeholders']:
                box, = plt.plot(
                    [0, col_chrom_len, col_chrom_len, 0, 0, col_chrom_len],
                    [0, row_chrom_len, 0, 0, row_chrom_len, row_chrom_len],
                    c='k',
                    lw=0.5
                )
                txt = plt.text(
                    0.5, 0.5,
                    'The requested region is too large\n'
                    'to display at this resolution.',
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=ax.transAxes
                )
                plotstate['placeholders'] = [box, txt]
        else:
            # remove placeholders if any and update
            while plotstate['placeholders']:
                plotstate['placeholders'].pop().remove()

            im.set_data(
                load_matrix(c, new_row_region, new_col_region, field, balanced, scale))

        im.set_extent(extent)
        ax.figure.canvas.draw_idle()

    binsize = c.info['bin-size']
    row_chrom_len = c.chromsizes[row_chrom]
    col_chrom_len = c.chromsizes[col_chrom]
    plotstate = {
        'placeholders': [],
        'prev_extent': get_extent(plt.gca())
    }
    plt.gcf().canvas.mpl_connect('button_release_event', update_heatmap)
    plt.show()


@cli.command()
@click.argument(
    "cool_uri",
    metavar="COOL_PATH")
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
@click.option(
    "--field",
    default='count',
    show_default=True,
    help="Pixel values to display.")
def show(cool_uri, range, range2, balanced, out, dpi, scale, force, zmin, zmax, cmap, field):
    """
    Display and browse a cooler in matplotlib.

    COOL_PATH : Path to a COOL file or Cooler URI.

    RANGE : The coordinates of the genomic region to display, in UCSC notation.
    Example: chr1:10,000,000-11,000,000

    """
    try:
        import matplotlib as mpl
        if out is not None:
            mpl.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("Install matplotlib to use cooler show", file=sys.stderr)
        sys.exit(1)

    c = Cooler(cool_uri)

    chromsizes = c.chromsizes
    row_region = range
    col_region = row_region if range2 is None else range2
    row_chrom, row_lo, row_hi = util.parse_region(row_region, chromsizes)
    col_chrom, col_lo, col_hi = util.parse_region(col_region, chromsizes)

    if ((get_matrix_size(c, row_region, col_region) >= MAX_MATRIX_SIZE_FILE)
        and not force):
        print(
            "The matrix of the selected region is too large. "
            "Try using lower resolution, selecting a smaller region, or use "
            "the '--force' flag to override this safety limit.",
            file=sys.stderr)
        sys.exit(1)

    plt.figure(figsize=(11, 10))
    plt.gcf().canvas.set_window_title('Contact matrix'.format())
    plt.title('')
    plt.imshow(
        load_matrix(c, row_region, col_region, field, balanced, scale),
        interpolation='none',
        extent=[col_lo, col_hi, row_hi, row_lo],
        vmin=zmin,
        vmax=zmax,
        cmap=cmap)

    # If plotting into a file, plot and quit
    plt.ylabel('{} coordinate'.format(row_chrom))
    plt.xlabel('{} coordinate'.format(col_chrom))
    cb = plt.colorbar()
    cb.set_label(
        {'linear': 'relative contact frequency',
         'log2'  : 'log 2 ( relative contact frequency )',
         'log10' : 'log 10 ( relative contact frequency )'}[scale])

    if out:
        plt.savefig(out, dpi=dpi)
    else:
        interactive(plt.gca(), c, row_chrom, col_chrom, field, balanced, scale)
