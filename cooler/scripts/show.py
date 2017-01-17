#!/usr/bin/env python
from __future__ import division, print_function
import argparse
import gzip
import sys

import numpy as np
import cooler, cooler.util
import matplotlib.pyplot as plt


def _str_to_num(string):
    if string is None:
        return string
    else:
        return float(string)


def main():
    parser = argparse.ArgumentParser(
        description="Display the contact matrix stored in a cooler file.")
    parser.add_argument(
        "cooler_file",
        help="Cooler file",
        metavar="COOLER_PATH")
    parser.add_argument(
        "--range", "-r",
        required=True,
        help="The coordinates of a displayed genomic region, in the UCSC notation."
        "Example: chr1:100,000,000-101,000,000")
    parser.add_argument(
        "--range2", "-r2",
        help="The coordinates of a genomic region shown along the columns. "
        "If omitted, the column range is the same as the row range. "
        "Use to display asymmetric matrices or trans interactions.")
    parser.add_argument(
        "--out", "-o",
        help="Save the image of the contact matrix to a file. "
        "If not specified, the matrix is displayed in an interactive window. "
        "The figure format is deduced from the extension of the file, "
        "the supported formats are png, jpg, svg, pdf, ps and eps.")
    parser.add_argument(
        "--dpi",
        help="The DPI of a figure, if saving to a file")
    parser.add_argument(
        "--raw",
        help="Show a raw (i.e. non-balanced) contact matrix. "
        "If not provided, display a balanced contact matrix.",
        action='store_true',
        default=False)
    parser.add_argument(
        "--linear", "-l",
        help="Do not perform a log10-transformation of the contact matrix"
        "If not provided, display log10 (contact matrix).",
        action='store_true',
        default=False)
    parser.add_argument(
        "--force", "-f",
        help="Force display very large matrices (>=10^8 pixels). "
        "Use at your own risk as it may cause performance issues.",
        action='store_true',
        default=False)
    parser.add_argument(
        "--zmin",
        help="The minimal displayed value of the matrix. If showing a log10 matrix, "
        "provide a log10 value as well."
        "To provide a negative value use a equal sign and quotes, e.g. -zmin='-0.5'"
        )
    parser.add_argument(
        "--zmax",
        help="The maximal displayed value of the matrix. If showing a log10 matrix, "
        "provide a log10 value as well."
        "To provide a negative value use a equal sign and quotes, e.g. -zmax='-0.5'"
        )
    parser.add_argument(
        "--cmap",
        help="The colormap used to display the contact matrix. "
        "See the full list at http://matplotlib.org/examples/color/colormaps_reference.html",
        default="YlOrRd"
    )
    args = vars(parser.parse_args())

    c = cooler.Cooler(args['cooler_file'])
    
    range_rows = args['range']
    range_cols = range_rows if args['range2'] is None else args['range2']
    nrows = c.extent(range_rows)[1] - c.extent(range_rows)[0]
    ncols = c.extent(range_cols)[1] - c.extent(range_cols)[0]
    if (ncols * nrows >= int(1e8)) and not args['force']:
        print(
            "The matrix of the selected region is too large. "
            "Try using lower resolution, selecting a smaller region, or use "
            "the '--force' flag to override this safety limit.")
        return

    mat = (c.matrix(balance=(not args['raw']))
            .fetch(range_rows, range_cols)
            .toarray())

    if not args['linear']:
        mat = np.log10(mat)

    chrm_row, lo_row, hi_row = cooler.util.parse_region_string(range_rows)
    chrm_col, lo_col, hi_col = cooler.util.parse_region_string(range_cols)
    vmin = _str_to_num(args['zmin'])
    vmax = _str_to_num(args['zmax'])

    plt.figure(figsize=(11,10))
    plt.gcf().canvas.set_window_title('Contact matrix'.format())
    plt.title('')
    plt.imshow(
        mat,
        interpolation='none',
        extent=[lo_col, hi_col, hi_row, lo_row],
        vmin=vmin,
        vmax=vmax,
        cmap=args['cmap'])
    plt.ylabel('{} coordinate'.format(chrm_row))
    plt.xlabel('{} coordinate'.format(chrm_col))
    plt.colorbar(label=(
        'relative contact frequency' if args['linear']
        else 'log 10 ( relative contact frequency )'))
    if args['out']:
        plt.savefig(args['out'], dpi=_str_to_num(args['dpi']))
    else:
        plt.show()
        

if __name__ == '__main__':
    main()
