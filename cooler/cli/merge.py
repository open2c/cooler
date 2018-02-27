# -*- coding: utf-8 -*-
from __future__ import division, print_function
from multiprocess import Pool

import click
from . import cli, logger
from ..api import Cooler
from ..io import create, CoolerMerger


@cli.command()
@click.argument(
    "out_path",
    type=click.Path(exists=False))
@click.argument(
    "in_paths",
    nargs=-1,
    type=click.Path(exists=True))
@click.option(
    "--chunksize", "-c",
    type=int,
    default=int(20e6),
    show_default=True)
def merge(out_path, in_paths, chunksize):
    """
    Merge multiple contact matrices with identical axes.

    Data columns merged:

        pixels/bin1_id, pixels/bin2_id, pixels/count

    Data columns preserved:

        chroms/name, chroms/length
        bins/chrom, bins/start, bins/end

    Additional columns in the the input files are not preserved in the output.

    """
    logger.info("Merging:\n{}".format('\n'.join(in_paths)))
    clrs = [Cooler(path) for path in in_paths]
    chromsizes = clrs[0].chromsizes
    bins = clrs[0].bins()[['chrom', 'start', 'end']][:]
    assembly = clrs[0].info.get('genome-assembly', None)
    iterator = CoolerMerger(clrs, maxbuf=chunksize)
    create(out_path, bins, iterator, assembly=assembly)
