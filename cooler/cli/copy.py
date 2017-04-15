# -*- coding: utf-8 -*-
from __future__ import division, print_function
import click
from . import cli
from ..io import ls, parse_cooler_uri


@cli.command()
@click.argument(
    "src_uri",
    type=str)
@click.argument(
    "dst_uri",
    type=str)
@click.option(
    "--append", "a",
    help="Append data to file if it already exists instead of overwriting.",
    is_flag=True,
    default=False,
    show_default=True)
def copy(src_uri, dst_uri, append):
    """
    Copy a cooler can from one file to another.

    SRC_URI : Path to source file or URI to source Cooler group
    DST_URI : Path to destination file or URI to destination Cooler group

    """
    src_path, src_group = parse_cooler_uri(src_uri)
    dst_path, dst_group = parse_cooler_uri(dst_uri)

    write_mode = 'a' if append else 'w'

    with h5py.File(src_path, 'r') as src, \
         h5py.File(dst_path, write_mode) as dest:

        src.copy(src_group, dest, dst_group)
