# -*- coding: utf-8 -*-
from __future__ import division, print_function
import h5py

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
    "--overwrite", "-w",
    help="Truncate and replace destination file if it already exists.",
    is_flag=True,
    default=False,
    show_default=True)
@click.option(
    "--link", "-l",
    help="If the source and destination file are the same, create a hard link "
         "to the source group instead of a true copy.",
    is_flag=True,
    default=False,
    show_default=True)
@click.option(
    "--soft-link", "-s",
    help="If the source and destination file are the same, create a soft link. "
         "If the destination file is different, create an external link. These "
         "link types are essentially text paths rather than pointers.",
    is_flag=True,
    default=False,
    show_default=True)
def copy(src_uri, dst_uri, overwrite, link, soft_link):
    """
    Copy a cooler can from one file to another.

    SRC_URI : Path to source file or URI to source Cooler group
    DST_URI : Path to destination file or URI to destination Cooler group

    See also
    --------
    h5copy tool from HDF5 suite

    """
    src_path, src_group = parse_cooler_uri(src_uri)
    dst_path, dst_group = parse_cooler_uri(dst_uri)

    if overwrite:
        write_mode = 'w'
    else:
        write_mode = 'r+'

    with h5py.File(src_path, 'r+') as src, \
         h5py.File(dst_path, write_mode) as dst:

        if dst_group in dst:
            click.confirm(
                "A group named '{}' already exists in '{}'. Overwite?".format(
                    dst_group, dst_path), 
                abort=True)
            del dst[dst_group]

        if src_path == dst_path:
            if link:
                src[dst_group] = src[src_group]
            elif soft_link:
                src[dst_group] = h5py.SoftLink(src_group)
        else:
            if link:
                click.echo("Can't hard link between two different files.")
                raise click.Abort
            elif soft_link:
                dst[dst_group] = h5py.ExternalLink(src_path, src_group)
            else:
                src.copy(src_group, dst, 
                         dst_group if dst_group != '/' else None)
