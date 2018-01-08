# -*- coding: utf-8 -*-
from __future__ import division, print_function
import os
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
    default=False)
@click.option(
    "--link", "-l",
    help="If the source and destination file are the same, create a hard link "
         "to the source group instead of a true copy.",
    is_flag=True,
    default=False)
@click.option(
    "--rename", "-m",
    help="If the source and destination file are the same, create a hard link "
         "to the source group and remove the original reference.",
    is_flag=True,
    default=False,)
@click.option(
    "--soft-link", "-s",
    help="If the source and destination file are the same, create a soft link. "
         "If the destination file is different, create an external link. This "
         "type of link uses a path rather than a pointer.",
    is_flag=True,
    default=False)
def copy(src_uri, dst_uri, overwrite, link, rename, soft_link):
    """
    Copy a Cooler from one file to another or within the same file.

    See also: h5copy, h5repack tools from HDF5 suite

    \b\bArguments:

    SRC_URI : Path to source file or URI to source Cooler group

    DST_URI : Path to destination file or URI to destination Cooler group

    """
    src_path, src_group = parse_cooler_uri(src_uri)
    dst_path, dst_group = parse_cooler_uri(dst_uri)

    if sum([link, rename, soft_link]) > 1:
        raise click.BadParameter(
            'Must provide at most one of: --link, --rename, --soft-link')

    if not os.path.isfile(dst_path) or overwrite:
        write_mode = 'w'
    else:
        write_mode = 'r+'

    with h5py.File(src_path, 'r+') as src, \
         h5py.File(dst_path, write_mode) as dst:

        if dst_group in dst and dst_group != '/':
            click.confirm(
                "A group named '{}' already exists in '{}'. Overwrite?".format(
                    dst_group, dst_path), 
                abort=True)
            del dst[dst_group]

        if src_path == dst_path:
            if link or rename:
                src[dst_group] = src[src_group]
                if rename:
                    del src[src_group]
            elif soft_link:
                src[dst_group] = h5py.SoftLink(src_group)
        else:
            if link:
                click.echo("Can't hard link between two different files.")
                raise click.Abort
            elif soft_link:
                dst[dst_group] = h5py.ExternalLink(src_path, src_group)
            else:
                if dst_group == '/':
                    for subgrp in src[src_group].keys():
                        src.copy(src_group + '/' + subgrp, dst, subgrp)
                    dst[dst_group].attrs.update(src[src_group].attrs)
                else:
                    src.copy(
                        src_group, dst, 
                        dst_group if dst_group != '/' else None)
