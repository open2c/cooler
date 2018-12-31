# -*- coding: utf-8 -*-
from __future__ import division, print_function

from . import cli
import click

from ..util import parse_cooler_uri
from .. import fileops


@cli.command()
@click.argument(
    "cool_path",
    metavar="COOL_PATH")
@click.option(
    "--long", "-l",
    help="Long listing format",
    is_flag=True)
def ls(cool_path, long):
    """
    List all coolers inside a file.

    """
    from ..api import Cooler
    for group_path in fileops.list_coolers(cool_path):
        uri = cool_path + '::' + group_path
        if long:
            binsize = Cooler(uri).binsize
            if binsize is None:
                s = '{}\t<variable>'.format(uri)
            else:
                s = '{}\t{:,}'.format(uri, binsize)
            click.echo(s)
        else:
            click.echo(uri)


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
def cp(src_uri, dst_uri, overwrite):
    """
    Copy a cooler from one file to another or within the same file.

    See also: h5copy, h5repack tools from HDF5 suite.

    """
    # \b\bArguments:
    # SRC_URI : Path to source file or URI to source Cooler group
    # DST_URI : Path to destination file or URI to destination Cooler group
    # if dst_group in dst and dst_group != '/':
    #     click.confirm(
    #         "A group named '{}' already exists in '{}'. Overwrite?".format(
    #             dst_group, dst_path),
    #         abort=True)
    #     del dst[dst_group]
    fileops.cp(src_uri, dst_uri, overwrite)


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
    "--soft", "-s",
    help="Creates a soft link rather than a hard link if the source and "
         "destination file are the same. Otherwise, creates an external link. "
         "This type of link uses a path rather than a pointer.",
    is_flag=True,
    default=False)
def ln(src_uri, dst_uri, overwrite):
    """
    Create a hard link to a cooler (rather than a true copy) in the same file.
    Also supports soft links (in the same file) or external links (different
    files).

    """
    fileops.ln(src_uri, dst_uri, overwrite, soft)


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
def mv(src_uri, dst_uri, overwrite):
    """
    Rename a cooler within the same file.

    """
    fileops.mv(src_uri, dst_uri, overwrite)


@cli.command()
@click.argument(
    "uri",
    type=str)
@click.option(
    "-L", "--level",
    type=int)
def tree(uri, level):
    """
    Display a file's data hierarchy.

    """
    t = fileops.pprint_data_tree(uri, level)
    click.echo(t)


@cli.command()
@click.argument(
    "uri",
    type=str)
@click.option(
    "-L", "--level",
    type=int)
def attrs(uri, level):
    """
    Display a file's attribute hierarchy.

    """
    t = fileops.pprint_attr_tree(uri, level)
    click.echo(t)
