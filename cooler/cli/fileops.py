# -*- coding: utf-8 -*-
from __future__ import division, print_function
import click
from . import cli
from ..io import parse_cooler_uri
from .. import io


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
    List all Coolers inside a COOL file.

    """
    from ..api import Cooler
    for group_path in io.ls(cool_path):
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
# @click.option(
#     "--link", "-l",
#     help="If the source and destination file are the same, create a hard link "
#          "to the source group instead of a true copy.",
#     is_flag=True,
#     default=False)
# @click.option(
#     "--rename", "-m",
#     help="If the source and destination file are the same, create a hard link "
#          "to the source group and remove the original reference.",
#     is_flag=True,
#     default=False,)
# @click.option(
#     "--soft-link", "-s",
#     help="If the source and destination file are the same, create a soft link. "
#          "If the destination file is different, create an external link. This "
#          "type of link uses a path rather than a pointer.",
#     is_flag=True,
#     default=False)
def cp(src_uri, dst_uri, overwrite):
    """
    Copy a Cooler from one file to another or within the same file.

    See also: h5copy, h5repack tools from HDF5 suite

    \b\bArguments:

    SRC_URI : Path to source file or URI to source Cooler group

    DST_URI : Path to destination file or URI to destination Cooler group

    """
    io.cp(src_uri, dst_uri, overwrite)


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
    Create a hard link to a Cooler (rather than a true copy) in the same file.
    Also supports soft links (in the same file) or external links (different
    files).

    """
    io.ln(src_uri, dst_uri, overwrite, soft)


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
    Rename a Cooler within the same file.

    """
    io.mv(src_uri, dst_uri, overwrite)


@cli.command()
@click.argument(
    "uri",
    type=str)
@click.option(
    "-L", "--level",
    type=int)
def tree(uri, level):
    """
    Display the data hierarchy.

    """
    t = io.fileops.pprint_data_tree(uri, level)
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
    Display the attribute hierarchy.

    """
    t = io.fileops.pprint_attr_tree(uri, level)
    click.echo(t)


# @cli.command()
# def rename_chroms():
#     pass
