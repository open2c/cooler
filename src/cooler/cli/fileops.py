import click

from .. import fileops
from . import cli


@cli.command()
@click.argument("cool_path", metavar="COOL_PATH")
@click.option("--long", "-l", help="Long listing format", is_flag=True)
def ls(cool_path, long):
    """
    List all coolers inside a file.

    COOL_PATH : Path to a file, or a URI to list only the coolers under a
    group, e.g. data.mcool::resolutions.

    """
    from ..api import Cooler
    from ..util import parse_cooler_uri

    file_path, root = parse_cooler_uri(cool_path)
    prefix = root.rstrip("/") + "/"
    for group_path in fileops.list_coolers(file_path):
        if root != "/" and group_path != root and not group_path.startswith(prefix):
            continue
        uri = file_path + "::" + group_path
        if long:
            binsize = Cooler(uri).binsize
            if binsize is None:
                s = f"{uri}\t<variable>"
            else:
                s = f"{uri}\t{binsize:,}"
            click.echo(s)
        else:
            click.echo(uri)


@cli.command()
@click.argument("src_uri", type=str)
@click.argument("dst_uri", type=str)
@click.option(
    "--overwrite",
    "-w",
    help="Replace the destination if it already exists. If the destination "
    "is the root of another file, the whole file is replaced.",
    is_flag=True,
    default=False,
)
def cp(src_uri, dst_uri, overwrite):
    """
    Copy a cooler from one file to another or within the same file.

    See also: h5copy, h5repack tools from HDF5 suite.

    """
    fileops.cp(src_uri, dst_uri, overwrite=overwrite)


@cli.command()
@click.argument("src_uri", type=str)
@click.argument("dst_uri", type=str)
@click.option(
    "--overwrite",
    "-w",
    help="Replace the destination if it already exists. If the destination "
    "is the root of another file, the whole file is replaced.",
    is_flag=True,
    default=False,
)
@click.option(
    "--soft",
    "-s",
    help="Creates a soft link rather than a hard link if the source and "
    "destination file are the same. Otherwise, creates an external link. "
    "This type of link uses a path rather than a pointer.",
    is_flag=True,
    default=False,
)
def ln(src_uri, dst_uri, overwrite, soft):
    """
    Create a hard link to a cooler (rather than a true copy) in the same file.
    Also supports soft links (in the same file) or external links (different
    files).

    """
    fileops.ln(src_uri, dst_uri, overwrite=overwrite, soft=soft)


@cli.command()
@click.argument("src_uri", type=str)
@click.argument("dst_uri", type=str)
@click.option(
    "--overwrite",
    "-w",
    help="Replace the destination if it already exists. If the destination "
    "is the root of another file, the whole file is replaced.",
    is_flag=True,
    default=False,
)
def mv(src_uri, dst_uri, overwrite):
    """
    Rename a cooler within the same file.

    """
    fileops.mv(src_uri, dst_uri, overwrite=overwrite)


@cli.command()
@click.argument("uri", type=str)
@click.option("-L", "--level", type=int)
def tree(uri, level):
    """
    Display a file's data hierarchy.

    """
    t = fileops.pprint_data_tree(uri, level)
    click.echo(t)


@cli.command()
@click.argument("uri", type=str)
@click.option("-L", "--level", type=int)
def attrs(uri, level):
    """
    Display a file's attribute hierarchy.

    """
    t = fileops.pprint_attr_tree(uri, level)
    click.echo(t)
