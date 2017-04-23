# -*- coding: utf-8 -*-
from __future__ import division, print_function
import click
from . import cli
from ..io import ls
from ..api import Cooler


@cli.command()
@click.argument(
    "cool_path",
    metavar="COOL_PATH")
@click.option(
    "--long", "-l",
    help="Long listing format",
    is_flag=True)
def list(cool_path, long):
    """
    List all Coolers inside a COOL file.

    """
    for group_path in ls(cool_path):
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
