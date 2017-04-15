# -*- coding: utf-8 -*-
from __future__ import division, print_function
import click
from . import cli
from ..io import ls


@cli.command()
@click.argument(
    "cool_path",
    metavar="COOL_PATH")
def list(cool_path):
    """
    List all cooler "cans" inside a cooler file.

    """
    for can in ls(cool_path):
        click.echo(cool_path + '::' + can)
