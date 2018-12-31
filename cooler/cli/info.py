# -*- coding: utf-8 -*-
from __future__ import division, print_function
import json
import sys

from ._util import exit_on_broken_pipe
from . import cli
import click

from ..api import Cooler
from ..util import attrs_to_jsonable


@cli.command()
@click.argument(
    "cool_uri",
    type=str,
    metavar="COOL_PATH")
@click.option(
    "--field", "-f",
    help="Print the value of a specific info field.",
    type=str)
@click.option(
    "--metadata", "-m",
    help="Print the user metadata in JSON format.",
    is_flag=True,
    default=False)
@click.option(
    "--out", "-o",
    help="Output file (defaults to stdout)")
@exit_on_broken_pipe(1)
def info(cool_uri, field, metadata, out):
    """
    Display a cooler's info and metadata.

    COOL_PATH : Path to a COOL file or cooler URI.

    """
    c = Cooler(cool_uri)

    # Write output
    if out is None:
        f = sys.stdout
    else:
        f = open(out, 'wt')

    if metadata:
        json.dump(c.info['metadata'], f, indent=4)
        print(end='\n', file=f)

    elif field is not None:
        try:
            result = c.info[field]
        except KeyError:
            print("Data field {} not found.".format(field))
            sys.exit(1)
        print(result, file=f)

    else:
        dct = c.info.copy()
        dct.pop('metadata', None)
        json.dump(attrs_to_jsonable(dct), f, indent=4)
        print(end='\n', file=f)

    f.flush()
