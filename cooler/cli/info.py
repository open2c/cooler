# -*- coding: utf-8 -*-
from __future__ import division, print_function
import json
import sys

import click
from . import cli
from ..api import Cooler


@cli.command()
@click.argument(
    "cool_path",
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
def info(cool_path, field, metadata, out):
    """
    Display file info and metadata.

    COOL_PATH : Path to a COOL file.

    """
    c = Cooler(cool_path)

    # Write output
    try:
        if out is None:
            f = sys.stdout
        else:
            f = open(out, 'wt')

        if metadata:
            json.dump(c.info['metadata'], f, indent=4)
        elif field is not None:
            try:
                result = c.info[field]
            except KeyError:
                print("Data field {} not found.".format(field))
                sys.exit(1)
            print(result, file=f)
        else:
            dct = c.info
            for field in dct.keys():
                if field != 'metadata':
                    print(field + '\t' + str(dct[field]), file=f)

    except OSError:
        pass
    finally:
        f.close()
