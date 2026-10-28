import json

import click

from ..api import Cooler
from ..util import attrs_to_jsonable
from . import cli
from ._util import exit_on_broken_pipe


@cli.command()
@click.argument("cool_uri", type=str, metavar="COOL_PATH")
@click.option(
    "--field", "-f", help="Print the value of a specific info field.", type=str
)
@click.option(
    "--metadata",
    "-m",
    help="Print the user metadata in JSON format.",
    is_flag=True,
    default=False,
)
@click.option("--out", "-o", help="Output file (defaults to stdout)")
@exit_on_broken_pipe(1)
def info(cool_uri, field, metadata, out):
    """
    Display a cooler's info and metadata.

    COOL_PATH : Path to a COOL file or cooler URI.

    """
    c = Cooler(cool_uri)

    if metadata:
        text = json.dumps(c.info["metadata"], indent=4)
    elif field is not None:
        if field not in c.info:
            raise click.BadParameter(
                f"No field '{field}' in the cooler's info.", param_hint="'--field'"
            )
        text = str(c.info[field])
    else:
        dct = c.info.copy()
        dct.pop("metadata", None)
        text = json.dumps(attrs_to_jsonable(dct), indent=4)

    if out is None:
        click.echo(text)
    else:
        with open(out, "w") as f:
            print(text, file=f)
