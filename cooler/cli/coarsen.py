import os.path as op

import click

from ..parallel import lock
from ..reduce import coarsen_cooler
from ..util import parse_cooler_uri
from . import cli
from ._util import parse_field_param


@cli.command()
@click.argument(
    "cool_uri",
    metavar="COOL_PATH"
)
@click.option(
    "--factor", "-k",
    help="Gridding factor. The contact matrix is coarsegrained by grouping "
    "each chromosomal contact block into k-by-k element tiles",
    type=int,
    default=2,
    show_default=True,
)
@click.option(
    "--nproc", "-n", "-p",
    help="Number of processes to use for batch processing chunks of pixels "
    "[default: 1, i.e. no process pool]",
    default=1,
    type=int,
)
@click.option(
    "--chunksize", "-c",
    help="Number of pixels allocated to each process",
    type=int,
    default=int(10e6),
    show_default=True,
)
@click.option(
    "--field",
    help="Specify the names of value columns to merge as '<name>'. "
    "Repeat the `--field` option for each one. "
    "Use '<name>,dtype=<dtype>' to specify the dtype. Include "
    "',agg=<agg>' to specify an aggregation function different from 'sum'.",
    type=str,
    multiple=True,
)
@click.option(
    "--out", "-o",
    required=True,
    help="Output file or URI"
)
@click.option(
    "--append", "-a",
    is_flag=True,
    default=False,
    help="Pass this flag to append the output cooler to an existing file "
         "instead of overwriting the file."
)
def coarsen(cool_uri, factor, nproc, chunksize, field, out, append):
    """
    Coarsen a cooler to a lower resolution.

    Works by pooling *k*-by-*k* neighborhoods of pixels and aggregating.
    Each chromosomal block is coarsened individually.

    COOL_PATH : Path to a COOL file or Cooler URI.

    """
    infile, _ = parse_cooler_uri(cool_uri)
    outfile, _ = parse_cooler_uri(out)
    same_file = op.realpath(infile) == op.realpath(outfile)

    if len(field):
        field_specifiers = [
            parse_field_param(arg, includes_colnum=False) for arg in field
        ]
        columns, _, dtypes, agg = zip(*field_specifiers)
        dtypes = {col: dt for col, dt in zip(columns, dtypes) if dt is not None}
        agg = {col: f for col, f in zip(columns, agg) if f is not None}
    else:
        # If no other fields are given, 'count' is implicitly chosen.
        # Default aggregation. Dtype will be inferred.
        columns, dtypes, agg = ["count"], None, None

    coarsen_cooler(
        cool_uri,
        out,
        factor,
        chunksize=chunksize,
        nproc=nproc,
        columns=columns,
        dtypes=dtypes,
        agg=agg,
        lock=lock if same_file else None,
        mode="a" if append else "w"
    )
