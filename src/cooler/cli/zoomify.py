import os.path as op
import shlex
import warnings
from math import ceil

import click

from .. import api, fileops
from .._reduce import (
    _AGGS,
    HIGLASS_TILE_DIM,
    preferred_sequence,
    zoomify_cooler,
)
from ..parallel import get_lock
from ..util import parse_cooler_uri
from . import cli, get_logger
from ._util import parse_field_param


def invoke_balance(args, resolutions, outfile):
    from .balance import balance as balance_cmd

    logger = get_logger(__name__)
    if args is None:
        args = []
    else:
        args = shlex.split(args)
    logger.debug(f"Balancing args: {args}")

    for res in resolutions:
        uri = outfile + "::resolutions/" + str(res)
        ctx = balance_cmd.make_context("balance", [uri, *args], resilient_parsing=True)
        name = ctx.params.get("name") or "weight"
        if name in api.Cooler(uri).bins() and not ctx.params.get("force"):
            continue
        logger.info(f"Balancing zoom level with bin size {res}")

        try:
            balance_cmd.main(args=[uri, *args], prog_name="cooler")
        except SystemExit as e:
            exit_code = e.code

            if exit_code is None:
                exit_code = 0

            if exit_code != 0:
                raise


@cli.command()
@click.argument("cool_uri", metavar="COOL_PATH")
@click.option(
    "--nproc",
    "-n",
    "-p",
    help="Number of processes to use for batch processing chunks of pixels "
    "[default: 1, i.e. no process pool]",
    default=1,
    type=int,
)
@click.option(
    "--chunksize",
    "-c",
    help="Number of pixels allocated to each process",
    type=int,
    default=int(10e6),
    show_default=True,
)
@click.option(
    "--resolutions",
    "-r",
    help="Comma-separated list of target resolutions. Use suffixes B or N to "
    "specify a progression: B for binary (geometric steps of factor 2), N for "
    "nice (geometric steps of factor 10 interleaved with steps of 2 and 5). "
    "Examples: 1000B=1000,2000,4000,8000,... 1000N=1000,2000,5000,10000,... "
    "5000N=5000,10000,25000,50000,... 4DN is an alias for 1000,2000,5000N "
    "[default: B]",
)
@click.option(
    "--balance",
    help="Apply balancing to each zoom level. Off by default.",
    is_flag=True,
    default=False,
)
@click.option(
    "--balance-args",
    help="Additional arguments to pass to cooler balance. "
    "To deal with space ambiguity, use quotes to pass multiple arguments, "
    "e.g. ``--balance-args '--nproc 8 --ignore-diags 3'``. Note that nproc "
    "for balancing must be specified independently of zoomify arguments.",
    type=str,
)
@click.option(
    "--base-uri",
    "-i",
    help="One or more additional base coolers to aggregate from, if needed.",
    multiple=True,
)
@click.option("--out", "-o", help="Output file or URI")
@click.option(
    "--field",
    help="Specify the names of value columns to merge as '<name>'. "
    "Repeat the `--field` option for each one. Append ':' and "
    "comma-separated properties to set the dtype and/or the aggregation, "
    "e.g. 'count:agg=max' or 'count:dtype=float,agg=max'. The aggregation "
    "is sum (default), min or max.",
    type=str,
    multiple=True,
)
def zoomify(
    cool_uri,
    nproc,
    chunksize,
    resolutions,
    balance,
    balance_args,
    field,
    base_uri,
    out,
):
    """
    Generate a multi-resolution cooler file by coarsening.

    COOL_PATH : Path to a COOL file or Cooler URI.

    """
    logger = get_logger(__name__)
    infile, _ = parse_cooler_uri(cool_uri)

    if out is None:
        outfile = op.splitext(infile)[0] + ".mcool"
    else:
        outfile, group = parse_cooler_uri(out)
        if group != "/":
            raise click.BadParameter(
                "zoomify writes a multi-resolution file: give a file path, not "
                "a group.",
                param_hint="'--out'",
            )
    for uri in [cool_uri, *base_uri]:
        if op.exists(outfile) and op.samefile(parse_cooler_uri(uri)[0], outfile):
            raise click.BadParameter(
                f"The output file '{outfile}' contains the input '{uri}', which "
                "would be overwritten. Give another output file.",
                param_hint="'--out'",
            )

    logger.info(f'Recursively aggregating "{cool_uri}"')
    logger.info(f'Writing to "{outfile}"')

    clr = api.Cooler(cool_uri)
    genome_length = clr.chromsizes.values.sum()

    # Determine the coarsest resolution based on fitting the entire genome
    # in a single 256 x 256 tile.
    if clr.binsize:
        maxres = ceil(genome_length / HIGLASS_TILE_DIM)
        curres = clr.binsize
    else:
        bins = clr.bins()[["start", "end"]][:]
        mean_fragsize = (bins["end"] - bins["start"]).mean()
        maxres = ceil(genome_length / mean_fragsize / HIGLASS_TILE_DIM)
        curres = 1

    # Default is to use a binary geometric progression
    if resolutions is None:
        resolutions = "b"

    # Parse and expand user-provided resolutions
    resolutions, rstring = [], resolutions
    for res in [s.strip().lower() for s in rstring.split(",")]:
        if ("n" in res or "b" in res) and maxres < curres:
            warnings.warn(
                "Map is already < 256 x 256. Provide resolutions "
                "explicitly if you want to coarsen more.",
                stacklevel=1,
            )
        if res == "n":
            r = preferred_sequence(curres, maxres, "nice")
        elif res == "b":
            r = preferred_sequence(curres, maxres, "binary")
        elif res == "4dn":
            r = [1000, 2000, *preferred_sequence(5000, maxres, "nice")]
        else:
            try:
                if res.endswith("n"):
                    r = preferred_sequence(int(res[:-1]), maxres, "nice")
                elif res.endswith("b"):
                    r = preferred_sequence(int(res[:-1]), maxres, "binary")
                else:
                    r = [int(res)]
            except ValueError:
                raise click.BadParameter(
                    f"Invalid resolution '{res}'.", param_hint="'--resolutions'"
                ) from None
        resolutions.extend(r)

    if len(field):
        field_specifiers = [
            parse_field_param(arg, includes_colnum=False, aggs=_AGGS) for arg in field
        ]
        columns, _, dtypes, agg = zip(*field_specifiers)
        columns = list(columns)
        dtypes = {col: dt for col, dt in zip(columns, dtypes) if dt is not None}
        agg = {col: f for col, f in zip(columns, agg) if f is not None}
    else:
        # If no other fields are given, 'count' is implicitly chosen.
        # Default aggregation. Dtype will be inferred.
        columns, dtypes, agg = ["count"], None, None

    try:
        zoomify_cooler(
            [cool_uri, *list(base_uri)],
            outfile,
            resolutions,
            chunksize,
            nproc=nproc,
            lock=get_lock(),
            columns=columns,
            dtypes=dtypes,
            agg=agg,
        )
    except ValueError as e:
        # e.g. resolutions that can't be made from the base resolutions
        raise click.ClickException(str(e)) from e

    if balance:
        # Balance every level in the output, including the base(s)
        output_resolutions = sorted(
            int(path.split("/")[-1]) for path in fileops.list_coolers(outfile)
        )
        invoke_balance(balance_args, output_resolutions, outfile)
