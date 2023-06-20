import shlex
import warnings
from math import ceil

import click

from .. import api
from ..parallel import lock
from ..reduce import (
    HIGLASS_TILE_DIM,
    legacy_zoomify,
    preferred_sequence,
    zoomify_cooler,
)
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
        if "weight" in api.Cooler(uri).bins():
            continue
        logger.info(f"Balancing zoom level with bin size {res}")

        try:
            balance_cmd.main(args=[uri, *args], prog_name="cooler")
        except SystemExit as e:
            # exc_info = sys.exc_info()
            exit_code = e.code

            if exit_code is None:
                exit_code = 0

            if exit_code != 0:
                raise e


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
    "e.g. --balance-args '--nproc 8 --ignore-diags 3'. Note that nproc for "
    "balancing must be specified independently of zoomify arguments.",
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
    "Repeat the `--field` option for each one. "
    "Use '<name>:dtype=<dtype>' to specify the dtype. Include "
    "',agg=<agg>' to specify an aggregation function different from 'sum'.",
    type=str,
    multiple=True,
)
@click.option(
    "--legacy",
    help="Use the legacy layout of integer-labeled zoom levels.",
    is_flag=True,
    default=False,
)
def zoomify(
    cool_uri,
    nproc,
    chunksize,
    resolutions,
    balance,
    balance_args,
    field,
    legacy,
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
        outfile = infile.replace(".cool", ".mcool")
    else:
        outfile, _ = parse_cooler_uri(out)

    logger.info(f'Recursively aggregating "{cool_uri}"')
    logger.info(f'Writing to "{outfile}"')

    if legacy:
        n_zooms, zoom_levels = legacy_zoomify(
            cool_uri, outfile, nproc, chunksize, lock=lock
        )

        if balance:
            from .balance import balance as balance_cmd

            if balance_args is None:
                balance_args = []
            else:
                balance_args = shlex.split(balance_args)
            logger.debug(f"Balancing args: {balance_args}")

            for level, res in reversed(list(zoom_levels.items())):
                uri = outfile + "::" + str(level)
                if level == str(n_zooms):
                    if "weight" in api.Cooler(uri).bins():
                        continue
                logger.info(f"Balancing zoom level {level}, bin size {res}")
                try:
                    balance_cmd.main(args=[uri, *balance_args], prog_name="cooler")
                except SystemExit as e:
                    # exc_info = sys.exc_info()
                    exit_code = e.code

                    if exit_code is None:
                        exit_code = 0

                    if exit_code != 0:
                        raise e

    else:
        clr = api.Cooler(cool_uri)
        genome_length = clr.chromsizes.values.sum()

        # Determine the coarsest resolution based on fitting the entire genome
        # in a single 256 x 256 tile.
        if clr.binsize:
            maxres = int(ceil(genome_length / HIGLASS_TILE_DIM))
            curres = clr.binsize
        else:
            bins = clr.bins()[["start", "end"]][:]
            mean_fragsize = (bins["end"] - bins["start"]).mean()
            maxres = int(ceil(genome_length / mean_fragsize / HIGLASS_TILE_DIM))
            curres = 1

        # Default is to use a binary geometric progression
        if resolutions is None:
            resolutions = "b"

        # Parse and expand user-provided resolutions
        resolutions, rstring = [], resolutions
        for res in [s.strip().lower() for s in rstring.split(",")]:
            if "n" in res or "b" in res and maxres < curres:
                warnings.warn(
                    "Map is already < 256 x 256. Provide resolutions "
                    "explicitly if you want to coarsen more."
                )
            if res == "n":
                r = preferred_sequence(curres, maxres, "nice")
            elif res == "b":
                r = preferred_sequence(curres, maxres, "binary")
            elif res == "4dn":
                r = [1000, 2000, *preferred_sequence(5000, maxres, "nice")]
            elif res.endswith("n"):
                res = int(res.split("n")[0])
                r = preferred_sequence(res, maxres, "nice")
            elif res.endswith("n"):
                res = int(res.split("b")[0])
                r = preferred_sequence(res, maxres, "binary")
            else:
                r = [int(res)]
            resolutions.extend(r)

        if len(field):
            field_specifiers = [
                parse_field_param(arg, includes_colnum=False) for arg in field
            ]
            columns, _, dtypes, agg = zip(*field_specifiers)
            columns = list(columns)
            dtypes = {col: dt for col, dt in zip(columns, dtypes) if dt is not None}
            agg = {col: f for col, f in zip(columns, agg) if f is not None}
        else:
            # If no other fields are given, 'count' is implicitly chosen.
            # Default aggregation. Dtype will be inferred.
            columns, dtypes, agg = ["count"], None, None

        # logger.info("Applying resolutions {}".format(resolutions))

        zoomify_cooler(
            [cool_uri, *list(base_uri)],
            outfile,
            resolutions,
            chunksize,
            nproc=nproc,
            lock=lock,
            columns=columns,
            dtypes=dtypes,
            agg=agg,
        )

        if balance:
            invoke_balance(balance_args, resolutions, outfile)
