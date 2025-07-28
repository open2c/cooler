import sys

import bioframe
import click
import h5py
import numpy as np
import pandas as pd
from multiprocess import Pool

from .._balance import balance_cooler
from ..api import Cooler
from ..util import bedslice, parse_cooler_uri
from . import cli, get_logger


@cli.command()
@click.argument("cool_uri", type=str, metavar="COOL_PATH")
@click.option(
    "--cis-only",
    help="Calculate weights against intra-chromosomal data only instead "
         "of genome-wide.",
    is_flag=True,
    default=False,
)
@click.option(
    "--trans-only",
    help="Calculate weights against inter-chromosomal data only instead "
         "of genome-wide.",
    is_flag=True,
    default=False,
)
@click.option(
    "--ignore-diags",
    help="Number of diagonals to ignore, including the main diagonal.",
    type=int,
    default=2,
    show_default=True,
)
@click.option(
    "--ignore-dist",
    help="Distance from the diagonal in bp to ignore.",
    type=int,
)
@click.option(
    "--mad-max",
    help="Ignore bins using the 'MAD-max' filter.",
    type=int,
    default=5,
    show_default=True,
)
@click.option(
    "--min-nnz",
    help="Ignore bins whose marginal number of nonzeros is less than this.",
    type=int,
    default=10,
    show_default=True,
)
@click.option(
    "--min-count",
    help="Ignore bins whose marginal count is less than this number.",
    type=int,
    default=0,
    show_default=True,
)
@click.option(
    "--blacklist",
    help="Path to a 3-column BED file of regions to mask out.",
    type=click.Path(exists=True),
)
@click.option(
    "--nproc",
    "-p",
    help="Number of processes to split the work between.",
    type=int,
    default=8,
    show_default=True,
)
@click.option(
    "--chunksize",
    "-c",
    help="Number of pixels handled by each worker process at a time.",
    type=int,
    default=int(10e6),
    show_default=True,
)
@click.option(
    "--tol",
    help="Threshold value of variance for convergence.",
    type=float,
    default=1e-5,
    show_default=True,
)
@click.option(
    "--max-iters",
    help="Maximum number of iterations if convergence not achieved.",
    type=int,
    default=200,
    show_default=True,
)
@click.option(
    "--name",
    help="Name of column to write to.",
    type=str,
    default="weight",
    show_default=True,
)
@click.option(
    "--force",
    "-f",
    help="Overwrite the target dataset if it exists.",
    is_flag=True,
    default=False,
)
@click.option(
    "--check",
    help="Check whether a 'weight' column already exists.",
    is_flag=True,
    default=False,
)
@click.option(
    "--stdout",
    help="Print weight column to stdout instead of saving to file.",
    is_flag=True,
    default=False,
)
@click.option(
    "--convergence-policy",
    help="What to do when balancing doesn't converge.",
    type=click.Choice(["store_final", "store_nan", "discard", "error"]),
    default="store_final",
    show_default=True,
)
def balance(
    cool_uri,
    nproc,
    chunksize,
    mad_max,
    min_nnz,
    min_count,
    blacklist,
    ignore_diags,
    tol,
    cis_only,
    trans_only,
    max_iters,
    name,
    force,
    check,
    stdout,
    convergence_policy,
    ignore_dist,
):
    """
    Out-of-core matrix balancing.

    Matrix must be symmetric. See the help for various filtering options to
    mask out poorly mapped bins.

    COOL_PATH : Path to a COOL file.
    """
    logger = get_logger(__name__)
    cool_path, group_path = parse_cooler_uri(cool_uri)

    if check:
        with h5py.File(cool_path, "r") as h5:
            grp = h5[group_path]
            if name not in grp["bins"]:
                click.echo(f"{cool_path}: No '{name}' column found.")
                sys.exit(1)
            else:
                click.echo(f"{cool_path}::{group_path} is balanced.")
                sys.exit(0)

    if cis_only and trans_only:
        raise click.UsageError(
            "Provide at most one of --cis-only and --trans-only flags"
        )

    with h5py.File(cool_path, "r+") as h5:
        grp = h5[group_path]
        if name in grp["bins"] and not stdout:
            if not force:
                print(
                    f"'{name}' column exists. Use --force to overwrite.",
                    file=sys.stderr,
                )
                sys.exit(1)
            else:
                del grp["bins"][name]

    logger.info(f'Balancing "{cool_uri}"')
    clr = Cooler(cool_uri)

    if blacklist is not None:
        with open(blacklist) as f:
            first_line = f.readline().strip()
            skiprows = 1 if first_line.startswith("track") else 0

        bad_regions = bioframe.read_table(
            blacklist,
            schema="bed3",
            skiprows=skiprows,
        )
        bad_regions = bad_regions.rename(
            columns={"chrom": "chrom", "start": "start", "end": "end"}
        )

        bins_grouped = clr.bins()[:].groupby("chrom", observed=True)
        chromsizes = clr.chromsizes

        bad_bins = []
        for _, reg in bad_regions.iterrows():
            result = bedslice(bins_grouped, chromsizes,
                              (reg.chrom, reg.start, reg.end))
            bad_bins.append(result.index.values)

        if bad_bins:
            bad_bins = np.concatenate(bad_bins)
        else:
            bad_bins = np.array([], dtype=int)
    else:
        bad_bins = None

    if ignore_dist is not None:
        ignore_diags = max(ignore_diags, int(np.ceil(ignore_dist / clr.binsize)))

    try:
        if nproc > 1:
            pool = Pool(nproc)
            map_ = pool.imap_unordered
        else:
            map_ = map

        bias, stats = balance_cooler(
            clr,
            chunksize=chunksize,
            cis_only=cis_only,
            trans_only=trans_only,
            tol=tol,
            min_nnz=min_nnz,
            min_count=min_count,
            blacklist=bad_bins,
            mad_max=mad_max,
            max_iters=max_iters,
            ignore_diags=ignore_diags,
            rescale_marginals=True,
            use_lock=False,
            map=map_,
        )

    finally:
        if nproc > 1:
            pool.close()

    if not np.all(stats["converged"]):
        logger.error("Iteration limit reached without convergence")
        if convergence_policy == "store_final":
            logger.error("Storing final result. Check log for convergence.")
        elif convergence_policy == "store_nan":
            logger.error("Saving weights as NaN.")
            bias[:] = np.nan
        elif convergence_policy == "discard":
            logger.error("Discarding result and aborting.")
            sys.exit(0)
        elif convergence_policy == "error":
            logger.error("Discarding result and aborting.")
            sys.exit(1)

    if stdout:
        pd.Series(bias).to_string(
            sys.stdout, header=False, index=False, na_rep="", float_format="%g"
        )
    else:
        with h5py.File(cool_path, "r+") as h5:
            grp = h5[group_path]
            h5opts = {"compression": "gzip", "compression_opts": 6,
                      "dtype": "float64"}
            grp["bins"].create_dataset(name, data=bias, **h5opts)
            grp["bins"][name].attrs.update(stats)
