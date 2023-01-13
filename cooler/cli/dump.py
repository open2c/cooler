import gzip
import sys

import click
import pandas as pd

from .. import api
from ..core import (
    CSRReader,
    DirectRangeQuery2D,
    FillLowerRangeQuery2D,
    region_to_extent,
)
from ..util import parse_region
from . import cli
from ._util import DelimitedTuple


def make_annotator(bins, balanced, join, annotate, one_based_ids, one_based_starts):
    def annotator(chunk):
        if annotate is not None:
            extra_fields = list(annotate)
            try:
                extra_cols = bins[extra_fields]
            except KeyError as e:
                print(f"Column not found:\n {e}")
                sys.exit(1)
            extra = api.annotate(
                chunk[["bin1_id", "bin2_id"]], extra_cols, replace=True
            )

        if balanced:
            df = api.annotate(chunk, bins[["weight"]])
            chunk["balanced"] = df["weight1"] * df["weight2"] * chunk["count"]

        if join:
            chunk = api.annotate(chunk, bins[["chrom", "start", "end"]], replace=True)

        if annotate is not None:
            chunk = pd.concat([chunk, extra], axis=1)

        if one_based_ids:
            for col in ["bin1_id", "bin2_id"]:
                if col in chunk.columns:
                    chunk[col] += 1

        if one_based_starts:
            for col in ["start1", "start2"]:
                if col in chunk.columns:
                    chunk[col] += 1

        return chunk

    return annotator


@cli.command()
@click.argument(
    "cool_uri",
    metavar="COOL_PATH"
)
@click.option(
    "--table", "-t",
    help="Which table to dump. Choosing 'chroms' or 'bins' will cause all "
    "pixel-related options to be ignored. Note that for coolers stored "
    "in symmetric-upper mode, 'pixels' only holds the upper triangle "
    "values of the matrix.",
    type=click.Choice(["chroms", "bins", "pixels"]),
    default="pixels",
    show_default=True,
)
@click.option(
    "--columns", "-c",
    help="Restrict output to a subset of columns, provided as a "
    "comma-separated list.",
    type=DelimitedTuple(sep=","),
)
@click.option(
    "--header", "-H",
    help="Print the header of column names as the first row.",
    is_flag=True,
    default=False,
    show_default=True,
)
@click.option(
    "--na-rep",
    help="Missing data representation. Default is empty ''.",
    default=""
)
@click.option(
    "--float-format",
    help="Format string for floating point numbers (e.g. '.12g', '03.2f').",
    default="g",
    show_default=True,
)
@click.option(
    "--range", "-r",
    help="The coordinates of a genomic region shown along the row dimension, "
    "in UCSC-style notation. (Example: chr1:10,000,000-11,000,000). "
    "If omitted, the entire contact matrix is printed.",
    type=str,
)
@click.option(
    "--range2", "-r2",
    type=str,
    help="The coordinates of a genomic region shown along the column dimension. "
    "If omitted, the column range is the same as the row range.",
)
@click.option(
    "--fill-lower", "-f",
    help="For coolers using 'symmetric-upper' storage, populate implicit areas "
    "of the genomic query box by generating lower triangle pixels. If not "
    "specified, only upper triangle pixels are reported. This option has no "
    "effect on coolers stored in 'square' mode.",
    is_flag=True,
    default=False,
    show_default=True,
)
@click.option(
    "--balanced/--no-balance", "-b",
    help="Apply balancing weights to data. This will print an extra column "
    "called `balanced`",
    is_flag=True,
    default=False,
    show_default=True,
)
@click.option(
    "--join",
    help="Print the full chromosome bin coordinates instead of bin IDs. "
    "This will replace the `bin1_id` column with `chrom1`, `start1`, and "
    "`end1`, and the `bin2_id` column with `chrom2`, `start2` and `end2`.",
    is_flag=True,
    default=False,
    show_default=True,
)
@click.option(
    "--annotate",
    help="Join additional columns from the bin table against the pixels. "
    "Provide a comma separated list of column names (no spaces). "
    "The merged columns will be suffixed by '1' and '2' accordingly.",
    type=DelimitedTuple(sep=","),
)
@click.option(
    "--one-based-ids",
    help="Print bin IDs as one-based rather than zero-based.",
    is_flag=True,
    default=False,
)
@click.option(
    "--one-based-starts",
    help="Print start coordinates as one-based rather than zero-based.",
    is_flag=True,
    default=False,
)
@click.option(
    "--chunksize", "-k",
    help="Sets the number of pixel records loaded from disk at one time. "
    "Can affect the performance of joins on high resolution datasets. ",
    type=int,
    default=1_000_000,
    show_default=True,
)
@click.option(
    "--out", "-o",
    help="Output text file If .gz extension is detected, file is written "
    "using zlib. Default behavior is to stream to stdout.",
)
def dump(
    cool_uri,
    table,
    columns,
    header,
    na_rep,
    float_format,
    range,
    range2,
    fill_lower,
    balanced,
    join,
    annotate,
    one_based_ids,
    one_based_starts,
    chunksize,
    out,
):
    """
    Dump a cooler's data to a text stream.

    COOL_PATH : Path to COOL file or cooler URI.

    """
    clr = api.Cooler(cool_uri)

    # Choose the output stream
    if out is None or out == "-":
        f = sys.stdout
    elif out.endswith(".gz"):
        f = gzip.open(out, "wt")
    else:
        f = open(out, "w")

    # Choose the source table
    if table == "chroms":
        selector = clr.chroms()
        if columns is not None:
            selector = selector[list(columns)]
        chunks = (selector[:],)

    elif table == "bins":
        selector = clr.bins()
        if columns is not None:
            selector = selector[list(columns)]
        chunks = (selector[:],)

    else:  # Pixel table

        # Load all the bins
        bins = clr.bins()[:]
        n_bins = len(bins)
        if chunksize is None:
            chunksize = len(bins)

        if balanced and "weight" not in bins.columns:
            print("Balancing weights not found", file=sys.stderr)
            sys.exit(1)

        h5 = clr.open("r")
        reader = CSRReader(h5['pixels'], h5['indexes/bin1_offset'][:])
        field = "count"

        if range:
            # User-specified bbox provided
            i0, i1 = region_to_extent(
                h5,
                clr._chromids,
                parse_region(range, clr.chromsizes),
                binsize=clr.binsize
            )
            if range2 is not None:
                j0, j1 = region_to_extent(
                    h5,
                    clr._chromids,
                    parse_region(range2, clr.chromsizes),
                    binsize=clr.binsize,
                )
            else:
                j0, j1 = i0, i1
            bbox = (i0, i1, j0, j1)
        else:
            # Dump everything
            bbox = (0, n_bins, 0, n_bins)

        if fill_lower and clr.storage_mode == "symmetric-upper":
            engine = FillLowerRangeQuery2D(reader, field, bbox, chunksize)
        else:
            engine = DirectRangeQuery2D(reader, field, bbox, chunksize)

        chunks = (
            pd.DataFrame(
                dct, columns=["bin1_id", "bin2_id", field],
            ) for dct in engine
        )

        if balanced or join or annotate:
            annotator = make_annotator(
                bins, balanced, join, annotate, one_based_ids, one_based_starts
            )
            chunks = map(annotator, chunks)

    if float_format is not None:
        float_format = "%" + float_format

    is_first_chunk = True
    for chunk in chunks:
        if is_first_chunk:
            if header:
                chunk[0:0].to_csv(
                    f, sep="\t", index=False, header=True, float_format=float_format
                )
            is_first_chunk = False

        chunk.to_csv(
            f,
            sep="\t",
            index=False,
            header=False,
            float_format=float_format,
            na_rep=na_rep,
        )

    else:
        f.flush()
