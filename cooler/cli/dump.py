# -*- coding: utf-8 -*-
from __future__ import division, print_function
from six.moves import map
import gzip
import sys

import numpy as np
import pandas as pd

from ._util import DelimitedTuple, exit_on_broken_pipe
from . import cli
import click

from ..util import parse_region
from ..core import region_to_extent
from .. import api


def _comes_before(a0, a1, b0, b1, strict=False):
    if a0 < b0: return a1 <= b0 if strict else a1 <= b1
    return False


def _contains(a0, a1, b0, b1, strict=False):
    if a0 > b0 or a1 < b1: return False
    if strict and (a0 == b0 or a1 == b1): return False
    return a0 <= b0 and a1 >= b1


def _prune_partition(edges, step):
    edges = np.asarray(edges)
    cumlen = np.r_[0, np.cumsum(np.diff(edges))]
    cuts = [step * i  for i in range(0, int(np.ceil(cumlen[-1] / step)))]
    cuts.append(cumlen[-1])
    return np.unique(np.searchsorted(cumlen, cuts))


class CSRReader(object):
    def __init__(self, h5, field, chunksize):
        self.h5 = h5
        self.field = field
        self.chunksize = chunksize
        self.bin1_selector = h5['pixels']['bin1_id']
        self.bin2_selector = h5['pixels']['bin2_id']
        self.data_selector = h5['pixels'][field]
        self.offset_selector = h5['indexes']['bin1_offset']

    def __call__(self, i0, i1, j0, j1, transpose=False):
        isempty = True

        bin1_selector = self.bin1_selector
        bin2_selector = self.bin2_selector
        data_selector = self.data_selector
        chunksize = self.chunksize

        if (i1 - i0 > 0) or (j1 - j0 > 0):

            # coarsegrain the offsets to extract a big chunk of rows at a time
            offsets = self.offset_selector[i0:i1 + 1]
            which_offsets = _prune_partition(offsets, chunksize)

            for o0, o1 in zip(which_offsets[:-1], which_offsets[1:]):

                # extract a chunk of rows
                slc = slice(offsets[o0], offsets[o1])
                bin2_extracted = bin2_selector[slc]
                data_extracted = data_selector[slc]

                i, j, v = [], [], []

                # filter each row
                delta = offsets[o0]
                for o in range(o0, o1):
                    # correct the offsets
                    lo = offsets[o] - delta
                    hi = offsets[o + 1] - delta

                    # this row
                    bin2 = bin2_extracted[lo:hi]

                    # filter for the range of j values we want
                    mask = (bin2 >= j0) & (bin2 < j1)
                    cols = bin2[mask]

                    # apply same mask for data
                    data = data_extracted[lo:hi][mask]

                    # shortcut for row data
                    rows = np.full(len(cols), i0 + o, dtype=bin1_selector.dtype)

                    i.append(rows)
                    j.append(cols)
                    v.append(data)

                if len(i):
                    isempty = False
                    i = np.concatenate(i, axis=0)
                    j = np.concatenate(j, axis=0)
                    v = np.concatenate(v, axis=0)
                    if transpose:
                        i, j = j, i
                    yield i, j, v

        if isempty:
            i = np.array([], dtype=bin1_selector.dtype)
            j = np.array([], dtype=bin2_selector.dtype)
            v = np.array([], dtype=data_selector.dtype)
            if transpose:
                i, j = j, i
            yield i, j, v


def query2d(triu_reader, i0, i1, j0, j1, duplex):
    # symmetric query
    if (i0, i1) == (j0, j1):
        for i, j, v in triu_reader(i0, i1, i0, i1):
            if duplex:
                nodiag = i != j
                i, j, v = np.r_[i, j[nodiag]], np.r_[j, i[nodiag]], np.r_[v, v[nodiag]]
            yield i, j, v

    # asymmetric query
    else:
        transpose = False
        if j0 < i0 or (i0 == j0 and i1 < j1):
            i0, i1, j0, j1 = j0, j1, i0, i1
            transpose = True

        # non-overlapping
        if _comes_before(i0, i1, j0, j1, strict=True):
            for i, j, v in triu_reader(i0, i1, j0, j1, transpose):
                yield i, j, v

        # partially overlapping
        elif _comes_before(i0, i1, j0, j1):
            for i, j, v in triu_reader(i0, j0, j0, i1, transpose):
                yield i, j, v
            for i, j, v in triu_reader(j0, i1, j0, i1, transpose):
                if duplex:
                    nodiag = i != j
                    i, j, v = np.r_[i, j[nodiag]], np.r_[j, i[nodiag]], np.r_[v, v[nodiag]]
                yield i, j, v
            for i, j, v in triu_reader(i0, i1, i1, j1, transpose):
                yield i, j, v

        # nested
        elif _contains(i0, i1, j0, j1):
            for i, j, v in triu_reader(i0, j0, j0, j1, transpose):
                yield i, j, v
            for j, i, v in triu_reader(j0, j1, j0, j1, transpose):
                if duplex:
                    nodiag = i != j
                    i, j, v = np.r_[i, j[nodiag]], np.r_[j, i[nodiag]], np.r_[v, v[nodiag]]
                yield i, j, v
            for j, i, v in triu_reader(j0, j1, j1, i1, transpose):
                yield i, j, v

        else:
            raise IndexError("This shouldn't happen")


def make_annotator(bins, balanced, join, annotate, one_based_ids, one_based_starts):
    def annotator(chunk):
        if annotate is not None:
            extra_fields = list(annotate)
            try:
                extra_cols = bins[extra_fields]
            except KeyError as e:
                print('Column not found:\n {}'.format(e))
                sys.exit(1)
            extra = api.annotate(
                chunk[['bin1_id', 'bin2_id']], extra_cols,
                replace=True)

        if balanced:
            df = api.annotate(chunk, bins[['weight']])
            chunk['balanced'] = df['weight1'] * df['weight2'] * chunk['count']

        if join:
            chunk = api.annotate(chunk, bins[['chrom', 'start', 'end']],
                replace=True)

        if annotate is not None:
            chunk = pd.concat([chunk, extra], axis=1)

        if one_based_ids:
            for col in ['bin1_id', 'bin2_id']:
                if col in chunk.columns:
                    chunk[col] += 1

        if one_based_starts:
            for col in ['start1', 'start2']:
                if col in chunk.columns:
                    chunk[col] += 1

        return chunk
    return annotator


@cli.command()
@click.argument(
    "cool_uri",
    metavar="COOL_PATH")
@click.option(
    "--table", "-t",
    help="Which table to dump. Choosing 'chroms' or 'bins' will cause all "
         "pixel-related options to be ignored. Note that for coolers stored "
         "in symmetric-upper mode, 'pixels' only holds the upper triangle "
         "values of the matrix.",
    type=click.Choice(['chroms', 'bins', 'pixels']),
    default='pixels',
    show_default=True)
@click.option(
    "--columns", "-c",
    help="Restrict output to a subset of columns, provided as a "
         "comma-separated list.",
    type=DelimitedTuple(sep=','))
@click.option(
    "--header", "-H",
    help="Print the header of column names as the first row.",
    is_flag=True,
    default=False,
    show_default=True)
@click.option(
    "--na-rep",
    help="Missing data representation. Default is empty ''.",
    default='')
@click.option(
    "--float-format",
    help="Format string for floating point numbers (e.g. '.12g', '03.2f').",
    default='g',
    show_default=True)
@click.option(
    "--range", "-r",
    help="The coordinates of a genomic region shown along the row dimension, "
         "in UCSC-style notation. (Example: chr1:10,000,000-11,000,000). "
         "If omitted, the entire contact matrix is printed.",
    type=str)
@click.option(
    "--range2", "-r2",
    type=str,
    help="The coordinates of a genomic region shown along the column dimension. "
         "If omitted, the column range is the same as the row range.")
@click.option(
    "--matrix", "-m",
    help="For coolers stored in symmetric-upper mode, ensure any empty areas of "
         "the genomic query window are populated by generating the "
         "lower-triangular pixels.",
    is_flag=True,
    default=False,
    show_default=True)
@click.option(
    "--balanced/--no-balance", "-b",
    help="Apply balancing weights to data. This will print an extra column "
         "called `balanced`",
    is_flag=True,
    default=False,
    show_default=True)
@click.option(
    "--join",
    help="Print the full chromosome bin coordinates instead of bin IDs. "
         "This will replace the `bin1_id` column with `chrom1`, `start1`, and "
         "`end1`, and the `bin2_id` column with `chrom2`, `start2` and `end2`.",
    is_flag=True,
    default=False,
    show_default=True)
@click.option(
    "--annotate",
    help="Join additional columns from the bin table against the pixels. "
         "Provide a comma separated list of column names (no spaces). "
         "The merged columns will be suffixed by '1' and '2' accordingly.",
    type=DelimitedTuple(sep=','))
@click.option(
    "--one-based-ids",
    help="Print bin IDs as one-based rather than zero-based.",
    is_flag=True,
    default=False)
@click.option(
    "--one-based-starts",
    help="Print start coordinates as one-based rather than zero-based.",
    is_flag=True,
    default=False)
@click.option(
    '--chunksize', "-k",
    help="Sets the amount of pixel data loaded from disk at one time. "
         "Can affect the performance of joins on high resolution datasets. "
         "Default is to load as many rows as there are bins.",
    type=int)
@click.option(
    "--out", "-o",
    help="Output text file If .gz extension is detected, file is written "
         "using zlib. Default behavior is to stream to stdout.")
@exit_on_broken_pipe(1)
def dump(cool_uri, table, columns, header, na_rep, float_format,
         range, range2, matrix, balanced, join, annotate, one_based_ids,
         one_based_starts, chunksize, out):
    """
    Dump a cooler's data to a text stream.

    COOL_PATH : Path to COOL file or cooler URI.

    """
    c = api.Cooler(cool_uri)

    # output stream
    if out is None or out == '-':
        f = sys.stdout
    elif out.endswith('.gz'):
        f = gzip.open(out, 'wt')
    else:
        f = open(out, 'wt')

    # choose the source
    if table == 'chroms':
        selector = c.chroms()
        if columns is not None:
            selector = selector[list(columns)]
        chunks = (selector[:],)
    elif table == 'bins':
        selector = c.bins()
        if columns is not None:
            selector = selector[list(columns)]
        chunks = (selector[:],)
    else:
        # load all the bins
        bins = c.bins()[:]
        if chunksize is None:
            chunksize = len(bins)

        if balanced and 'weight' not in bins.columns:
            print('Balancing weights not found', file=sys.stderr)
            sys.exit(1)

        h5 = c.open('r')
        if range:
            i0, i1 = region_to_extent(
                h5,
                c._chromids,
                parse_region(range, c.chromsizes),
                binsize=c.binsize)
            if range2 is not None:
                j0, j1 = region_to_extent(
                    h5,
                    c._chromids,
                    parse_region(range2, c.chromsizes),
                    binsize=c.binsize)
            else:
                j0, j1 = i0, i1

            triu_reader = CSRReader(h5, 'count', chunksize)
            if matrix and c.storage_mode == 'symmetric-upper':
                selector = query2d(triu_reader, i0, i1, j0, j1, duplex=True)
            else:
                selector = triu_reader(i0, i1, j0, j1, transpose=False)

            chunks = (
                pd.DataFrame(
                    {'bin1_id': i , 'bin2_id': j, 'count': v},
                    columns=['bin1_id', 'bin2_id', 'count'])
                        for i, j, v in selector)
        else:
            selector = c.pixels()
            if columns is not None:
                selector = selector[list(columns)]
            n = len(selector)
            edges = np.arange(0, n + chunksize, chunksize)
            edges[-1] = n

            if matrix and c.storage_mode == 'symmetric-upper':
                def _select(lo, hi):
                    df = selector[lo:hi]
                    dfT = df.copy()
                    dfT['bin1_id'], dfT['bin2_id'] = df['bin2_id'], df['bin1_id']
                    return pd.concat([df, dfT])
                chunks = (_select(lo, hi) for lo, hi in zip(edges[:-1], edges[1:]))
            else:
                chunks = (selector[lo:hi] for lo, hi in zip(edges[:-1], edges[1:]))

        if balanced or join or annotate:
            annotator = make_annotator(bins, balanced, join, annotate,
                                       one_based_ids, one_based_starts)
            chunks = map(annotator, chunks)

    first = True
    if float_format is not None:
        float_format = '%' + float_format

    for chunk in chunks:
        if first:
            if header:
                chunk[0:0].to_csv(
                    f,
                    sep='\t',
                    index=False,
                    header=True,
                    float_format=float_format)
            first = False

        chunk.to_csv(
            f,
            sep='\t',
            index=False,
            header=False,
            float_format=float_format,
            na_rep=na_rep)

    else:
        f.flush()
