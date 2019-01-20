# -*- coding: utf-8 -*-
from __future__ import division, print_function
import numpy as np
import pandas as pd
from scipy import sparse
from cooler.util import parse_region_string
import cooler
import pytest


def test_region_string_parser():
    # UCSC-style names
    assert parse_region_string('chr21') == ('chr21', None, None)
    assert parse_region_string('chr21:1000-2000') == ('chr21', 1000, 2000)
    assert parse_region_string('chr21:1,000-2,000') == ('chr21', 1000, 2000)

    # Ensembl style names
    assert parse_region_string('6') == ('6', None, None)
    assert parse_region_string('6:1000-2000') == ('6', 1000, 2000)
    assert parse_region_string('6:1,000-2,000') == ('6', 1000, 2000)

    # FASTA style names
    assert parse_region_string('gb|accession|locus') == ('gb|accession|locus', None, None)
    assert parse_region_string('gb|accession|locus:1000-2000') == ('gb|accession|locus', 1000, 2000)
    assert parse_region_string('gb|accession|locus:1,000-2,000') == ('gb|accession|locus', 1000, 2000)

    # Punctuation in names (aside from :)
    assert parse_region_string('name-with-hyphens-') == ('name-with-hyphens-', None, None)
    assert parse_region_string('GL000207.1') == ('GL000207.1', None, None)
    assert parse_region_string('GL000207.1:1000-2000') == ('GL000207.1', 1000, 2000)

    # Trailing dash
    assert parse_region_string('chr21:1000-') == ('chr21', 1000, None)

    # Humanized units
    assert parse_region_string('6:1kb-2kb') == ('6', 1000, 2000)
    assert parse_region_string('6:1k-2000') == ('6', 1000, 2000)
    assert parse_region_string('6:1kb-2M') == ('6', 1000, 2000000)
    assert parse_region_string('6:1Gb-') == ('6', 1000000000, None)

    with pytest.raises(ValueError):
        parse_region_string('chr1:2,000-1,000')  # reverse selection
        parse_region_string('chr1::1000-2000')  # more than one colon


def test_selector1d():
    slicer = lambda fields, lo, hi: (lo, hi)
    fetcher = lambda x: x
    nmax = 50

    s = cooler.core.RangeSelector1D(None, slicer, fetcher, nmax)
    assert s[30] == (30, 31)
    assert s[10:20] == (10, 20)
    assert s[:20] == (0, 20)
    assert s[10:] == (10, nmax)
    assert s[:] == (0, nmax)
    assert s[:nmax] == (0, nmax)
    assert s[:-10] == (0, nmax-10)
    assert s[1:1] == (1, 1)
    with pytest.raises(IndexError):
        s[:, :]
    with pytest.raises(ValueError):
        s[::2]
    #assert_raises(TypeError, lambda : s['blah'])
    assert s.shape == (nmax,)

    # FIXME - questionable behavior
    assert s[30:20] == (30, 20)  # lo > hi
    assert s[nmax+10:nmax+30] == (nmax+10, nmax+30)  # lo > nmax
    assert s[10.0] == (10, 11)  # accepting floats
    #assert s[10.1] == (10.1, 11.1)  # not casting
    #assert s[nmax+10] == (nmax+10, nmax+11)


def test_selector2d():
    slicer = lambda field, i0, i1, j0, j1: (i0, i1, j0, j1)
    fetcher = lambda x: x
    nmax = 50

    s = cooler.core.RangeSelector2D(None, slicer, fetcher, (nmax, nmax))
    assert s[30] == (30, 31, 0, nmax)
    assert s[10:20, 10:20] == (10, 20, 10, 20)
    assert s[:] == (0, nmax, 0, nmax)
    with pytest.raises(IndexError):
        s[:, :, :]
    with pytest.raises(ValueError):
        s[::2, :]
    assert s.shape == (nmax, nmax)


def test_region_to_extent(mock_cooler):
    chromID_lookup = pd.Series({'chr1': 0, 'chr2': 1})
    binsize = 100

    region = ('chr1', 159, 402)
    first, last = 1, 4
    assert cooler.api.region_to_extent(
        mock_cooler, chromID_lookup, region, binsize) == (first, last+1)
    assert cooler.api.region_to_extent(
        mock_cooler, chromID_lookup, region, None) == (first, last+1)

    region = ('chr1', 159, 400)
    first, last = 1, 3
    assert cooler.api.region_to_extent(
        mock_cooler, chromID_lookup, region, binsize) == (first, last+1)
    assert cooler.api.region_to_extent(
        mock_cooler, chromID_lookup, region, None) == (first, last+1)


def test_slice_matrix(mock_cooler):
    slices = [
        (0, 10, 0, 10),
        (0, 10, 10, 20),
        (5, 15, 10, 20),
        (10, 20, 5, 15),
        (1, 1, 5, 15),
        (1, 1, 1, 1),
    ]
    for i0, i1, j0, j1 in slices:
        triu_reader = cooler.core.CSRReader(mock_cooler, 'count', max_chunk=10)

        # triangular query
        index = triu_reader.index_col(i0, i1, j0, j1)
        i, j, v = triu_reader.query(i0, i1, j0, j1)
        assert len(index) == len(i)

        # rectangular query
        i, j, v = cooler.core.query_rect(triu_reader.query, i0, i1, j0, j1)
        mat = sparse.coo_matrix((v, (i-i0, j-j0)), (i1-i0, j1-j0)).toarray()
        r = sparse.coo_matrix((
            (mock_cooler['pixels/count'],
            (mock_cooler['pixels/bin1_id'], mock_cooler['pixels/bin2_id']))
        ), (mock_cooler.attrs['nbins'],) * 2)
        r_full = r.toarray() + r.toarray().T
        assert np.allclose(r_full[i0:i1, j0:j1], mat)
