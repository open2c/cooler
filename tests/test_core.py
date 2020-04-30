from __future__ import division, print_function
from scipy import sparse
import numpy as np
import pandas as pd
import pytest

from cooler import core


def test_selector1d():
    slicer = lambda fields, lo, hi: (lo, hi)  # noqa
    fetcher = lambda x: x  # noqa
    nmax = 50

    s = core.RangeSelector1D(None, slicer, fetcher, nmax)
    assert s[30] == (30, 31)
    assert s[10:20] == (10, 20)
    assert s[:20] == (0, 20)
    assert s[10:] == (10, nmax)
    assert s[:] == (0, nmax)
    assert s[:nmax] == (0, nmax)
    assert s[:-10] == (0, nmax - 10)
    assert s[1:1] == (1, 1)
    with pytest.raises(IndexError):
        s[:, :]
    with pytest.raises(ValueError):
        s[::2]
    # assert_raises(TypeError, lambda : s['blah'])
    assert s.shape == (nmax,)

    # FIXME - questionable behavior
    assert s[30:20] == (30, 20)  # lo > hi
    assert s[nmax + 10 : nmax + 30] == (nmax + 10, nmax + 30)  # lo > nmax
    assert s[10.0] == (10, 11)  # accepting floats
    # assert s[10.1] == (10.1, 11.1)  # not casting
    # assert s[nmax+10] == (nmax+10, nmax+11)


def test_selector2d():
    slicer = lambda field, i0, i1, j0, j1: (i0, i1, j0, j1)  # noqa
    fetcher = lambda x: x  # noqa
    nmax = 50

    s = core.RangeSelector2D(None, slicer, fetcher, (nmax, nmax))
    assert s[30] == (30, 31, 0, nmax)
    assert s[10:20, 10:20] == (10, 20, 10, 20)
    assert s[:] == (0, nmax, 0, nmax)
    with pytest.raises(IndexError):
        s[:, :, :]
    with pytest.raises(ValueError):
        s[::2, :]
    assert s.shape == (nmax, nmax)


def test_region_to_extent(mock_cooler):
    chromID_lookup = pd.Series({"chr1": 0, "chr2": 1})
    binsize = 100

    region = ("chr1", 159, 402)
    first, last = 1, 4
    assert core.region_to_extent(
        mock_cooler, chromID_lookup, region, binsize
    ) == (first, last + 1)
    assert core.region_to_extent(mock_cooler, chromID_lookup, region, None) == (
        first,
        last + 1,
    )

    region = ("chr1", 159, 400)
    first, last = 1, 3
    assert core.region_to_extent(
        mock_cooler, chromID_lookup, region, binsize
    ) == (first, last + 1)
    assert core.region_to_extent(mock_cooler, chromID_lookup, region, None) == (
        first,
        last + 1,
    )


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
        triu_reader = core.CSRReader(mock_cooler, "count", max_chunk=10)

        # triangular query
        index = triu_reader.index_col(i0, i1, j0, j1)
        i, j, v = triu_reader.query(i0, i1, j0, j1)
        assert len(index) == len(i)

        # rectangular query
        i, j, v = core.query_rect(triu_reader.query, i0, i1, j0, j1)
        mat = sparse.coo_matrix((v, (i - i0, j - j0)), (i1 - i0, j1 - j0)).toarray()
        r = sparse.coo_matrix(
            (
                (
                    mock_cooler["pixels/count"],
                    (mock_cooler["pixels/bin1_id"], mock_cooler["pixels/bin2_id"]),
                )
            ),
            (mock_cooler.attrs["nbins"],) * 2,
        )
        r_full = r.toarray() + r.toarray().T
        assert np.allclose(r_full[i0:i1, j0:j1], mat)
