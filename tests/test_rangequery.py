"""
Property-based tests of 2D range queries against dense NumPy references.

Random dense matrices are stored as CSR pixel tables (upper triangle or full
square), queried with random bounding boxes and chunk sizes, and compared with
plain NumPy slicing of the dense matrix.
"""

import os.path as op
import sys
import tempfile

import hypothesis.strategies as st
import numpy as np
import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis.extra.numpy import arrays
from scipy import sparse

import cooler
from cooler import core
from cooler.core._rangequery import concat, frame_slice_from_dict

MAX_BINS = 25


@st.composite
def dense_matrices(draw, symmetric):
    n = draw(st.integers(1, MAX_BINS))
    # Mostly zeros, so the pixel tables are sparse
    mat = draw(arrays(np.int32, (n, n), elements=st.sampled_from([0, 0, 0, 1, 2, 7])))
    if symmetric:
        mat = np.triu(mat) + np.triu(mat, 1).T
    return mat


@st.composite
def bboxes(draw, n):
    i0 = draw(st.integers(0, n))
    i1 = draw(st.integers(i0, n))
    j0 = draw(st.integers(0, n))
    j1 = draw(st.integers(j0, n))
    return i0, i1, j0, j1


def csr_pixels(mat, upper):
    """Pixel table (dict of columns) and bin1 offsets for a dense matrix."""
    stored = np.triu(mat) if upper else mat
    # np.nonzero returns coordinates in row-major order, i.e. CSR order
    i, j = np.nonzero(stored)
    pixels = {
        "bin1_id": i.astype(np.int64),
        "bin2_id": j.astype(np.int64),
        "count": stored[i, j],
    }
    offsets = np.searchsorted(i, np.arange(len(mat) + 1))
    return pixels, offsets


def check_result(query, pixels, expected):
    i0, i1, j0, j1 = query.bbox
    result = query.get()

    # Every pixel is inside the bounding box and reported exactly once
    bin1, bin2 = result["bin1_id"], result["bin2_id"]
    assert np.all((bin1 >= i0) & (bin1 < i1) & (bin2 >= j0) & (bin2 < j1))
    assert len(set(zip(bin1.tolist(), bin2.tolist()))) == len(bin1)

    # Values match the dense reference, for dense and sparse outputs
    np.testing.assert_array_equal(query.to_array(), expected)
    assert query.to_array().shape == (i1 - i0, j1 - j0)
    arr = query.to_sparse_array()
    assert isinstance(arr, sparse.coo_array)
    np.testing.assert_array_equal(arr.toarray(), expected)
    mat = query.to_sparse_matrix()
    assert isinstance(mat, sparse.coo_matrix)
    np.testing.assert_array_equal(mat.toarray(), expected)

    # The index points back at the stored pixel each value came from
    if query.return_index:
        np.testing.assert_array_equal(
            pixels["count"][result["__index"]], result["count"]
        )


@settings(max_examples=200, deadline=None)
@given(data=st.data(), chunksize=st.integers(1, 50))
def test_fill_lower_query_symmetric_upper(data, chunksize):
    mat = data.draw(dense_matrices(symmetric=True))
    bbox = data.draw(bboxes(len(mat)))
    pixels, offsets = csr_pixels(mat, upper=True)
    reader = core.CSRReader(pixels, offsets)

    i0, i1, j0, j1 = bbox
    query = core.FillLowerRangeQuery2D(
        reader, "count", bbox, chunksize, return_index=True
    )
    check_result(query, pixels, mat[i0:i1, j0:j1])


@settings(max_examples=200, deadline=None)
@given(data=st.data(), chunksize=st.integers(1, 50))
def test_direct_query_symmetric_upper(data, chunksize):
    mat = data.draw(dense_matrices(symmetric=True))
    bbox = data.draw(bboxes(len(mat)))
    pixels, offsets = csr_pixels(mat, upper=True)
    reader = core.CSRReader(pixels, offsets)

    # Only what is stored: the upper triangle
    i0, i1, j0, j1 = bbox
    query = core.DirectRangeQuery2D(reader, "count", bbox, chunksize, return_index=True)
    check_result(query, pixels, np.triu(mat)[i0:i1, j0:j1])


@settings(max_examples=200, deadline=None)
@given(data=st.data(), chunksize=st.integers(1, 50))
def test_direct_query_square(data, chunksize):
    mat = data.draw(dense_matrices(symmetric=False))
    bbox = data.draw(bboxes(len(mat)))
    pixels, offsets = csr_pixels(mat, upper=False)
    reader = core.CSRReader(pixels, offsets)

    i0, i1, j0, j1 = bbox
    query = core.DirectRangeQuery2D(reader, "count", bbox, chunksize, return_index=True)
    check_result(query, pixels, mat[i0:i1, j0:j1])


@settings(max_examples=50, deadline=None)
@given(data=st.data(), symmetric=st.booleans())
def test_matrix_selector(data, symmetric):
    mat = data.draw(dense_matrices(symmetric=symmetric))
    n = len(mat)
    pixels, _ = csr_pixels(mat, upper=symmetric)
    bins = pd.DataFrame(
        {"chrom": "chr1", "start": np.arange(n) * 10, "end": np.arange(1, n + 1) * 10}
    )

    with tempfile.TemporaryDirectory() as tmp:
        path = op.join(tmp, "test.cool")
        cooler.create_cooler(
            path, bins, pd.DataFrame(pixels), symmetric_upper=symmetric
        )
        clr = cooler.Cooler(path)

        for _ in range(5):
            i0, i1, j0, j1 = data.draw(bboxes(n))
            expected = mat[i0:i1, j0:j1]

            sel = clr.matrix(balance=False)
            np.testing.assert_array_equal(sel[i0:i1, j0:j1], expected)

            sel = clr.matrix(balance=False, sparse=True)
            np.testing.assert_array_equal(sel[i0:i1, j0:j1].toarray(), expected)

            # Pixel output only includes stored pixels (no filled-in lower
            # triangle), each exactly once
            stored = np.triu(mat) if symmetric else mat
            df = clr.matrix(balance=False, as_pixels=True)[i0:i1, j0:j1]
            assert not df.duplicated(["bin1_id", "bin2_id"]).any()
            dense = np.zeros((i1 - i0, j1 - j0), dtype=mat.dtype)
            dense[df["bin1_id"] - i0, df["bin2_id"] - j0] = df["count"]
            np.testing.assert_array_equal(dense, stored[i0:i1, j0:j1])


@settings(max_examples=100, deadline=None)
@given(data=st.data(), chunksize=st.integers(1, 20))
def test_chunked_and_dask_outputs(data, chunksize):
    mat = data.draw(dense_matrices(symmetric=True))
    bbox = data.draw(bboxes(len(mat)))
    pixels, offsets = csr_pixels(mat, upper=True)
    reader = core.CSRReader(pixels, offsets)
    query = core.FillLowerRangeQuery2D(reader, "count", bbox, chunksize)
    expected = query.to_frame().sort_values(["bin1_id", "bin2_id"])

    def same_pixels(df):
        df = df.sort_values(["bin1_id", "bin2_id"])
        np.testing.assert_array_equal(df.to_numpy(), expected.to_numpy())

    # The chunks partition the result
    chunks = [query.get_chunk(i) for i in range(query.n_chunks)]
    assert len(chunks) == query.n_chunks == len(query.tasks)
    if chunks:
        same_pixels(frame_slice_from_dict(concat(*chunks), "count"))
    else:
        assert expected.empty

    # Dask outputs match the eager result
    dask = pytest.importorskip("dask")
    (delayed_chunks,) = dask.compute(query.to_delayed())
    assert len(delayed_chunks) == query.n_chunks
    same_pixels(query.to_dask_frame().compute())


def test_csr_reader_direct_calls():
    mat = np.array([[1, 2, 0], [2, 0, 3], [0, 3, 4]], dtype=np.int32)
    pixels, offsets = csr_pixels(mat, upper=True)
    reader = core.CSRReader(pixels, offsets)
    bbox = (0, 3, 0, 3)

    # Without a row span, the reader covers all rows of the bounding box
    full = reader("count", bbox)
    explicit = reader("count", bbox, row_span=(0, 3))
    for key in explicit:
        np.testing.assert_array_equal(full[key], explicit[key])

    # An empty row span returns empty columns with the pixel table's dtypes
    empty = reader("count", bbox, row_span=(1, 1), return_index=True)
    assert set(empty) == {"bin1_id", "bin2_id", "count", "__index"}
    assert all(len(col) == 0 for col in empty.values())
    assert empty["count"].dtype == pixels["count"].dtype

    # Chunks are indexed from 0 to n_chunks - 1
    query = core.DirectRangeQuery2D(reader, "count", bbox, chunksize=1)
    query.get_chunk(query.n_chunks - 1)
    with pytest.raises(IndexError):
        query.get_chunk(query.n_chunks)
    with pytest.raises(IndexError):
        query.get_chunk(-1)


def test_dask_output_needs_dask(monkeypatch):
    path = op.join(op.dirname(__file__), "data", "toy.symm.upper.2.cool")
    clr = cooler.Cooler(path)
    with clr.open("r") as h5:
        reader = core.CSRReader(h5["pixels"], h5["indexes/bin1_offset"][:])
        query = core.DirectRangeQuery2D(reader, "count", (0, 4, 0, 4), 10)
        monkeypatch.setitem(sys.modules, "dask", None)
        monkeypatch.setitem(sys.modules, "dask.dataframe", None)
        with pytest.raises(ImportError, match="dask"):
            query.to_delayed()
        with pytest.raises(ImportError, match="dask"):
            query.to_dask_frame()
