from io import BytesIO

import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from cooler import core
from cooler.core._rangequery import _comes_before, _contains
from cooler.core._selectors import _IndexingMixin


def make_hdf5_table(mode):
    s = BytesIO()
    f = h5py.File(s, mode)
    h5opts = dict(compression='gzip', compression_opts=6, maxshape=(None,))
    grp = f.create_group('table')
    grp.create_dataset(
        'chrom',
        data=np.array(['chr1', 'chr1', 'chr1', 'chr2', 'chr2'], dtype='S'),
        **h5opts
    )
    grp.create_dataset(
        'start',
        data=[0, 10, 20, 0, 10],
        **h5opts
    )
    grp.create_dataset(
        'end',
        data=[10, 20, 32, 10, 21],
        **h5opts
    )
    grp.create_dataset(
        'value',
        data=[1.1, 2.0, 3.0, 4.0, 5.0],
        **h5opts
    )
    f.flush()
    return f


def test_get():
    f = make_hdf5_table('a')
    out = core.get(f['table'], 0, 3, ['chrom', 'value'])
    assert isinstance(out, pd.DataFrame)
    assert len(out.columns) == 2
    assert out['chrom'].astype('U').tolist() == ['chr1', 'chr1', 'chr1']
    assert np.allclose(out['value'].values, [1.1, 2.0, 3.0])

    out = core.get(f['table'], 0, 3, 'value')
    assert isinstance(out, pd.Series)
    assert np.allclose(out.values, [1.1, 2.0, 3.0])

    out = core.get(f['table'], 0, 3, 'value', as_dict=True)
    assert isinstance(out, dict)
    assert np.allclose(out['value'], [1.1, 2.0, 3.0])

    out = core.get(f['table'])
    assert len(out) == 5
    assert len(out.columns) == 4

    out = core.get(f['table'], lo=None)
    assert len(out) == 5
    assert len(out.columns) == 4

    out = core.get(f['table'], lo=3)
    assert len(out) == 2
    assert len(out.columns) == 4


def test_put():
    f = make_hdf5_table('a')

    # append
    df = pd.DataFrame({
        'chrom': ['chr3', 'chr3'],
        'start': [0, 20],
        'end': [20, 40],
        'value': [4.0, 5.0],
    })
    core.put(f['table'], df, lo=5)
    f.flush()
    out = core.get(f['table'])
    assert len(out) == 7

    # insert a categorical column
    s = pd.Series(pd.Categorical(out['chrom'], ordered=True), index=out.index)
    s.name = 'chrom_enum'
    core.put(f['table'], s)
    assert h5py.check_dtype(enum=f['table/chrom_enum'].dtype)
    out = core.get(f['table'])
    assert len(out.columns) == 5
    assert pd.api.types.is_categorical_dtype(out['chrom_enum'].dtype)
    out = core.get(f['table'], convert_enum=False)
    assert len(out.columns) == 5
    assert pd.api.types.is_integer_dtype(out['chrom_enum'].dtype)

    # don't convert categorical to enum
    s.name = 'chrom_string'
    core.put(f['table'], s, store_categories=False)
    out = core.get(f['table'])
    assert len(out.columns) == 6
    assert not pd.api.types.is_categorical_dtype(out['chrom_string'].dtype)

    # scalar input
    core.put(f['table'], {'foo': 42})
    out = core.get(f['table'])
    assert len(out.columns) == 7
    assert (out['foo'] == 42).all()


def test_delete():
    f = make_hdf5_table('a')
    core.delete(f['table'])
    assert len(f['table'].keys()) == 0

    f = make_hdf5_table('a')
    core.delete(f['table'], ['chrom'])
    assert len(f['table'].keys()) == 3

    f = make_hdf5_table('a')
    core.delete(f['table'], 'chrom')
    assert len(f['table'].keys()) == 3


def test_region_to_offset_extent(mock_cooler):
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
    assert core.region_to_offset(
        mock_cooler, chromID_lookup, region, binsize
    ) == first
    assert core.region_to_offset(
        mock_cooler, chromID_lookup, region, None
    ) == first

    region = ("chr1", 159, 400)
    first, last = 1, 3
    assert core.region_to_extent(
        mock_cooler, chromID_lookup, region, binsize
    ) == (first, last + 1)
    assert core.region_to_extent(mock_cooler, chromID_lookup, region, None) == (
        first,
        last + 1,
    )
    assert core.region_to_offset(
        mock_cooler, chromID_lookup, region, binsize
    ) == first
    assert core.region_to_offset(
        mock_cooler, chromID_lookup, region, None
    ) == first


def test_interval_ops():
    assert _comes_before(1, 5, 6, 10)
    assert not _comes_before(6, 10, 1, 5)
    assert _comes_before(1, 5, 6, 10, strict=True)
    assert _comes_before(1, 5, 5, 10, strict=True)
    assert _comes_before(1, 5, 3, 10)
    assert not _comes_before(1, 5, 3, 10, strict=True)

    assert _contains(1, 10, 3, 5)
    assert _contains(1, 10, 3, 5, strict=True)
    assert _contains(1, 10, 3, 10)
    assert not _contains(1, 10, 3, 10, strict=True)
    assert not _contains(1, 5, 6, 10)


def test_indexing_mixin():

    class Impl(_IndexingMixin):
        def __init__(self, shape):
            self._shape = shape

        def __getitem__(self, key):
            s1, s2 = self._unpack_index(key)
            i0, i1 = self._process_slice(s1, self._shape[0])
            j0, j1 = self._process_slice(s2, self._shape[1])
            return i0, i1, j0, j1

    obj = Impl((10, 10))

    # row scalar
    assert obj[5] == (5, 6, 0, 10)
    assert obj[5, ] == (5, 6, 0, 10)

    # row slice
    assert obj[:] == (0, 10, 0, 10)
    assert obj[1:5] == (1, 5, 0, 10)
    assert obj[:-2] == (0, 8, 0, 10)
    assert obj[-2:] == (8, 10, 0, 10)

    # slice + scalar
    assert obj[1:5, 3] == (1, 5, 3, 4)
    assert obj[2, 1:5] == (2, 3, 1, 5)
    assert obj[2, 0:-2] == (2, 3, 0, 8)
    assert obj[-2, 0:-2] == (8, 9, 0, 8)

    # row + col scalar query
    assert obj[5, 5] == (5, 6, 5, 6)

    # row + col slices
    assert obj[:, :] == (0, 10, 0, 10)
    assert obj[1:5, :] == (1, 5, 0, 10)
    assert obj[:, 2:3] == (0, 10, 2, 3)
    assert obj[1:5, 2:3] == (1, 5, 2, 3)

    with pytest.raises(IndexError):
        obj[10]

    with pytest.raises(TypeError):
        obj[{}]

    # with pytest.raises(TypeError):
    #     obj[4.5]


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


    slicer = lambda fields, lo, hi: pd.DataFrame(  # noqa
        np.zeros((hi - lo, len(fields))),
        columns=fields
    )
    fetcher = lambda x: list(map(int, x.split(':')))  # noqa
    nmax = 50
    sel = core.RangeSelector1D(['a', 'b', 'c'], slicer, fetcher, nmax)
    assert sel.columns.tolist() == ['a', 'b', 'c']
    assert list(sel.keys()) == ['a', 'b', 'c']
    assert isinstance(sel.dtypes, pd.Series)
    assert 'a' in sel
    assert len(sel) == 50
    assert len(sel[['a', 'b']].columns) == 2
    assert len(sel[['a']].columns) == 1
    assert np.all(sel[5] == 0)
    assert np.all(sel[5, ] == 0)
    assert len(sel.fetch('5:10')) == 5

    # some things are broken here
    series_view = sel['a']
    assert len(series_view) == 50
    assert series_view.shape == (50,)
    # series_view.columns ???


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


    slicer = lambda field, i0, i1, j0, j1: ( # noqa
        np.zeros((i1 - i0, j1 - j0))
    )
    fetcher = lambda x, y=None: (0, 10, 0, 10)  # noqa
    nmax = 50
    sel = core.RangeSelector2D('count', slicer, fetcher, (nmax, nmax))
    assert sel.shape == (50, 50)
    assert len(sel) == 50
    assert sel[:10, 5:10].shape == (10, 5)
    assert sel.fetch('0:10', '0:10').shape == (10, 10)


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
        r = sparse.coo_matrix(
            (
                (
                    mock_cooler["pixels/count"],
                    (mock_cooler["pixels/bin1_id"], mock_cooler["pixels/bin2_id"]),
                )
            ),
            (mock_cooler.attrs["nbins"],) * 2,
        )
        r_triu = r.toarray()
        r_fill = r.toarray() + r.toarray().T

        reader = core.CSRReader(
            mock_cooler["pixels"],
            mock_cooler["indexes"]["bin1_offset"]
        )

        # query of data in storage (upper triangle)
        query = core.DirectRangeQuery2D(
            reader,
            field="count",
            bbox=(i0, i1, j0, j1),
            chunksize=10,
            return_index=True
        )
        arr_triu = query.to_array()
        assert np.allclose(r_triu[i0:i1, j0:j1], arr_triu)

        # query with filled-in lower triangle
        query = core.FillLowerRangeQuery2D(
            reader,
            field="count",
            bbox=(i0, i1, j0, j1),
            chunksize=10,
            return_index=True
        )
        arr_fill = query.to_array()
        assert np.allclose(r_fill[i0:i1, j0:j1], arr_fill)


def test_csr_reader():
    pass


def test_query_rect():
    pass
