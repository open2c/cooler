# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
from pandas.api.types import is_categorical
import pandas as pd
import numpy as np
import h5py
import six


def get(grp, lo=0, hi=None, fields=None, convert_enum=True, as_dict=False):
    """
    Query a range of rows from a table as a dataframe.

    A table is an HDF5 group containing equal-length 1D datasets serving as
    columns.

    Parameters
    ----------
    grp : ``h5py.Group`` or any dict-like of array-likes
        Handle to an HDF5 group containing only 1D datasets or any similar
        collection of 1D datasets or arrays
    lo, hi : int, optional
        Range of rows to select from the table.
    fields : str or sequence of str, optional
        Column or list of columns to query. Defaults to all available columns.
        A single string returns a Series instead of a DataFrame.
    convert_enum : bool, optional
        Whether to convert HDF5 enum datasets into ``pandas.Categorical``
        columns instead of plain integer columns. Default is True.
    kwargs : optional
        Options to pass to ``pandas.DataFrame`` or ``pandas.Series``.

    Returns
    -------
    DataFrame or Series

    Notes
    -----
    HDF5 ASCII datasets are converted to Unicode.

    """
    series = False
    if fields is None:
        fields = list(grp.keys())
    elif isinstance(fields, six.string_types):
        fields = [fields]
        series = True

    data = {}
    for field in fields:
        dset = grp[field]

        if convert_enum:
            dt = h5py.check_dtype(enum=dset.dtype)
        else:
            dt = None

        if dt is not None:
            data[field] = pd.Categorical.from_codes(
                dset[lo:hi],
                sorted(dt, key=dt.__getitem__),
                ordered=True)
        elif dset.dtype.type == np.string_:
            data[field] = dset[lo:hi].astype('U')
        else:
            data[field] = dset[lo:hi]

    if as_dict:
        return data

    if data and lo is not None:
        index = np.arange(lo, lo + len(next(iter(data.values()))))
    else:
        index = None

    if series:
        return pd.Series(
            data[fields[0]],
            index=index,
            name=field)
    else:
        return pd.DataFrame(
            data,
            columns=fields,
            index=index)


def put(grp, df, lo=0, store_categories=True, h5opts=None):
    """
    Store a dataframe into a column-oriented table store.

    A table is an HDF5 group containing equal-length 1D datasets serving as
    columns.

    Parameters
    ----------
    h5 : ``h5py.Group``
        Handle to an HDF5 group containing only 1D datasets or any similar
        collection of 1D datasets or arrays
    df : DataFrame or Series
        Data columns to write to the HDF5 group
    lo : int, optional
        Row offset for data to be stored.
    store_categories : bool, optional
        Whether to convert ``pandas.Categorical`` columns into HDF5 enum
        datasets instead of plain integer datasets. Default is True.
    h5opts : dict, optional
        HDF5 dataset filter options to use (compression, shuffling,
        checksumming, etc.). Default is to use autochunking and GZIP
        compression, level 6.

    Notes
    -----
    Categorical data must be ASCII compatible.

    """
    if h5opts is None:
        h5opts = dict(compression='gzip', compression_opts=6)

    if isinstance(df, pd.Series):
        df = df.to_frame()

    fields = df.keys()
    for field, data in six.iteritems(df):

        if np.isscalar(data):
            data = np.array([data])
            dtype = data.dtype
            fillvalue = None
        elif is_categorical(data):
            if store_categories:
                cats = data.cat.categories
                enum = (data.cat.codes.dtype,
                        dict(zip(cats, range(len(cats)))))
                data = data.cat.codes
                dtype = h5py.special_dtype(enum=enum)
                fillvalue = -1
            else:
                data = data.cat.codes
                dtype = data.dtype
                fillvalue = -1
        else:
            data = np.asarray(data)
            if data.dtype in (object, str, bytes):
                dtype = np.dtype('S')
                data = np.array(data, dtype=dtype)
                fillvalue = None
            else:
                dtype = data.dtype
                fillvalue = None

        hi = lo + len(data)
        try:
            dset = grp[field]
        except KeyError:
            dset = grp.create_dataset(
                field,
                shape=(hi,),
                dtype=dtype,
                maxshape=(None,),
                fillvalue=fillvalue,
                **h5opts)
        if hi > len(dset):
            dset.resize((hi,))

        dset[lo:hi] = data


def delete(grp, fields=None):
    """
    Delete columns from a table.

    A table is an HDF5 group containing equal-length 1D datasets serving as
    columns.

    Parameters
    ----------
    grp : ``h5py.Group``
        Handle to an HDF5 group containing only 1D datasets or any similar
        collection of 1D datasets or arrays
    fields : str or sequence of str, optional
        Column or list of columns to query. Defaults to all available columns.
        A single string returns a Series instead of a DataFrame.

    Notes
    -----
    Deleting objects leaves "holes" in HDF5 files and doesn't shrink the file.
    You will need to repack or copy the file contents to reclaim space.
    See the h5repack tool.

    """
    if fields is None:
        fields = list(grp.keys())
    elif isinstance(fields, six.string_types):
        fields = [fields]
    for field in fields:
        if field in grp.keys():
            del grp[field]


def _region_to_extent(h5, chrom_ids, region, binsize):
    chrom, start, end = region
    cid = chrom_ids[chrom]
    if binsize is not None:
        chrom_offset = h5['indexes']['chrom_offset'][cid]
        yield chrom_offset + int(np.floor(start/binsize))
        yield chrom_offset + int(np.ceil(end/binsize))
    else:
        chrom_lo = h5['indexes']['chrom_offset'][cid]
        chrom_hi = h5['indexes']['chrom_offset'][cid + 1]
        chrom_bins = h5['bins']['start'][chrom_lo:chrom_hi]
        yield chrom_lo + np.searchsorted(chrom_bins, start, 'right') - 1
        yield chrom_lo + np.searchsorted(chrom_bins, end, 'left')


def region_to_offset(h5, chrom_ids, region, binsize=None):
    return next(_region_to_extent(h5, chrom_ids, region, binsize))


def region_to_extent(h5, chrom_ids, region, binsize=None):
    return tuple(_region_to_extent(h5, chrom_ids, region, binsize))


class CSRReader(object):
    """
    Retrieves data from a 2D range query on the pixel table of a Cooler.

    Parameters
    ----------
    h5 : ``h5py.File`` or ``h5py.Group``
        Root node of a cooler tree.
    field : str
        Column of the pixel table to query.
    max_chunk : int
        Size of largest chunk to read into memory in a single disk fetch.
        Increase this to increase performance for large queries at the cost of
        memory usage.

    """
    def __init__(self, h5, field, max_chunk):
        self.h5 = h5
        self.field = field
        self.max_chunk = max_chunk

    def index_col(self, i0, i1, j0, j1):
        """Retrieve pixel table row IDs corresponding to query rectangle."""
        edges = self.h5['indexes']['bin1_offset'][i0:i1 + 1]
        index = []
        for lo1, hi1 in zip(edges[:-1], edges[1:]):
            if hi1 - lo1 > 0:
                bin2 = self.h5['pixels']['bin2_id'][lo1:hi1]
                mask = (bin2 >= j0) & (bin2 < j1)
                index.append(lo1 + np.flatnonzero(mask))
        if not index:
            return np.array([], dtype=int)
        else:
            return np.concatenate(index, axis=0)


    def query(self, i0, i1, j0, j1):
        """Retrieve sparse matrix data inside a query rectangle."""
        h5 = self.h5
        field = self.field

        i, j, v = [], [], []
        if (i1 - i0 > 0) or (j1 - j0 > 0):
            edges = h5['indexes']['bin1_offset'][i0:i1 + 1]
            data = h5['pixels'][field]
            p0, p1 = edges[0], edges[-1]

            if (p1 - p0) < self.max_chunk:
                all_bin2 = h5['pixels']['bin2_id'][p0:p1]
                all_data = data[p0:p1]
                dtype = all_bin2.dtype
                for row_id, lo, hi in zip(range(i0, i1),
                                          edges[:-1] - p0,
                                          edges[1:]  - p0):
                    bin2 = all_bin2[lo:hi]
                    mask = (bin2 >= j0) & (bin2 < j1)
                    cols = bin2[mask]

                    i.append(np.full(len(cols), row_id, dtype=dtype))
                    j.append(cols)
                    v.append(all_data[lo:hi][mask])
            else:
                for row_id, lo, hi in zip(range(i0, i1), edges[:-1], edges[1:]):
                    bin2 = h5['pixels']['bin2_id'][lo:hi]
                    mask = (bin2 >= j0) & (bin2 < j1)
                    cols = bin2[mask]
                    dtype = bin2.dtype

                    i.append(np.full(len(cols), row_id, dtype=dtype))
                    j.append(cols)
                    v.append(data[lo:hi][mask])

        if not i:
            i = np.array([], dtype=int)
            j = np.array([], dtype=int)
            v = np.array([])
        else:
            i = np.concatenate(i, axis=0)
            j = np.concatenate(j, axis=0)
            v = np.concatenate(v, axis=0)

        return i, j, v


def _check_bounds(lo, hi, N):
    if hi > N:
        raise IndexError('slice index ({}) out of range'.format(hi))
    if lo < 0:
        raise IndexError('slice index ({}) out of range'.format(lo))


def _comes_before(a0, a1, b0, b1, strict=False):
    if a0 < b0: return a1 <= b0 if strict else a1 <= b1
    return False


def _contains(a0, a1, b0, b1, strict=False):
    if a0 > b0 or a1 < b1: return False
    if strict and (a0 == b0 or a1 == b1): return False
    return a0 <= b0 and a1 >= b1


def query_rect(triu_reader, i0, i1, j0, j1, duplex=True):
    """
    Process a 2D range query on a symmetric matrix using a reader that
    retrieves only upper triangle pixels from the matrix.

    This function is responsible for filling in the missing data in the query
    rectangle and for ensuring that diagonal elements are not duplicated.
    Removing duplicates is important because various sparse matrix constructors
    sum duplicates instead of ignoring them.


    Parameters
    ----------

    triu_reader : callable
        Callable that takes a query rectangle but only returns elements from the
        upper triangle of the parent matrix.

    i0, i1, j0, j1 : int
        Bounding matrix coordinates of the query rectangle. Assumed to be within
        the bounds of the parent matrix.


    Returns
    -------

    i, j, v : 1D arrays


    Details
    -------

    Query cases to consider based on the axes ranges (i0, i1) and (j0, j1):

    1. they are identical
    2. different and non-overlapping
    3. different but partially overlapping
    4. different but one is nested inside the other

    - (1) requires filling in the lower triangle.
    - (3) and (4) require splitting the selection into instances of (1) and (2).

    In some cases, the input axes ranges are swapped to retrieve the data,
    then the final result is transposed.

    """

    # symmetric query
    if (i0, i1) == (j0, j1):
        i, j, v = triu_reader(i0, i1, i0, i1)
        if duplex:
            nodiag = i != j
            i, j, v = np.r_[i, j[nodiag]], np.r_[j, i[nodiag]], np.r_[v, v[nodiag]]

    # asymmetric query
    else:
        transpose = False
        if j0 < i0 or (i0 == j0 and i1 < j1):
            i0, i1, j0, j1 = j0, j1, i0, i1
            transpose = True

        # non-overlapping
        if _comes_before(i0, i1, j0, j1, strict=True):
            i, j, v = triu_reader(i0, i1, j0, j1)

        # partially overlapping
        elif _comes_before(i0, i1, j0, j1):
            ix, jx, vx = triu_reader(i0, j0, j0, i1)
            iy, jy, vy = triu_reader(j0, i1, j0, i1)
            iz, jz, vz = triu_reader(i0, i1, i1, j1)
            if duplex:
                nodiag = iy != jy
                iy, jy, vy = np.r_[iy, jy[nodiag]], np.r_[jy, iy[nodiag]], np.r_[vy, vy[nodiag]]
            i, j, v = np.r_[ix, iy, iz], np.r_[jx, jy, jz], np.r_[vx, vy, vz]

        # nested
        elif _contains(i0, i1, j0, j1):
            ix, jx, vx = triu_reader(i0, j0, j0, j1)
            jy, iy, vy = triu_reader(j0, j1, j0, j1)
            jz, iz, vz = triu_reader(j0, j1, j1, i1)
            if duplex:
                nodiag = iy != jy
                iy, jy, vy = np.r_[iy, jy[nodiag]], np.r_[jy, iy[nodiag]], np.r_[vy, vy[nodiag]]
            i, j, v = np.r_[ix, iy, iz], np.r_[jx, jy, jz], np.r_[vx, vy, vz]

        else:
            raise IndexError("This shouldn't happen")

        if transpose:
            i, j = j, i

    return i, j, v



class _IndexingMixin(object):

    def _unpack_index(self, key):
        if isinstance(key, tuple):
            if len(key) == 2:
                row, col = key
            elif len(key) == 1:
                row, col = key[0], slice(None)
            else:
                raise IndexError('invalid number of indices')
        else:
            row, col = key, slice(None)
        return row, col

    def _isintlike(self, num):
        try:
            int(num)
        except (TypeError, ValueError):
            return False
        return True

    def _process_slice(self, s, nmax):
        if isinstance(s, slice):
            if s.step not in (1, None):
                raise ValueError('slicing with step != 1 not supported')
            i0, i1 = s.start, s.stop
            if i0 is None:
                i0 = 0
            elif i0 < 0:
                i0 = nmax + i0
            if i1 is None:
                i1 = nmax
            elif i1 < 0:
                i1 = nmax + i1
            return i0, i1
        elif self._isintlike(s):
            if s < 0:
                s += nmax
            if s >= nmax:
                raise IndexError('index is out of bounds')
            return int(s), int(s + 1)
        else:
            raise TypeError('expected slice or scalar')


class RangeSelector1D(_IndexingMixin):
    """
    Selector for out-of-core tabular data. Provides DataFrame-like selection of
    columns and list-like access to rows.

    Examples
    --------

    Passing a column name or list of column names as subscript returns a new
    selector.

    >>> sel[ ['A', 'B'] ]  # doctest: +SKIP
    >>> sel['C']

    Passing a scalar or slice as subscript invokes the slicer.

    >>> sel[0]  # doctest: +SKIP
    >>> sel['A'][50:100]

    Calling the fetch method invokes the fetcher to parse the input into an
    integer range and then invokes the slicer.

    >>> sel.fetch('chr3:10,000,000-12,000,000') # doctest: +SKIP
    >>> sel.fetch(('chr3', 10000000, 12000000))

    """
    def __init__(self, fields, slicer, fetcher, nmax):
        self.fields = fields
        self._slice = slicer
        self._fetch = fetcher
        self._shape = (nmax,)

    @property
    def shape(self):
        return self._shape

    @property
    def columns(self):
        return self._slice(self.fields, 0, 0).columns

    @property
    def dtypes(self):
        return self._slice(self.fields, 0, 0).dtypes

    def keys(self):
        return list(self.columns)

    def __len__(self):
        return self._shape[0]

    def __contains__(self, key):
        return key in self.columns

    def __getitem__(self, key):
        # requesting a subset of columns
        if isinstance(key, (list, str)):
            return self.__class__(
                key, self._slice, self._fetch, self._shape[0])

        # requesting an interval of rows
        if isinstance(key, tuple):
            if len(key) == 1:
                key = key[0]
            else:
                raise IndexError('too many indices for table')
        lo, hi = self._process_slice(key, self._shape[0])
        return self._slice(self.fields, lo, hi)

    def fetch(self, *args, **kwargs):
        if self._fetch is not None:
            lo, hi = self._fetch(*args, **kwargs)
            return self._slice(self.fields, lo, hi)
        else:
            raise NotImplementedError


class RangeSelector2D(_IndexingMixin):
    """
    Selector for out-of-core sparse matrix data. Supports 2D scalar and slice
    subscript indexing.

    """
    def __init__(self, field, slicer, fetcher, shape):
        self.field = field
        self._slice = slicer
        self._fetch = fetcher
        self._shape = shape

    @property
    def shape(self):
        return self._shape

    def __len__(self):
        return self._shape[0]

    def __getitem__(self, key):
        s1, s2 = self._unpack_index(key)
        i0, i1 = self._process_slice(s1, self._shape[0])
        j0, j1 = self._process_slice(s2, self._shape[1])
        return self._slice(self.field, i0, i1, j0, j1)

    def fetch(self, *args, **kwargs):
        if self._fetch is not None:
            i0, i1, j0, j1 = self._fetch(*args, **kwargs)
            return self._slice(self.field, i0, i1, j0, j1)
        else:
            raise NotImplementedError
