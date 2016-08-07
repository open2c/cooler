# -*- coding: utf-8 -*-
from __future__ import division, print_function
import numpy as np


class IndexingMixin(object):

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


class RangeSelector1D(IndexingMixin):
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

    Iterate over the table in chunks of a given size using the slicer.

    >>> for chunk in sel.iterchunks(size=1000):  # doctest: +SKIP
    >>>     ...

    """
    def __init__(self, fields, slicer, fetcher, nmax):
        self.fields = fields
        self._slice = slicer
        self._fetch = fetcher
        self._shape = (nmax,)

    @property
    def shape(self):
        return self._shape

    def __len__(self):
        return self._shape[0]

    def __getitem__(self, key):
        if isinstance(key, (list, str)):
            return self.__class__(key, self._slice, self._fetch, self._shape[0])

        if isinstance(key, tuple):
            if len(key) == 1:
                key = key[0]
            else:
                raise IndexError('too many indices for table')
        lo, hi = self._process_slice(key, self._shape[0])
        return self._slice(self.fields, lo, hi)

    def iterchunks(self, lo=0, hi=None, size=None):
        lo, hi = self._process_slice(slice(lo, hi), self._shape[0])
        if size is None:
            size = hi - lo
        for i in range(lo, hi, size):
            yield self._slice(self.fields, i, i+size)

    def fetch(self, *args, **kwargs):
        if self._fetch is not None:
            lo, hi = self._fetch(*args, **kwargs)
            return self._slice(self.fields, lo, hi)
        else:
            raise NotImplementedError

    # def to_dask(self):
    #     pass


class RangeSelector2D(IndexingMixin):
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


def iter_dataspans(h5, i0, i1, j0, j1):
    if (i1 - i0 > 0) or (j1 - j0 > 0):
        edges = h5['indexes']['bin1_offset'][i0:i1+1]
        for lo1, hi1 in zip(edges[:-1], edges[1:]):
            bin2 = h5['pixels']['bin2_id'][lo1:hi1]
            lo2 = lo1 + np.searchsorted(bin2, j0)
            hi2 = lo1 + np.searchsorted(bin2, j1)
            yield lo2, hi2


def iter_rowspans_with_colmask(h5, i0, i1, j0, j1):
    if (i1 - i0 > 0) or (j1 - j0 > 0):
        edges = h5['indexes']['bin1_offset'][i0:i1+1]
        for lo, hi in zip(edges[:-1], edges[1:]):
            bin2 = h5['pixels']['bin2_id'][lo:hi]
            mask = (bin2 >= j0) & (bin2 < j1)
            yield lo, hi, mask


def slice_triu_csr(h5, field, i0, i1, j0, j1, max_query):
    i, j, v = [], [], []
    if (i1 - i0 > 0) or (j1 - j0 > 0):
        edges = h5['indexes']['bin1_offset'][i0:i1 + 1]
        data = h5['pixels'][field]
        p0, p1 = edges[0], edges[-1]
        ptr = 0
        indptr = [ptr]

        if (p1 - p0) < max_query:
            all_bin2 = h5['pixels']['bin2_id'][p0:p1]
            all_data = data[p0:p1]
            for row_id, lo, hi in zip(range(i0, i1), edges[:-1] - p0, edges[1:] - p0):
                bin2 = all_bin2[lo:hi]
                mask = (bin2 >= j0) & (bin2 < j1)
                cols = bin2[mask]
                ptr += len(v[-1])
                # ind.append(np.arange(p0+lo, p0+hi))
                indptr.append(ptr)
                j.append(cols)
                v.append(all_data[lo:hi][mask])
        else:
            for row_id, lo, hi in zip(range(i0, i1), edges[:-1], edges[1:]):
                bin2 = h5['pixels']['bin2_id'][lo:hi]
                mask = (bin2 >= j0) & (bin2 < j1)
                cols = bin2[mask]
                ptr += len(v[-1])
                # ind.append(np.arange(lo, hi))
                indptr.append(ptr)
                j.append(cols)
                v.append(data[lo:hi][mask])

    indptr = np.array(indptr)
    if not i:
        j = np.array([], dtype=np.int32)
        v = np.array([])
    else:
        j = np.concatenate(j, axis=0)
        v = np.concatenate(v, axis=0)

    return indptr, j, v


def slice_triu_coo(h5, field, i0, i1, j0, j1, max_query):
    i, j, v = [], [], []
    if (i1 - i0 > 0) or (j1 - j0 > 0):
        edges = h5['indexes']['bin1_offset'][i0:i1 + 1]
        data = h5['pixels'][field]
        p0, p1 = edges[0], edges[-1]

        if (p1 - p0) < max_query:
            all_bin2 = h5['pixels']['bin2_id'][p0:p1]
            all_data = data[p0:p1]
            for row_id, lo, hi in zip(range(i0, i1), edges[:-1] - p0, edges[1:] - p0):
                bin2 = all_bin2[lo:hi]
                mask = (bin2 >= j0) & (bin2 < j1)
                cols = bin2[mask]
                # ind.append(np.arange(p0+lo, p0+hi))
                i.append(np.full(len(cols), row_id, dtype=np.int32))
                j.append(cols)
                v.append(all_data[lo:hi][mask])
        else:
            for row_id, lo, hi in zip(range(i0, i1), edges[:-1], edges[1:]):
                bin2 = h5['pixels']['bin2_id'][lo:hi]
                mask = (bin2 >= j0) & (bin2 < j1)
                cols = bin2[mask]
                # ind.append(np.arange(lo, hi))
                i.append(np.full(len(cols), row_id, dtype=np.int32))
                j.append(cols)
                v.append(data[lo:hi][mask])

    if not i:
        i = np.array([], dtype=np.int32)
        j = np.array([], dtype=np.int32)
        v = np.array([])
    else:
        i = np.concatenate(i, axis=0)
        j = np.concatenate(j, axis=0)
        v = np.concatenate(v, axis=0)

    return i, j, v


def query_triu(h5, field, i0, i1, j0, j1, max_query):
    i, j, v = slice_triu_coo(h5, field, i0, i1, j0, j1, max_query)
    edges = h5['indexes']['bin1_offset'][i0:i1 + 1]
    index = np.concatenate([np.arange(lo, hi) for lo, hi in
                            zip(edges[:-1], edges[1:])], axis=0)
    return index, i, j, v


def query_symmetric(h5, field, i0, i1, j0, j1, max_query):
    # Query cases to consider wrt the axes ranges (i0, i1) and (j0 j1):
    # 1. they are identical
    # 2. different and non-overlapping
    # 3. different but partially overlapping
    # 4. different but one inside another other
    #
    # (1) requires filling in the lower triangle.
    # (3) and (4) require splitting the selection into instances of (1) and (2).
    #
    # In some cases, the input axes ranges are swapped to retrieve the data,
    # then the final result is transposed.
    n_bins = h5.attrs['nbins']
    _check_bounds(i0, i1, n_bins)
    _check_bounds(j0, j1, n_bins)

    # symmetric query
    if (i0, i1) == (j0, j1):
        i, j, v = slice_triu_coo(h5, field, i0, i1, i0, i1, max_query)
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
            i, j, v = slice_triu_coo(h5, field, i0, i1, j0, j1, max_query)

        # partially overlapping
        elif _comes_before(i0, i1, j0, j1):
            ix, jx, vx = slice_triu_coo(h5, field, i0, j0, j0, i1, max_query)
            iy, jy, vy = slice_triu_coo(h5, field, j0, i1, j0, i1, max_query)
            iz, jz, vz = slice_triu_coo(h5, field, i0, i1, i1, j1, max_query)
            nodiag = iy != jy
            iy, jy, vy = np.r_[iy, jy[nodiag]], np.r_[jy, iy[nodiag]], np.r_[vy, vy[nodiag]]
            i, j, v = np.r_[ix, iy, iz], np.r_[jx, jy, jz], np.r_[vx, vy, vz]

        # nested
        elif _contains(i0, i1, j0, j1):
            ix, jx, vx = slice_triu_coo(h5, field, i0, j0, j0, j1, max_query)
            iy, jy, vy = slice_triu_coo(h5, field, j0, j1, j0, j1, max_query)
            jz, iz, vz = slice_triu_coo(h5, field, j0, j1, j1, i1, max_query)
            nodiag = iy != jy
            iy, jy, vy = np.r_[iy, jy[nodiag]], np.r_[jy, iy[nodiag]], np.r_[vy, vy[nodiag]]
            i, j, v = np.r_[ix, iy, iz], np.r_[jx, jy, jz], np.r_[vx, vy, vz]

        else:
            raise IndexError("This shouldn't happen")

        if transpose:
            i, j = j, i

    return i, j, v
