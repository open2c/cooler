# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals
import numpy as np


class IndexMixin(object):

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


class Sliceable1D(IndexMixin):
    def __init__(self, slicer, fetcher, nmax):
        self._slice = slicer
        self._fetch = fetcher
        self._shape = (nmax,)
    @property
    def shape(self):
        return self._shape
    def __getitem__(self, key):
        if isinstance(key, tuple):
            if len(key) == 1:
                key = key[0]
            else:
                raise IndexError('too many indices for table')
        lo, hi = self._process_slice(key, self._shape[0])
        return self._slice(lo, hi)
    def fetch(self, *args, **kwargs):
        if self._fetch is not None:
            lo, hi = self._fetch(*args, **kwargs)
            return self._slice(lo, hi)
        else:
            raise NotImplementedError


class Sliceable2D(IndexMixin):
    def __init__(self, slicer, fetcher, shape):
        self._slice = slicer
        self._fetch = fetcher
        self._shape = shape
    @property
    def shape(self):
        return self._shape
    def __getitem__(self, key):
        s1, s2 = self._unpack_index(key)
        i0, i1 = self._process_slice(s1, self._shape[0])
        j0, j1 = self._process_slice(s2, self._shape[1])
        return self._slice(i0, i1, j0, j1)
    def fetch(self, *args, **kwargs):
        if self._fetch is not None:
            i0, i1, j0, j1 = self._fetch(*args, **kwargs)
            return self._slice(i0, i1, j0, j1)
        else:
            raise NotImplementedError


def _region_to_extent(h5, chrom_ids, region, binsize):
    chrom, start, end = region
    chrom_id = chrom_ids.at[chrom]
    if binsize is not None:
        chrom_offset = h5['indexes']['chrom_offset'][chrom_id]
        yield chrom_offset + int(np.floor(start/binsize))
        yield chrom_offset + int(np.ceil(end/binsize))
    else:
        chrom_lo = h5['indexes']['chrom_offset'][chrom_id]
        chrom_hi = h5['indexes']['chrom_offset'][chrom_id + 1]
        chrom_bins = h5['bins']['start'][chrom_lo:chrom_hi]
        yield chrom_lo + np.searchsorted(chrom_bins, start, 'right') - 1
        yield chrom_lo + np.searchsorted(chrom_bins, end, 'left')


def region_to_offset(h5, chrom_ids, region, binsize=None):
    return next(_region_to_extent(h5, chrom_ids, region, binsize))


def region_to_extent(h5, chrom_ids, region, binsize=None):
    return tuple(_region_to_extent(h5, chrom_ids, region, binsize))


def bin1_to_pixel(h5, bin_id):
    return h5['indexes']['bin1_offset'][bin_id]


def iter_dataspans(h5, i0, i1, j0, j1):
    if (i1 - i0 > 0) or (j1 - j0 > 0):
        edges = h5['indexes']['bin1_offset'][i0:i1+1]
        for lo1, hi1 in zip(edges[:-1], edges[1:]):
            bin2 = h5['matrix']['bin2_id'][lo1:hi1]
            lo2 = lo1 + np.searchsorted(bin2, j0)
            hi2 = lo1 + np.searchsorted(bin2, j1)
            yield lo2, hi2


def iter_rowspans_with_colmask(h5, i0, i1, j0, j1):
    if (i1 - i0 > 0) or (j1 - j0 > 0):
        edges = h5['indexes']['bin1_offset'][i0:i1+1]
        for lo, hi in zip(edges[:-1], edges[1:]):
            bin2 = h5['matrix']['bin2_id'][lo:hi]
            mask = (bin2 >= j0) & (bin2 < j1)
            yield lo, hi, mask


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


def slice_triu_as_table(h5, field, i0, i1, j0, j1):
    bin1 = h5['matrix']['bin1_id']
    bin2 = h5['matrix']['bin2_id']
    data = h5['matrix'][field]
    ind, i, j, v = [], [], [], []
    for lo, hi in iter_dataspans(h5, i0, i1, j0, j1):
        ind.append(np.arange(lo, hi))
        i.append(bin1[lo:hi])
        j.append(bin2[lo:hi])
        v.append(data[lo:hi])
    if not i:
        ind = np.array([], dtype=np.int32)
        i = np.array([], dtype=np.int32)
        j = np.array([], dtype=np.int32)
        v = np.array([])
    else:
        ind = np.concatenate(ind, axis=0)
        i = np.concatenate(i, axis=0)
        j = np.concatenate(j, axis=0)
        v = np.concatenate(v, axis=0)
    return ind, i, j, v


def slice_triu_coo(h5, field, i0, i1, j0, j1):
    edges = h5['indexes']['bin1_offset'][i0:i1+1]
    i, j, v = [], [], []
    if (i1 - i0 > 0) or (j1 - j0 > 0):
        edges = h5['indexes']['bin1_offset'][i0:i1+1]
        data = h5['matrix'][field]
        for row_id, lo, hi in zip(range(i0, i1), edges[:-1], edges[1:]):
            bin2 = h5['matrix']['bin2_id'][lo:hi]
            mask = (bin2 >= j0) & (bin2 < j1)
            cols = bin2[mask]
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


def slice_triu_csr(h5, field, i0, i1, j0, j1):
    edges = h5['indexes']['bin1_offset'][i0:i1+1]
    j, v = [], []
    if (i1 - i0 > 0) or (j1 - j0 > 0):
        edges = h5['indexes']['bin1_offset'][i0:i1+1]
        data = h5['matrix'][field]
        ptr = 0
        indptr = [ptr]
        for row_id, lo, hi in zip(range(i0, i1), edges[:-1], edges[1:]):
            bin2 = h5['matrix']['bin2_id'][lo:hi]
            mask = (bin2 >= j0) & (bin2 < j1)
            j.append(bin2[mask])
            v.append(data[lo:hi][mask])
            ptr += len(v[-1])
            indptr.append(ptr)
    indptr = np.array(indptr)
    if not indptr:
        j = np.array([], dtype=np.int32)
        v = np.array([])
    else:
        j = np.concatenate(j, axis=0)
        v = np.concatenate(v, axis=0)
    return indptr, j, v


def slice_matrix(h5, field, i0, i1, j0, j1):
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

    if (i0, i1) == (j0, j1):
        i, j, v = slice_triu_coo(h5, field, i0, i1, i0, i1)
        i, j, v = np.r_[i, j], np.r_[j, i], np.r_[v, v]
    else:
        transpose = False
        if j0 < i0 or (i0 == j0 and i1 < j1):
            i0, i1, j0, j1 = j0, j1, i0, i1
            transpose = True

        if _comes_before(i0, i1, j0, j1, strict=True):
            i, j, v = slice_triu_coo(h5, field, i0, i1, j0, j1)
        elif _comes_before(i0, i1, j0, j1):
            ix, jx, vx = slice_triu_coo(h5, field, i0, j0, j0, i1)
            iy, jy, vy = slice_triu_coo(h5, field, j0, i1, j0, i1)
            iz, jz, vz = slice_triu_coo(h5, field, i0, i1, i1, j1)
            iy, jy, vy = np.r_[iy, jy], np.r_[jy, iy], np.r_[vy, vy]
            i, j, v = np.r_[ix, iy, iz], np.r_[jx, jy, jz], np.r_[vx, vy, vz]
        elif _contains(i0, i1, j0, j1):
            ix, jx, vx = slice_triu_coo(h5, field, i0, j0, j0, j1)
            iy, jy, vy = slice_triu_coo(h5, field, j0, j1, j0, j1)
            jz, iz, vz = slice_triu_coo(h5, field, j0, j1, j1, i1)
            iy, jy, vy = np.r_[iy, jy], np.r_[jy, iy], np.r_[vy, vy]
            i, j, v = np.r_[ix, iy, iz], np.r_[jx, jy, jz], np.r_[vx, vy, vz]
        else:
            raise IndexError("This shouldn't happen")

        if transpose:
            i, j = j, i

    # Remove duplicates coming from main diagonal entries
    # http://stackoverflow.com/questions/28677162/
    # ignoring-duplicate-entries-in-sparse-matrix
    ij = np.c_[i, j]
    idx = np.unique(ij.view(ij.dtype.descr * 2), return_index=True)[1]

    return i[idx], j[idx], v[idx]
