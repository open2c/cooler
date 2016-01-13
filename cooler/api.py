# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals
from scipy.sparse import coo_matrix
import numpy as np
import pandas

from .models import (Sliceable1D, Sliceable2D, slice_matrix, 
                     region_to_offset, region_to_extent)
from .util import parse_region
from .io import open_hdf5


def get(h5, table_name, fields=None, lo=None, hi=None):
    if fields is None:
        fields = h5[table_name].keys()
    return pandas.DataFrame(
        {field: h5[table_name][field][lo:hi] for field in fields},
        columns=fields,
        index=np.arange(lo, hi),
        **kwargs)


def info(h5):
    return dict(h5.attrs.items())


def chromtable(h5, lo=None, hi=None):
    names = h5['scaffolds']['name'][lo:hi].astype('U')
    lengths = h5['scaffolds']['length'][lo:hi]
    if lo is not None:
        index = np.arange(lo, lo+len(names))
    else:
        index = None
    return pandas.DataFrame({
            'name': names,
            'id': np.arange(len(names)),
            'length': lengths,
        }, columns=['name', 'id', 'length'],
           index=index)


def bintable(h5, lo=None, hi=None):
    chrom_ids = h5['bins']['chrom_id'][lo:hi]
    names = h5['scaffolds']['name'][:].astype('U')
    chroms = names[chrom_ids]
    starts = h5['bins']['start'][lo:hi]
    ends = h5['bins']['end'][lo:hi]
    if lo is not None:
        index = np.arange(lo, lo+len(chroms))
    else:
        index = None
    return pandas.DataFrame({
            'chrom': chroms,
            'start': starts,
            'end': ends,
        }, columns=['chrom', 'start', 'end'],
           index=index)


def pixeltable(h5, field, lo=None, hi=None, join=True):
    bin1 = h5['matrix']['bin1_id'][lo:hi]
    bin2 = h5['matrix']['bin2_id'][lo:hi]
    values = h5['matrix'][field][lo:hi]
    if lo is not None:
        index = np.arange(lo, lo+len(bin1))
    else:
        index = None
    df = pandas.DataFrame({
            'bin1_id': bin1,
            'bin2_id': bin2,
            field: values,
        },
        columns=['bin1_id', 'bin2_id', field],
        index=index)

    if join:
        bins = bintable(h5, lo, hi)
        df = (pandas.merge(bins, 
                           df, 
                           left_index=True,
                           right_on='bin2_id')
                    .drop('bin2_id', axis=1))
        df = (pandas.merge(bins, 
                           df, 
                           left_index=True, 
                           right_on='bin1_id', 
                           suffixes=('1', '2'))
                    .drop('bin1_id', axis=1))

    return df


def matrix(h5, field, i0, i1, j0, j1):
    i, j, v = slice_matrix(h5, field, i0, i1, j0, j1)
    return coo_matrix((v, (i-i0, j-j0)), (i1-i0, j1-j0))


class Cooler(object):
    def __init__(self, fp):
        self.fp = fp
        with open_hdf5(self.fp) as h5:
            self._chromtable = chromtable(h5)
            self._chromtable.index = self._chromtable['name']
            self._info = info(h5)

    def offset(self, region):
        with open_hdf5(self.fp) as h5:
            return region_to_offset(h5, self._chromtable, parse_region(region))

    def extent(self, region):
        with open_hdf5(self.fp) as h5:
            return region_to_offset(h5, self._chromtable, parse_region(region))

    @property
    def info(self):
        with open_hdf5(self.fp) as h5:
            return info(h5)

    def get(self, table_name, fields=None):
        if fields is None:
            with open_hdf5(self.fp) as h5:
                return h5[table_name].keys()
        def _slice(lo, hi):
            with open_hdf5(self.fp) as h5:
                return get(h5, table_name, fields, lo, hi)
        return Sliceable1D(_slice, None)

    def chromtable(self):
        def _slice(lo, hi):
            with open_hdf5(self.fp) as h5:
                return chromtable(h5, lo, hi)
        return Sliceable1D(_slice, None, self._info['nchroms'])
        
    def bintable(self):
        def _slice(lo, hi):
            with open_hdf5(self.fp) as h5:
                return bintable(h5, lo, hi)
        def _fetch(region):
            with open_hdf5(self.fp) as h5:
                return region_to_extent(h5, self._chromtable,
                                        parse_region(region))
        return Sliceable1D(_slice, _fetch, self._info['nbins'])

    def pixeltable(self, join=True):
        def _slice(lo, hi):
            with open_hdf5(self.fp) as h5:
                return pixeltable(h5, 'count', lo, hi, join)
        def _fetch(region):
            with open_hdf5(self.fp) as h5:
                lo, hi = tuple(
                    map(bin_to_pixel, 
                        region_to_extent(h5, self._chromtable, 
                                         parse_region(region))))
                return lo, hi
        return Sliceable1D(_slice, _fetch, self._info['nnz'])

    def matrix(self, field):
        def _slice(i0, i1, j0, j1):
            with open_hdf5(self.fp) as h5:
                return matrix(h5, field, i0, i1, j0, j1)
        def _fetch(region, region2=None):
            with open_hdf5(self.fp) as h5:
                if region2 is None:
                    region2 = region
                region1, region2 = parse_region(region), parse_region(region2)
                i0, i1 = region_to_extent(h5, self._chromtable, region1)
                j0, j1 = region_to_extent(h5, self._chromtable, region2)
                return i0, i1, j0, j1
        return Sliceable2D(_slice, _fetch, (self._info['nbins'],) * 2)
