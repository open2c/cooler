# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals
import json

from scipy.sparse import coo_matrix
import numpy as np
import pandas

from .models import (Sliceable1D, Sliceable2D, slice_matrix, 
                     region_to_offset, region_to_extent)
from .util import parse_region
from .io import open_hdf5


def get(h5, table_name, lo=0, hi=None, fields=None, **kwargs):
    """
    Fetch the raw columns of a table.

    Parameters
    ----------
    h5 : ``h5py.File`` or ``h5py.Group``
        Open handle to cooler file.
    table_name : str
        Name of HDF5 Group.
    lo, hi : int, optional
        Range of rows to select from the table.
    fields : sequence of str, optional
        Selection of columns to query. Defaults to all available columns.
    kwargs : optional
        Options to pass to ``pandas.DataFrame``.

    Returns
    -------
    DataFrame

    """
    if fields is None:
        fields = list(h5[table_name].keys())
    data = {field: h5[table_name][field][lo:hi] for field in fields}
    if lo is not None:
        index = np.arange(lo, lo + len(next(iter(data.values()))))
    else:
        index = None
    return pandas.DataFrame(
        data,
        columns=fields,
        index=np.arange(lo, hi),
        **kwargs)


def info(h5):
    """
    File and user metadata dict.

    Parameters
    ----------
    h5 : ``h5py.File`` or ``h5py.Group``
        Open handle to cooler file.

    Returns
    -------
    dict

    """
    d = dict(h5.attrs.items())
    d['metadata'] = json.loads(d.get('metadata', '{}'))
    return d


def chromtable(h5, lo=0, hi=None):
    """
    Table describing the chromosomes/scaffolds/contigs used.
    They appear in the same order they occur in the heatmap.

    Parameters
    ----------
    h5 : ``h5py.File`` or ``h5py.Group``
        Open handle to cooler file.
    lo, hi : int, optional
        Range of rows to select from the table.

    Returns
    -------
    DataFrame

    """
    names = h5['scaffolds']['name'][lo:hi].astype('U')
    lengths = h5['scaffolds']['length'][lo:hi]
    if lo is not None:
        index = np.arange(lo, lo+len(names))
    else:
        index = None
    return pandas.DataFrame({
            'name': names,
            'length': lengths,
        }, columns=['name', 'length'],
           index=index)


def bintable(h5, lo=0, hi=None):
    """
    Table describing the genomic bins that make up the axes of the heatmap.

    Parameters
    ----------
    h5 : ``h5py.File`` or ``h5py.Group``
        Open handle to cooler file.
    lo, hi : int, optional
        Range of rows to select from the table.

    Returns
    -------
    DataFrame

    """
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


def pixeltable(h5, lo=0, hi=None, fields=None, join=True):
    """
    Table describing the nonzero upper triangular pixels of the Hi-C contact
    heatmap.

    Parameters
    ----------
    h5 : ``h5py.File`` or ``h5py.Group``
        Open handle to cooler file.
    lo, hi : int, optional
        Range of rows to select from the table.
    fields : sequence of str, optional
        Subset of columns to select from table.
    join : bool, optional
        Whether or not to expand bin ID columns to their full bin description 
        (chrom, start, end). Default is True.

    Returns
    -------
    DataFrame

    """
    if fields is None:
        fields = set(h5['matrix'].keys())
        fields.remove('bin1_id')
        fields.remove('bin2_id')

    bin1 = h5['matrix']['bin1_id'][lo:hi]
    bin2 = h5['matrix']['bin2_id'][lo:hi]
    if lo is not None:
        index = np.arange(lo, lo+len(bin1))
    else:
        index = None
    data = {
        'bin1_id': bin1,
        'bin2_id': bin2,
    }
    data.update({field: h5['matrix'][field][lo:hi] for field in fields})

    df = pandas.DataFrame(
        data,
        columns=['bin1_id', 'bin2_id'] + list(fields),
        index=index)

    if join:
        bins = bintable(h5, bin2.min(), bin2.max()+1)
        df = (pandas.merge(bins, 
                           df, 
                           left_index=True,
                           right_on='bin2_id')
                    .drop('bin2_id', axis=1))
        bins = bintable(h5, bin1.min(), bin1.max()+1)
        df = (pandas.merge(bins,
                           df, 
                           left_index=True, 
                           right_on='bin1_id', 
                           suffixes=('1', '2'))
                    .drop('bin1_id', axis=1))

    return df


def matrix(h5, i0, i1, j0, j1, field=None):
    """
    Range query on the Hi-C contact heatmap.

    Parameters
    ----------
    h5 : ``h5py.File`` or ``h5py.Group``
        Open handle to cooler file.
    i0, i1 : int, optional
        Bin range along the 0th (row) axis of the heatap.
    j0, j1 : int, optional
        Bin range along the 1st (col) axis of the heatap.
    field : str, optional
        Which column of the pixel table to fill the matrix with. By default,
        the 'count' column is used.

    Returns
    -------
    coo_matrix (use the ``toarray()`` method to convert to a numpy ``ndarray``.)

    """
    if field is None: field = 'count'
    i, j, v = slice_matrix(h5, field, i0, i1, j0, j1)
    return coo_matrix((v, (i-i0, j-j0)), (i1-i0, j1-j0))


class Cooler(object):
    def __init__(self, fp):
        self.fp = fp
        with open_hdf5(self.fp) as h5:
            self._chromtable = chromtable(h5)
            self._chromtable['id'] = self._chromtable.index
            self._chromtable.index = self._chromtable['name']
            self._info = info(h5)

    def offset(self, region):
        with open_hdf5(self.fp) as h5:
            return region_to_offset(h5, self._chromtable,
                parse_region(region, self._chromtable['length']))

    def extent(self, region):
        with open_hdf5(self.fp) as h5:
            return region_to_offset(h5, self._chromtable,
                parse_region(region, self._chromtable['length']))

    @property
    def info(self):
        with open_hdf5(self.fp) as h5:
            return info(h5)

    def get(self, table_name, fields=None):
        with open_hdf5(self.fp) as h5:
            nmax = h5[table_name][next(iter(h5[table_name].keys()))]
        def _slice(lo, hi):
            with open_hdf5(self.fp) as h5:
                return get(h5, table_name, lo, hi, fields)
        return Sliceable1D(_slice, None, nmax)

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

    def pixeltable(self, fields=None, join=True):
        def _slice(lo, hi):
            with open_hdf5(self.fp) as h5:
                return pixeltable(h5, lo, hi, fields, join)
        def _fetch(region):
            with open_hdf5(self.fp) as h5:
                i0, i1 = region_to_extent(h5, self._chromtable,
                    parse_region(region, self._chromtable['length']))
                lo = h5['indexes']['bin1_offset'][i0]
                hi = h5['indexes']['bin1_offset'][i1]
                return lo, hi
        return Sliceable1D(_slice, _fetch, self._info['nnz'])

    def matrix(self, field=None):
        def _slice(i0, i1, j0, j1):
            with open_hdf5(self.fp) as h5:
                return matrix(h5, i0, i1, j0, j1, field)
        def _fetch(region, region2=None):
            with open_hdf5(self.fp) as h5:
                if region2 is None:
                    region2 = region
                region1 = parse_region(region, self._chromtable['length'])
                region2 = parse_region(region2, self._chromtable['length'])
                i0, i1 = region_to_extent(h5, self._chromtable, region1)
                j0, j1 = region_to_extent(h5, self._chromtable, region2)
                return i0, i1, j0, j1
        return Sliceable2D(_slice, _fetch, (self._info['nbins'],) * 2)
