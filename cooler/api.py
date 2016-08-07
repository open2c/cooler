# -*- coding: utf-8 -*-
from __future__ import division, print_function
import json
import six

from scipy.sparse import coo_matrix
import numpy as np
import pandas
import h5py

from .core import (RangeSelector1D, RangeSelector2D,
                   query_symmetric, query_triu)
from .util import parse_region
from .io import open_hdf5


def get(h5, table_name, lo=0, hi=None, fields=None, **kwargs):
    """
    Query a range of rows from a table as a dataframe.

    A table is an HDF5 group containing equal-length 1D datasets serving as
    columns.

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
    grp = h5[table_name]
    series = False
    if fields is None:
        fields = list(grp.keys())
    elif isinstance(fields, six.string_types):
        fields = [fields]
        series = True

    data = {}
    for field in fields:
        dset = grp[field]
        dt = h5py.check_dtype(enum=dset.dtype)
        if dt is not None:
            data[field] = pandas.Categorical.from_codes(
                dset[lo:hi],
                sorted(dt, key=dt.__getitem__),
                ordered=True)
        elif dset.dtype.type == np.string_:
            data[field] = dset[lo:hi].astype('U')
        else:
            data[field] = dset[lo:hi]

    if data and lo is not None:
        index = np.arange(lo, lo + len(next(iter(data.values()))))
    else:
        index = None

    if series:
        return pandas.Series(
            data[fields[0]],
            index=index,
            name=field,
            **kwargs)
    else:
        return pandas.DataFrame(
            data,
            columns=fields,
            index=index,
            **kwargs)


def _join_bins_and_pixels(h5, df):
    if 'bin2_id' in df:
        bin2 = df['bin2_id']
        bins = get(h5, 'bins', bin2.min(), bin2.max()+1,
                   ['chrom', 'start', 'end'])
        df = (pandas.merge(bins,
                           df,
                           left_index=True,
                           right_on='bin2_id')
                    .drop('bin2_id', axis=1))
    if 'bin1_id' in df:
        bin1 = df['bin1_id']
        bins = get(h5, 'bins', bin1.min(), bin1.max()+1,
                   ['chrom', 'start', 'end'])
        df = (pandas.merge(bins,
                           df,
                           left_index=True,
                           right_on='bin1_id',
                           suffixes=('1', '2'))
                    .drop('bin1_id', axis=1))
    return df


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


class Cooler(object):
    """
    Convenient interface to a cooler HDF5 file.

    * Metadata is accessible as a dictionary through the ``info`` property.

    * Data table range queries are provided through table selectors,
      ``chroms``, ``bins``, and ``pixels``, and are returned as DataFrames
      or Series.

    * Matrix range queries are provided via a matrix selector, ``matrix``,
      and are returned as ``scipy.sparse.coo_matrix`` arrays.

    Parameters
    ----------
    fp : str, h5py.File or h5py.Group
        File path or open handle to the root HDF5 group of a cooler. If ``fp``
        is a file path, the file will be opened temporarily in read-only mode
        when performing operations: such ``Cooler`` objects can be serialized
        and used with multiprocessing, for example.

    """
    def __init__(self, fp):
        self.fp = fp
        with open_hdf5(self.fp) as h5:
            _chromtable = chroms(h5)
            _chromtable['id'] = _chromtable.index
            _chromtable.index = _chromtable['name']
            self._chromlens = _chromtable['length']
            self._chromids = _chromtable['id']
            self._info = info(h5)

    def _get_index(self, name):
        with open_hdf5(self.fp) as h5:
            return h5['indexes'][name][:]

    def offset(self, region):
        """ Bin ID containing the left end of a genomic region

        Parameters
        ----------
        region : str or tuple
            Genomic range

        Returns
        -------
        int

        Examples
        --------
        >>> c.offset('chr3')  # doctest: +SKIP
        1311

        """
        with open_hdf5(self.fp) as h5:
            return region_to_offset(
                h5, self._chromids,
                parse_region(region, self._chromlens))

    def extent(self, region):
        """ Bin IDs containing the left and right ends of a genomic region

        Parameters
        ----------
        region : str or tuple
            Genomic range

        Returns
        -------
        2-tuple of ints

        Examples
        --------
        >>> c.extent('chr3')  # doctest: +SKIP
        (1311, 2131)

        """
        with open_hdf5(self.fp) as h5:
            return region_to_extent(
                h5, self._chromids,
                parse_region(region, self._chromlens))

    @property
    def info(self):
        """ File information and metadata

        Returns
        -------
        dict

        """
        with open_hdf5(self.fp) as h5:
            return info(h5)

    @property
    def shape(self):
        return (self.info['nbins'],) * 2

    def chroms(self):
        """ Chromosome table selector

        Returns
        -------
        Table selector

        """
        def _slice(fields, lo, hi):
            with open_hdf5(self.fp) as h5:
                return chroms(h5, lo, hi, fields)

        return RangeSelector1D(None, _slice, None, self._info['nchroms'])

    def bins(self):
        """ Bin table selector

        Returns
        -------
        Table selector

        """

        def _slice(fields, lo, hi):
            with open_hdf5(self.fp) as h5:
                return bins(h5, lo, hi, fields)

        def _fetch(region):
            with open_hdf5(self.fp) as h5:
                return region_to_extent(h5, self._chromids,
                                        parse_region(region, self._chromlens))

        return RangeSelector1D(None, _slice, _fetch, self._info['nbins'])

    def pixels(self, join=False):
        """ Pixel table selector

        Parameters
        ----------
        join : bool, optional
            Whether to expand bin ID columns into chrom, start, and end columns.
            Default is ``False``.

        Returns
        -------
        Table selector

        """

        def _slice(fields, lo, hi):
            with open_hdf5(self.fp) as h5:
                return pixels(h5, lo, hi, fields, join)

        def _fetch(region):
            with open_hdf5(self.fp) as h5:
                i0, i1 = region_to_extent(
                    h5, self._chromids,
                    parse_region(region, self._chromlens))
                lo = h5['indexes']['bin1_offset'][i0]
                hi = h5['indexes']['bin1_offset'][i1]
                return lo, hi

        return RangeSelector1D(None, _slice, _fetch, self._info['nnz'])

    def matrix(self, field=None, as_pixels=False, join=False, max_query=500000000):
        """ Contact matrix selector

        Parameters
        ----------
        field : str, optional
            Which column of the pixel table to fill the matrix selection with.
            By default, the 'count' column is used.
        as_pixels : bool, optional
            Instead of a complete rectangular sparse matrix, return a DataFrame
            containing the corresponding rows from the pixel table.
            Default is False.
        join : bool, optional
            Whether to expand bin ID columns. False by default. Only applies if
            ``as_pixels`` is True.

        Returns
        -------
        Matrix selector

        """

        def _slice(field, i0, i1, j0, j1):
            with open_hdf5(self.fp) as h5:
                return matrix(h5, i0, i1, j0, j1, field, as_pixels, join, max_query)

        def _fetch(region, region2=None):
            with open_hdf5(self.fp) as h5:
                if region2 is None:
                    region2 = region
                region1 = parse_region(region, self._chromlens)
                region2 = parse_region(region2, self._chromlens)
                i0, i1 = region_to_extent(h5, self._chromids, region1)
                j0, j1 = region_to_extent(h5, self._chromids, region2)
                return i0, i1, j0, j1

        return RangeSelector2D(field, _slice, _fetch, (self._info['nbins'],) * 2)


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
    d = {}
    for k, v in h5.attrs.items():
        if isinstance(v, six.string_types):
            try:
                v = json.loads(v)
            except ValueError:
                pass
        d[k] = v
    return d


def chroms(h5, lo=0, hi=None, fields=None):
    """
    Table describing the chromosomes/scaffolds/contigs used.
    They appear in the same order they occur in the heatmap.

    Parameters
    ----------
    h5 : ``h5py.File`` or ``h5py.Group``
        Open handle to cooler file.
    lo, hi : int, optional
        Range of rows to select from the table.
    fields : sequence of str, optional
        Subset of columns to select from table.

    Returns
    -------
    DataFrame

    """
    if fields is None:
        fields = (pandas.Index(['name', 'length'])
                        .append(pandas.Index(h5['chroms'].keys()))
                        .drop_duplicates())
    return get(h5, 'chroms', lo, hi, fields)


def bins(h5, lo=0, hi=None, fields=None):
    """
    Table describing the genomic bins that make up the axes of the heatmap.

    Parameters
    ----------
    h5 : ``h5py.File`` or ``h5py.Group``
        Open handle to cooler file.
    lo, hi : int, optional
        Range of rows to select from the table.
    fields : sequence of str, optional
        Subset of columns to select from table.

    Returns
    -------
    DataFrame

    """
    if fields is None:
        fields = (pandas.Index(['chrom', 'start', 'end'])
                        .append(pandas.Index(h5['bins'].keys()))
                        .drop_duplicates())
    return get(h5, 'bins', lo, hi, fields)


def pixels(h5, lo=0, hi=None, fields=None, join=True):
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
        fields = (pandas.Index(['bin1_id', 'bin2_id', 'count'])
                        .append(pandas.Index(h5['pixels'].keys()))
                        .drop_duplicates())

    df = get(h5, 'pixels', lo, hi, fields)

    if join:
        df = _join_bins_and_pixels(h5, df)

    return df


def matrix(h5, i0, i1, j0, j1, field=None, as_pixels=False, join=True, max_query=500000000):
    """
    Two-dimensional range query on the Hi-C contact heatmap.
    Returns either a rectangular sparse ``coo_matrix`` or a dataframe of upper
    triangle pixels.

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
    as_pixels: bool, optional
        Return a dataframe of the corresponding rows from the pixel table
        instead of a rectangular sparse matrix. False by default.
    join : bool, optional
        If returning pixels, specifies whether to expand the bin ID columns.
        Has no effect when requesting a rectangular matrix. Default is True.

    Returns
    -------
    coo_matrix (use the ``toarray()`` method to convert to a numpy ``ndarray``.)

    """
    if field is None:
        field = 'count'

    if as_pixels:
        ind, i, j, v = query_triu(h5, field, i0, i1, j0, j1, max_query)
        cols = ['bin1_id', 'bin2_id', field]
        df = pandas.DataFrame(dict(zip(cols, [i, j, v])),
                              columns=cols, index=ind)
        if join:
            df = _join_bins_and_pixels(h5, df)

        return df
    else:
        i, j, v = query_symmetric(h5, field, i0, i1, j0, j1, max_query)
        return coo_matrix((v, (i-i0, j-j0)), (i1-i0, j1-j0))
