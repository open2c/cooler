# -*- coding: utf-8 -*-
from __future__ import division, print_function
import json
import six
import os

from scipy.sparse import coo_matrix
import numpy as np
import pandas
import h5py

from .core import (get, region_to_offset, region_to_extent, RangeSelector1D, 
                   RangeSelector2D, TriuReader, query_rect)
from .util import parse_region, open_hdf5, closing_hdf5
from .io import parse_cooler_uri


class Cooler(object):
    """
    Convenient interface to a Cooler in a COOL file.

    Notes
    -----
    * Metadata is accessible as a dictionary through the ``info`` property.

    * Data table range queries are provided through table selectors:
      ``chroms``, ``bins``, and ``pixels``, which return DataFrames
      or Series.

    * Matrix range queries are provided via a matrix selector, ``matrix``,
      which return NumPy arrays or SciPy sparse ``coo_matrix``.

    """
    def __init__(self, store, root=None, **kwargs):
        """
        Parameters
        ----------
        store : str, h5py.File or h5py.Group
            Path to a COOL file, Cooler URI, or open handle to the root HDF5 
            group of a Cooler.
        root : str, optional
            HDF5 Group path to root of cooler group if ``store`` is a file. 
            This option is deprecated. Instead, use a Cooler URI of the form
            "file_path::group_path".
        kwargs : optional
            Options to be passed to h5py.File() upon every access. By default,
            the file is opened with the default driver and mode='r'.

        Notes
        -----
        If ``store`` is a file path, the file will be opened temporarily in
        when performing operations. This allows ``Cooler`` objects to be
        serialized for multiprocess and distributed computations.
        
        """
        if isinstance(store, six.string_types):
            if root is None:
                self.filename, self.root = parse_cooler_uri(store)
            elif h5py.is_hdf5(store):
                with open_hdf5(store, **kwargs) as h5:
                    self.filename = h5.file.filename
                    self.root = root
            else:
                raise ValueError('Not a valid path to a Cooler file')
            self.uri = self.filename + '::' + self.root
            self.store = self.filename
            self.open_kws = kwargs
        else:
            # Assume an open HDF5 handle, ignore open_kws
            self.filename = store.file.filename
            self.root = store.name
            self.uri = self.filename + '::' + self.root
            self.store = store.file
            self.open_kws = {}
        
        with open_hdf5(self.store, **self.open_kws) as h5:
            grp = h5[self.root]
            _ct = chroms(grp)
            _ct['name'] = _ct['name'].astype(object)
            self._chromsizes = _ct.set_index('name')['length']
            self._chromids = dict(zip(_ct['name'], range(len(_ct))))
            self._info = info(grp)

    def _load_dset(self, path):
        with open_hdf5(self.store, **self.open_kws) as h5:
            grp = h5[self.root]
            return grp[path][:]

    def _load_attrs(self, path):
        with open_hdf5(self.store, **self.open_kws) as h5:
            grp = h5[self.root]
            return dict(grp[path].attrs)

    def open(self, mode='r', **kwargs):
        """ Open the HDF5 group containing the Cooler with h5py

        Functions as a context manager. Any ``open_kws`` passed during 
        construction are ignored.

        Parameters
        ----------
        mode : str
            r (readonly) or r+ (read/write), default: 'r'

        Additional keywords
            See h5py.File

        """
        grp = h5py.File(self.filename, mode, **kwargs)[self.root]
        return closing_hdf5(grp)

    @property
    def binsize(self):
        """ Resolution in base pairs if uniform else None """
        return self._info['bin-size']

    @property
    def chromsizes(self):
        """ Ordered mapping of reference sequences to their lengths in bp """
        return self._chromsizes

    @property
    def chromnames(self):
        """ List of reference sequence names """
        return list(self._chromsizes.index)

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
        with open_hdf5(self.store, **self.open_kws) as h5:
            grp = h5[self.root]
            return region_to_offset(
                grp, self._chromids,
                parse_region(region, self._chromsizes))

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
        with open_hdf5(self.store, **self.open_kws) as h5:
            grp = h5[self.root]
            return region_to_extent(
                grp, self._chromids,
                parse_region(region, self._chromsizes))

    @property
    def info(self):
        """ File information and metadata

        Returns
        -------
        dict

        """
        with open_hdf5(self.store, **self.open_kws) as h5:
            grp = h5[self.root]
            return info(grp)

    @property
    def shape(self):
        return (self._info['nbins'],) * 2

    def chroms(self, **kwargs):
        """ Chromosome table selector

        Returns
        -------
        Table selector

        """
        def _slice(fields, lo, hi):
            with open_hdf5(self.store, **self.open_kws) as h5:
                grp = h5[self.root]
                return chroms(grp, lo, hi, fields, **kwargs)

        return RangeSelector1D(None, _slice, None, self._info['nchroms'])

    def bins(self, **kwargs):
        """ Bin table selector

        Returns
        -------
        Table selector

        """

        def _slice(fields, lo, hi):
            with open_hdf5(self.store, **self.open_kws) as h5:
                grp = h5[self.root]
                return bins(grp, lo, hi, fields, **kwargs)

        def _fetch(region):
            with open_hdf5(self.store, **self.open_kws) as h5:
                grp = h5[self.root]
                return region_to_extent(grp, self._chromids,
                                        parse_region(region, self._chromsizes))

        return RangeSelector1D(None, _slice, _fetch, self._info['nbins'])

    def pixels(self, join=False, **kwargs):
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
            with open_hdf5(self.store, **self.open_kws) as h5:
                grp = h5[self.root]
                return pixels(grp, lo, hi, fields, join, **kwargs)

        def _fetch(region):
            with open_hdf5(self.store, **self.open_kws) as h5:
                grp = h5[self.root]
                i0, i1 = region_to_extent(
                    grp, self._chromids,
                    parse_region(region, self._chromsizes))
                lo = grp['indexes']['bin1_offset'][i0]
                hi = grp['indexes']['bin1_offset'][i1]
                return lo, hi

        return RangeSelector1D(None, _slice, _fetch, self._info['nnz'])

    def matrix(self, field=None, balance=True, sparse=False, as_pixels=False, 
               join=False, ignore_index=True, max_chunk=500000000):
        """ Contact matrix selector

        Parameters
        ----------
        field : str, optional
            Which column of the pixel table to fill the matrix with. By default,
            the 'count' column is used.
        balance : bool, optional
            Whether to apply pre-calculated matrix balancing weights to the
            selection. Default is True.
        sparse: bool, optional
            Return a scipy.sparse.coo_matrix instead of a dense 2D numpy array.
        as_pixels: bool, optional
            Return a DataFrame of the corresponding rows from the pixel table
            instead of a rectangular sparse matrix. False by default.
        join : bool, optional
            If requesting pixels, specifies whether to expand the bin ID columns
            into (chrom, start, end). Has no effect when requesting a 
            rectangular matrix. Default is True.
        ignore_index : bool, optional
            If requesting pixels, don't populate the index column with the pixel
            IDs to improve performance. Default is True.

        Returns
        -------
        Matrix selector

        """

        def _slice(field, i0, i1, j0, j1):
            with open_hdf5(self.store, **self.open_kws) as h5:
                grp = h5[self.root]
                return matrix(grp, i0, i1, j0, j1, field, balance, sparse,
                    as_pixels, join, ignore_index, max_chunk)

        def _fetch(region, region2=None):
            with open_hdf5(self.store, **self.open_kws) as h5:
                grp = h5[self.root]
                if region2 is None:
                    region2 = region
                region1 = parse_region(region, self._chromsizes)
                region2 = parse_region(region2, self._chromsizes)
                i0, i1 = region_to_extent(grp, self._chromids, region1)
                j0, j1 = region_to_extent(grp, self._chromids, region2)
                return i0, i1, j0, j1

        return RangeSelector2D(field, _slice, _fetch, (self._info['nbins'],) * 2)

    def __repr__(self):
        if isinstance(self.store, six.string_types):
            filename = os.path.basename(self.store)
            container = '{}::{}'.format(filename, self.root)
        else:
            container = repr(self.store)
        return '<Cooler "{}">'.format(container)


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


def chroms(h5, lo=0, hi=None, fields=None, **kwargs):
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
    return get(h5['chroms'], lo, hi, fields, **kwargs)


def bins(h5, lo=0, hi=None, fields=None, **kwargs):
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
    out = get(h5['bins'], lo, hi, fields, **kwargs)

    import numbers
    if ('convert_enum' in kwargs and kwargs['convert_enum'] and
            issubclass(out['chrom'].dtype.type, numbers.Integral)):
        chromnames = chroms(h5, fields='name')
        out['chrom'] = pandas.Categorical.from_codes(
            out['chrom'], chromnames, ordered=True)
    return out


def pixels(h5, lo=0, hi=None, fields=None, join=True, **kwargs):
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

    df = get(h5['pixels'], lo, hi, fields, **kwargs)

    if join:
        bins = get(h5['bins'], 0, None, ['chrom', 'start', 'end'], **kwargs)
        df = annotate(df, bins)

    return df


def annotate(pixels, bins, replace=True):
    """
    Add bin annotations to a data frame of pixels.

    This is done by performing a relational "join" against the bin IDs of a 
    table that describes properties of the genomic bins. New columns will be 
    appended on the left of the output data frame.

    Parameters
    ----------
    pixels : DataFrame
        A data frame containing columns named ``bin1_id`` and/or ``bin2_id``.
        If columns ``bin1_id`` and ``bin2_id`` are both present in ``pixels``,
        the adjoined columns will be suffixed with '1' and '2' accordingly.
    bins : DataFrame or DataFrame view
        Data structure that contains a full description of the genomic bins of
        the contact matrix, where the index corresponds to bin IDs.
    replace : bool, optional
        Whether to remove the original ``bin1_id`` and ``bin2_id`` columns from
        the output. Default is True.

    Returns
    -------
    DataFrame

    """
    columns = pixels.columns
    ncols = len(columns)

    if 'bin1_id' in columns:
        if len(bins) > len(pixels):
            bin1 = pixels['bin1_id']
            lo = bin1.min()
            hi = bin1.max() + 1
            lo = 0 if np.isnan(lo) else lo
            hi = 0 if np.isnan(hi) else hi
            right = bins[lo:hi]
        else:
            right = bins[:]

        pixels = pixels.merge(
            right,
            how='left',
            left_on='bin1_id',
            right_index=True)

    if 'bin2_id' in columns:
        if len(bins) > len(pixels):
            bin2 = pixels['bin2_id']
            lo = bin2.min()
            hi = bin2.max() + 1
            lo = 0 if np.isnan(lo) else lo
            hi = 0 if np.isnan(hi) else hi
            right = bins[lo:hi]
        else:
            right = bins[:]

        pixels = pixels.merge(
            right,
            how='left',
            left_on='bin2_id',
            right_index=True,
            suffixes=('1', '2'))

    # rearrange columns
    pixels = pixels[list(pixels.columns[ncols:]) + list(pixels.columns[:ncols])]

    # drop bin IDs
    if replace:
        cols_to_drop = [col for col in ('bin1_id', 'bin2_id') if col in columns]
        pixels = pixels.drop(cols_to_drop, axis=1)

    return pixels


def matrix(h5, i0, i1, j0, j1, field=None, balance=True, sparse=False,
           as_pixels=False, join=True, ignore_index=True, max_chunk=500000000):
    """
    Two-dimensional range query on the Hi-C contact heatmap.
    Depending on the options, returns either a 2D NumPy array, a rectangular
    sparse ``coo_matrix``, or a data frame of upper triangle pixels.

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
    balance : bool, optional
        Whether to apply pre-calculated matrix balancing weights to the
        selection. Default is True.
    sparse: bool, optional
        Return a scipy.sparse.coo_matrix instead of a dense 2D numpy array.
    as_pixels: bool, optional
        Return a DataFrame of the corresponding rows from the pixel table
        instead of a rectangular sparse matrix. False by default.
    join : bool, optional
        If requesting pixels, specifies whether to expand the bin ID columns
        into (chrom, start, end). Has no effect when requesting a rectangular
        matrix. Default is True.
    ignore_index : bool, optional
        If requesting pixels, don't populate the index column with the pixel
        IDs to improve performance. Default is True.

    Returns
    -------
    ndarray, coo_matrix or DataFrame

    Notes
    -----
    Use the ``toarray()`` method to convert to a sparse matrix to a dense
    NumPy array.

    """
    if field is None:
        field = 'count'

    triu_reader = TriuReader(h5, field, max_chunk)

    if isinstance(balance, str):
        name = balance
    elif balance:
        name = 'weight'

    if balance and name not in h5['bins']:
        raise ValueError(
            "No column 'bins/{}' found. Use ``cooler.ice`` to ".format(name) +
            "calculate balancing weights or set balance=False.")

    if as_pixels:
        index = None if ignore_index else triu_reader.index_col(i0, i1, j0, j1)
        i, j, v = triu_reader.query(i0, i1, j0, j1)

        cols = ['bin1_id', 'bin2_id', field]
        df = pandas.DataFrame(dict(zip(cols, [i, j, v])),
                              columns=cols, index=index)

        if balance:
            weights = Cooler(h5).bins()[[name]]
            df2 = annotate(df, weights)
            df['balanced'] = df2[name+'1'] * df2[name+'2'] * df2[field]

        if join:
            bins = Cooler(h5).bins()[['chrom', 'start', 'end']]
            df = annotate(df, bins)

        return df

    elif sparse:
        i, j, v = query_rect(triu_reader.query, i0, i1, j0, j1)
        mat = coo_matrix((v, (i-i0, j-j0)), (i1-i0, j1-j0))

        if balance:
            weights = h5['bins'][name]
            bias1 = weights[i0:i1]
            bias2 = bias1 if (i0, i1) == (j0, j1) else weights[j0:j1]
            mat.data = bias1[mat.row] * bias2[mat.col] * mat.data

        return mat

    else:
        i, j, v = query_rect(triu_reader.query, i0, i1, j0, j1)
        arr = coo_matrix((v, (i-i0, j-j0)), (i1-i0, j1-j0)).toarray()

        if balance:
            weights = h5['bins'][name]
            bias1 = weights[i0:i1]
            bias2 = bias1 if (i0, i1) == (j0, j1) else weights[j0:j1]
            arr = arr * np.outer(bias1, bias2)

        return arr
