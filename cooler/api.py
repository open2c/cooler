# -*- coding: utf-8 -*-
from __future__ import division, print_function
import json
import six
import os

from pandas.api.types import is_integer_dtype
from scipy.sparse import coo_matrix
import numpy as np
import pandas as pd
import h5py

from .core import (get, region_to_offset, region_to_extent, RangeSelector1D,
                   RangeSelector2D, CSRReader, query_rect)
from .util import parse_cooler_uri, parse_region, open_hdf5, closing_hdf5
from .fileops import list_coolers


__all__ = ['Cooler', 'annotate']


# The 4DN data portal and hic2cool store these weight vectors in divisive form
_4DN_DIVISIVE_WEIGHTS = {'KR', 'VC', 'VC_SQRT'}


class Cooler(object):
    """
    A convenient interface to a cooler data collection.

    Parameters
    ----------
    store : str, :py:class:`h5py.File` or :py:class:`h5py.Group`
        Path to a cooler file, URI string, or open handle to the root HDF5
        group of a cooler data collection.
    root : str, optional [deprecated]
        HDF5 Group path to root of cooler group if ``store`` is a file.
        This option is deprecated. Instead, use a URI string of the form
        :file:`<file_path>::<group_path>`.
    kwargs : optional
        Options to be passed to :py:class:`h5py.File()` upon every access.
        By default, the file is opened with the default driver and mode='r'.

    Notes
    -----
    If ``store`` is a file path, the file will be opened temporarily in
    when performing operations. This allows :py:class:`Cooler` objects to be
    serialized for multiprocess and distributed computations.

    Metadata is accessible as a dictionary through the :py:attr:`info` property.

    Table selectors, created using :py:meth:`chroms`, :py:meth:`bins`, and
    :py:meth:`pixels`, perform range queries over table rows,
    returning :py:class:`pd.DataFrame` and :py:class:`pd.Series`.

    A matrix selector, created using :py:meth:`matrix`, performs 2D matrix
    range queries, returning :py:class:`numpy.ndarray` or
    :py:class:`scipy.sparse.coo_matrix`.

    """
    def __init__(self, store, root=None, **kwargs):
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
        self._refresh()

    def _refresh(self):
        try:
            with open_hdf5(self.store, **self.open_kws) as h5:
                grp = h5[self.root]
                _ct = chroms(grp)
                _ct['name'] = _ct['name'].astype(object)
                self._chromsizes = _ct.set_index('name')['length']
                self._chromids = dict(zip(_ct['name'], range(len(_ct))))
                self._info = info(grp)
                mode = self._info.get('storage-mode', u"symmetric-upper")
                self._is_symm_upper = mode == u"symmetric-upper"
        except KeyError:
            err_msg = "No cooler found at: {}.".format(self.store)
            listing = list_coolers(self.store)
            if len(listing):
                err_msg += (" Coolers found in {}. ".format(listing) +
                            "Use '::' to specify a group path")
            raise KeyError(err_msg)

    def _load_dset(self, path):
        with open_hdf5(self.store, **self.open_kws) as h5:
            grp = h5[self.root]
            return grp[path][:]

    def _load_attrs(self, path):
        with open_hdf5(self.store, **self.open_kws) as h5:
            grp = h5[self.root]
            return dict(grp[path].attrs)

    def open(self, mode='r', **kwargs):
        """ Open the HDF5 group containing the Cooler with :py:mod:`h5py`

        Functions as a context manager. Any ``open_kws`` passed during
        construction are ignored.

        Parameters
        ----------
        mode : str, optional [default: 'r']
            * ``'r'`` (readonly)
            * ``'r+'`` or ``'a'`` (read/write)

        Notes
        -----
            For other parameters, see :py:class:`h5py.File`.

        """
        grp = h5py.File(self.filename, mode, **kwargs)[self.root]
        return closing_hdf5(grp)

    @property
    def storage_mode(self):
        """Indicates whether ordinary sparse matrix encoding is used
        (``"square"``) or whether a symmetric matrix is encoded by storing only
        the upper triangular elements (``"symmetric-upper"``).
        """
        return self._info.get('storage-mode', u"symmetric-upper")

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
               join=False, ignore_index=True, divisive_weights=None,
               max_chunk=500000000):
        """ Contact matrix selector

        Parameters
        ----------
        field : str, optional
            Which column of the pixel table to fill the matrix with. By default,
            the 'count' column is used.
        balance : bool, optional
            Whether to apply pre-calculated matrix balancing weights to the
            selection. Default is True and uses a column named 'weight'.
            Alternatively, pass the name of the bin table column containing
            the desired balancing weights. Set to False to return untransformed
            counts.
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
        divisive_weights : bool, optional
            Force balancing weights to be interpreted as divisive (True) or
            multiplicative (False). Weights are always assumed to be
            multiplicative by default unless named KR, VC or SQRT_VC, in which
            case they are assumed to be divisive by default.

        Returns
        -------
        Matrix selector

        Notes
        -----
        If ``as_pixels=True``, only data explicitly stored in the pixel table
        will be returned: if the cooler's storage mode is symmetric-upper,
        lower triangular elements will not be generated. If ``as_pixels=False``,
        those missing non-zero elements will automatically be filled in.

        """
        if balance in _4DN_DIVISIVE_WEIGHTS and divisive_weights is None:
            divisive_weights = True

        def _slice(field, i0, i1, j0, j1):
            with open_hdf5(self.store, **self.open_kws) as h5:
                grp = h5[self.root]
                return matrix(grp, i0, i1, j0, j1, field, balance, sparse,
                    as_pixels, join, ignore_index, divisive_weights, max_chunk,
                    self._is_symm_upper)

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
    h5 : :py:class:`h5py.File` or :py:class:`h5py.Group`
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
    h5 : :py:class:`h5py.File` or :py:class:`h5py.Group`
        Open handle to cooler file.
    lo, hi : int, optional
        Range of rows to select from the table.
    fields : sequence of str, optional
        Subset of columns to select from table.

    Returns
    -------
    :py:class:`DataFrame`

    """
    if fields is None:
        fields = (pd.Index(['name', 'length'])
                        .append(pd.Index(h5['chroms'].keys()))
                        .drop_duplicates())
    return get(h5['chroms'], lo, hi, fields, **kwargs)


def bins(h5, lo=0, hi=None, fields=None, **kwargs):
    """
    Table describing the genomic bins that make up the axes of the heatmap.

    Parameters
    ----------
    h5 : :py:class:`h5py.File` or :py:class:`h5py.Group`
        Open handle to cooler file.
    lo, hi : int, optional
        Range of rows to select from the table.
    fields : sequence of str, optional
        Subset of columns to select from table.

    Returns
    -------
    :py:class:`DataFrame`

    """
    if fields is None:
        fields = (pd.Index(['chrom', 'start', 'end'])
                        .append(pd.Index(h5['bins'].keys()))
                        .drop_duplicates())

    # If convert_enum is not explicitly set to False, chrom IDs will get
    # converted to categorical chromosome names, provided the ENUM header
    # exists in bins/chrom. Otherwise, they will return as integers.
    out = get(h5['bins'], lo, hi, fields, **kwargs)

    # Handle the case where the ENUM header doesn't exist but we want to
    # convert integer chrom IDs to categorical chromosome names.
    if (is_integer_dtype(out['chrom'].dtype)
            and kwargs.get('convert_enum', True)):
        chromnames = chroms(h5, fields='name')
        out['chrom'] = pd.Categorical.from_codes(
            out['chrom'], chromnames, ordered=True)

    return out


def pixels(h5, lo=0, hi=None, fields=None, join=True, **kwargs):
    """
    Table describing the nonzero upper triangular pixels of the Hi-C contact
    heatmap.

    Parameters
    ----------
    h5 : :py:class:`h5py.File` or :py:class:`h5py.Group`
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
    :py:class:`DataFrame`

    """
    if fields is None:
        fields = (pd.Index(['bin1_id', 'bin2_id'])
                        .append(pd.Index(h5['pixels'].keys()))
                        .drop_duplicates())

    df = get(h5['pixels'], lo, hi, fields, **kwargs)

    if join:
        bins = get(h5['bins'], 0, None, ['chrom', 'start', 'end'], **kwargs)
        df = annotate(df, bins, replace=True)

    return df


def annotate(pixels, bins, replace=False):
    """
    Add bin annotations to a data frame of pixels.

    This is done by performing a relational "join" against the bin IDs of a
    table that describes properties of the genomic bins. New columns will be
    appended on the left of the output data frame.

    .. versionchanged:: 0.8.0
       The default value of ``replace`` changed to False.

    Parameters
    ----------
    pixels : :py:class:`DataFrame`
        A data frame containing columns named ``bin1_id`` and/or ``bin2_id``.
        If columns ``bin1_id`` and ``bin2_id`` are both present in ``pixels``,
        the adjoined columns will be suffixed with '1' and '2' accordingly.
    bins : :py:class:`DataFrame` or DataFrame selector
        Data structure that contains a full description of the genomic bins of
        the contact matrix, where the index corresponds to bin IDs.
    replace : bool, optional
        Remove the original ``bin1_id`` and ``bin2_id`` columns from the
        output. Default is False.

    Returns
    -------
    :py:class:`DataFrame`

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
           as_pixels=False, join=True, ignore_index=True, divisive_weights=False,
           max_chunk=500000000, is_upper=True):
    """
    Two-dimensional range query on the Hi-C contact heatmap.
    Depending on the options, returns either a 2D NumPy array, a rectangular
    sparse ``coo_matrix``, or a data frame of pixels.

    Parameters
    ----------
    h5 : :py:class:`h5py.File` or :py:class:`h5py.Group`
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
        selection. Default is True and uses a column named 'weight'.
        Alternatively, pass the name of the bin table column containing the
        desired balancing weights. Set to False to return untransformed counts.
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
    If ``as_pixels=True``, only data explicitly stored in the pixel table
    will be returned: if the cooler's storage mode is symmetric-upper,
    lower triangular elements will not be generated. If ``as_pixels=False``,
    those missing non-zero elements will automatically be filled in.

    """
    if field is None:
        field = 'count'

    if isinstance(balance, str):
        name = balance
    elif balance:
        name = 'weight'

    if balance and name not in h5['bins']:
        raise ValueError(
            "No column 'bins/{}' found. Use ``cooler.balance_cooler`` to ".format(name) +
            "calculate balancing weights or set balance=False.")

    if as_pixels:
        reader = CSRReader(h5, field, max_chunk)
        index = None if ignore_index else reader.index_col(i0, i1, j0, j1)
        i, j, v = reader.query(i0, i1, j0, j1)

        cols = ['bin1_id', 'bin2_id', field]
        df = pd.DataFrame(dict(zip(cols, [i, j, v])),
                              columns=cols, index=index)

        if balance:
            weights = Cooler(h5).bins()[[name]]
            df2 = annotate(df, weights, replace=False)
            if divisive_weights:
                df2[name+'1'] = 1 / df2[name+'1']
                df2[name+'2'] = 1 / df2[name+'2']
            df['balanced'] = df2[name+'1'] * df2[name+'2'] * df2[field]

        if join:
            bins = Cooler(h5).bins()[['chrom', 'start', 'end']]
            df = annotate(df, bins, replace=True)

        return df

    elif sparse:
        reader = CSRReader(h5, field, max_chunk)
        if is_upper:
            i, j, v = query_rect(reader.query, i0, i1, j0, j1, duplex=True)
        else:
            i, j, v = reader.query(i0, i1, j0, j1)
        mat = coo_matrix((v, (i-i0, j-j0)), (i1-i0, j1-j0))

        if balance:
            weights = h5['bins'][name]
            bias1 = weights[i0:i1]
            bias2 = bias1 if (i0, i1) == (j0, j1) else weights[j0:j1]
            if divisive_weights:
                bias1 = 1 / bias1
                bias2 = 1 / bias2
            mat.data = bias1[mat.row] * bias2[mat.col] * mat.data

        return mat

    else:
        reader = CSRReader(h5, field, max_chunk)
        if is_upper:
            i, j, v = query_rect(reader.query, i0, i1, j0, j1, duplex=True)
        else:
            i, j, v = reader.query(i0, i1, j0, j1)
        arr = coo_matrix((v, (i-i0, j-j0)), (i1-i0, j1-j0)).toarray()

        if balance:
            weights = h5['bins'][name]
            bias1 = weights[i0:i1]
            bias2 = bias1 if (i0, i1) == (j0, j1) else weights[j0:j1]
            if divisive_weights:
                bias1 = 1 / bias1
                bias2 = 1 / bias2
            arr = arr * np.outer(bias1, bias2)

        return arr
