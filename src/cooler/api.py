from __future__ import annotations

import json
import os

import h5py
import numpy as np
import pandas as pd
from pandas.api.types import is_integer_dtype
from scipy.sparse import coo_matrix

from .core import (
    CSRReader,
    DirectRangeQuery2D,
    FillLowerRangeQuery2D,
    RangeSelector1D,
    RangeSelector2D,
    get,
    region_to_extent,
    region_to_offset,
)
from .fileops import list_coolers
from .util import closing_hdf5, open_hdf5, parse_cooler_uri, parse_region

__all__ = ["Cooler", "annotate"]


# The 4DN data portal and hic2cool store these weight vectors in divisive form
_4DN_DIVISIVE_WEIGHTS = {"KR", "VC", "VC_SQRT"}


class Cooler:
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

    Metadata is accessible as a dictionary through the :py:attr:`info`
    property.

    Table selectors, created using :py:meth:`chroms`, :py:meth:`bins`, and
    :py:meth:`pixels`, perform range queries over table rows,
    returning :py:class:`pd.DataFrame` and :py:class:`pd.Series`.

    A matrix selector, created using :py:meth:`matrix`, performs 2D matrix
    range queries, returning :py:class:`numpy.ndarray` or
    :py:class:`scipy.sparse.coo_matrix`.

    """

    def __init__(self, store: str | h5py.Group, root: str | None = None, **kwargs):
        if isinstance(store, str):
            if root is None:
                self.filename, self.root = parse_cooler_uri(store)
            elif h5py.is_hdf5(store):
                with open_hdf5(store, **kwargs) as h5:
                    self.filename = h5.file.filename
                    self.root = root
            else:
                raise ValueError("Not a valid path to a Cooler file")
            self.uri = self.filename + "::" + self.root
            self.store = self.filename
            self.open_kws = kwargs
        else:
            # Assume an open HDF5 handle, ignore open_kws
            self.filename = store.file.filename
            self.root = store.name
            self.uri = self.filename + "::" + self.root
            self.store = store.file
            self.open_kws = {}
        self._refresh()

    def _refresh(self) -> None:
        try:
            with open_hdf5(self.store, **self.open_kws) as h5:
                grp = h5[self.root]
                _ct = chroms(grp)
                _ct["name"] = _ct["name"].astype(object)
                self._chromsizes = _ct.set_index("name")["length"]
                self._chromids = dict(zip(_ct["name"], range(len(_ct))))
                self._info = info(grp)
                mode = self._info.get("storage-mode", "symmetric-upper")
                self._is_symm_upper = mode == "symmetric-upper"
        except KeyError:
            err_msg = f"No cooler found at: {self.store}."
            listing = list_coolers(self.store)
            if len(listing):
                err_msg += (
                    f" Coolers found in {listing}. "
                    + "Use '::' to specify a group path"
                )
            raise KeyError(err_msg) from None

    def _load_dset(self, path: str) -> np.ndarray:
        with open_hdf5(self.store, **self.open_kws) as h5:
            grp = h5[self.root]
            return grp[path][:]

    def _load_attrs(self, path: str) -> dict:
        with open_hdf5(self.store, **self.open_kws) as h5:
            grp = h5[self.root]
            return dict(grp[path].attrs)

    def open(self, mode: str = "r", **kwargs) -> h5py.Group:
        """Open the HDF5 group containing the Cooler with :py:mod:`h5py`

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
    def storage_mode(self) -> str:
        """Indicates whether ordinary sparse matrix encoding is used
        (``"square"``) or whether a symmetric matrix is encoded by storing only
        the upper triangular elements (``"symmetric-upper"``).
        """
        return self._info.get("storage-mode", "symmetric-upper")

    @property
    def binsize(self) -> int | None:
        """Resolution in base pairs if uniform else None"""
        return self._info["bin-size"]

    @property
    def chromsizes(self) -> pd.Series:
        """Ordered mapping of reference sequences to their lengths in bp"""
        return self._chromsizes

    @property
    def chromnames(self) -> list[str]:
        """List of reference sequence names"""
        return list(self._chromsizes.index)

    def offset(self, region: str | tuple[str, int, int]) -> int:
        """Bin ID containing the left end of a genomic region

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
                grp,
                self._chromids,
                parse_region(region, self._chromsizes),
                self.binsize,
            )

    def extent(self, region: str | tuple[str, int, int]) -> tuple[int, int]:
        """Bin IDs containing the left and right ends of a genomic region

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
                grp,
                self._chromids,
                parse_region(region, self._chromsizes),
                self.binsize,
            )

    @property
    def info(self) -> dict:
        """File information and metadata

        Returns
        -------
        dict

        """
        with open_hdf5(self.store, **self.open_kws) as h5:
            grp = h5[self.root]
            return info(grp)

    @property
    def shape(self) -> tuple[int, int]:
        return (self._info["nbins"],) * 2

    def chroms(self, **kwargs) -> RangeSelector1D:
        """Chromosome table selector

        Returns
        -------
        Table selector

        """

        def _slice(fields, lo, hi):
            with open_hdf5(self.store, **self.open_kws) as h5:
                grp = h5[self.root]
                return chroms(grp, lo, hi, fields, **kwargs)

        return RangeSelector1D(None, _slice, None, self._info["nchroms"])

    def bins(self, **kwargs) -> RangeSelector1D:
        """Bin table selector

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
                return region_to_extent(
                    grp,
                    self._chromids,
                    parse_region(region, self._chromsizes),
                    self.binsize,
                )

        return RangeSelector1D(None, _slice, _fetch, self._info["nbins"])

    def pixels(self, join: bool = False, **kwargs) -> RangeSelector1D:
        """Pixel table selector

        Parameters
        ----------
        join : bool, optional
            Whether to expand bin ID columns into chrom, start, and end
            columns. Default is ``False``.

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
                    grp,
                    self._chromids,
                    parse_region(region, self._chromsizes),
                    self.binsize,
                )
                lo = grp["indexes"]["bin1_offset"][i0]
                hi = grp["indexes"]["bin1_offset"][i1]
                return lo, hi

        return RangeSelector1D(None, _slice, _fetch, self._info["nnz"])

    def matrix(
        self,
        field: str | None = None,
        balance: bool | str = True,
        sparse: bool = False,
        as_pixels: bool = False,
        join: bool = False,
        ignore_index: bool = True,
        divisive_weights: bool | None = None,
        chunksize: int = 10000000,
    ) -> RangeSelector2D:
        """Contact matrix selector

        Parameters
        ----------
        field : str, optional
            Which column of the pixel table to fill the matrix with. By
            default, the 'count' column is used.
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
            If requesting pixels, specifies whether to expand the bin ID
            columns into (chrom, start, end). Has no effect when requesting a
            rectangular matrix. Default is True.
        ignore_index : bool, optional
            If requesting pixels, don't populate the index column with the
            pixel IDs to improve performance. Default is True.
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
        lower triangular elements will not be generated. If
        ``as_pixels=False``, those missing non-zero elements will
        automatically be filled in.

        """
        if balance in _4DN_DIVISIVE_WEIGHTS and divisive_weights is None:
            divisive_weights = True

        def _slice(field, i0, i1, j0, j1):
            with open_hdf5(self.store, **self.open_kws) as h5:
                grp = h5[self.root]
                return matrix(
                    grp,
                    i0,
                    i1,
                    j0,
                    j1,
                    field,
                    balance,
                    sparse,
                    as_pixels,
                    join,
                    ignore_index,
                    divisive_weights,
                    chunksize,
                    self._is_symm_upper,
                )

        def _fetch(region, region2=None):
            with open_hdf5(self.store, **self.open_kws) as h5:
                grp = h5[self.root]
                if region2 is None:
                    region2 = region
                region1 = parse_region(region, self._chromsizes)
                region2 = parse_region(region2, self._chromsizes)
                i0, i1 = region_to_extent(grp, self._chromids, region1, self.binsize)
                j0, j1 = region_to_extent(grp, self._chromids, region2, self.binsize)
                return i0, i1, j0, j1

        return RangeSelector2D(field, _slice, _fetch, (self._info["nbins"],) * 2)

    def __repr__(self) -> str:
        if isinstance(self.store, str):
            filename = os.path.basename(self.store)
            container = f"{filename}::{self.root}"
        else:
            container = repr(self.store)
        return f'<Cooler "{container}">'


def info(h5: h5py.Group) -> dict:
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
        if isinstance(v, str):
            try:
                v = json.loads(v)
            except ValueError:
                pass
        d[k] = v
    return d


def chroms(
    h5: h5py.Group,
    lo: int = 0,
    hi: int | None = None,
    fields: list[str] | None = None,
    **kwargs,
) -> pd.DataFrame:
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
        fields = (
            pd.Index(["name", "length"])
            .append(pd.Index(h5["chroms"].keys()))
            .drop_duplicates()
        )
    return get(h5["chroms"], lo, hi, fields, **kwargs)


def bins(
    h5: h5py.Group,
    lo: int = 0,
    hi: int | None = None,
    fields: list[str] | None = None,
    **kwargs,
) -> pd.DataFrame:
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
        fields = (
            pd.Index(["chrom", "start", "end"])
            .append(pd.Index(h5["bins"].keys()))
            .drop_duplicates()
        )

    # If convert_enum is not explicitly set to False, chrom IDs will get
    # converted to categorical chromosome names, provided the ENUM header
    # exists in bins/chrom. Otherwise, they will return as integers.
    out = get(h5["bins"], lo, hi, fields, **kwargs)

    # Handle the case where the ENUM header doesn't exist but we want to
    # convert integer chrom IDs to categorical chromosome names.
    if "chrom" in fields:
        convert_enum = kwargs.get("convert_enum", True)
        if isinstance(fields, str):
            chrom_col = out
        else:
            chrom_col = out["chrom"]

        if is_integer_dtype(chrom_col.dtype) and convert_enum:
            chromnames = chroms(h5, fields="name")
            chrom_col = pd.Categorical.from_codes(chrom_col, chromnames, ordered=True)
            if isinstance(fields, str):
                out = pd.Series(chrom_col, out.index)
            else:
                out["chrom"] = chrom_col

    return out


def pixels(
    h5: h5py.Group,
    lo: int = 0,
    hi: int | None = None,
    fields: list[str] | None = None,
    join: bool = True,
    **kwargs,
) -> pd.DataFrame:
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
        fields = (
            pd.Index(["bin1_id", "bin2_id"])
            .append(pd.Index(h5["pixels"].keys()))
            .drop_duplicates()
        )

    df = get(h5["pixels"], lo, hi, fields, **kwargs)

    if join:
        bins = get(h5["bins"], 0, None, ["chrom", "start", "end"], **kwargs)
        df = annotate(df, bins, replace=True)

    return df


def annotate(
    pixels: pd.DataFrame, bins: pd.DataFrame | RangeSelector1D, replace: bool = False
) -> pd.DataFrame:
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

    # End-inclusive slicer for bins
    if isinstance(bins, RangeSelector1D):

        def _loc_slice(sel, beg, end):
            # slicing a range selector is end-exclusive like iloc
            return sel[beg : end + 1 if end is not None else None]

    else:

        def _loc_slice(df, beg, end):
            # loc slicing a dataframe is end-inclusive
            return df.loc[beg:end]

    # Extract the required bin ranges from the bin table.
    # NOTE: Bin IDs in the pixel table may be uint. Avoid using these for
    # indexing - they can easily get cast to float and cause problems.
    anns = []

    # Select bin annotations that correspond to the bin1 IDs in the pixels df
    if "bin1_id" in columns:
        bin1 = pixels["bin1_id"].to_numpy().astype(np.int64, copy=False, casting="safe")
        if len(bin1) == 0:
            bmin = bmax = 0
        elif len(bins) > len(pixels):
            bmin, bmax = bin1.min(), bin1.max()
        else:
            bmin, bmax = 0, None
        ann1 = _loc_slice(bins, bmin, bmax)
        anns.append(
            ann1.iloc[bin1 - ann1.index[0]]
            .rename(columns=lambda x: x + "1")
            .reset_index(drop=True)
        )

    # Select bin annotations that correspond to the bin2 IDs in the pixels df
    if "bin2_id" in columns:
        bin2 = pixels["bin2_id"].to_numpy().astype(np.int64, copy=False, casting="safe")
        if len(bin2) == 0:
            bmin = bmax = 0
        elif len(bins) > len(pixels):
            bmin, bmax = bin2.min(), bin2.max()
        else:
            bmin, bmax = 0, None
        ann2 = _loc_slice(bins, bmin, bmax)
        anns.append(
            ann2.iloc[bin2 - ann2.index[0]]
            .rename(columns=lambda x: x + "2")
            .reset_index(drop=True)
        )

    # Drop original bin IDs if not wanted
    if replace:
        cols_to_drop = [col for col in ("bin1_id", "bin2_id") if col in columns]
        pixels = pixels.drop(cols_to_drop, axis=1)

    # Concatenate bin annotations with pixels
    out = pd.concat([*anns, pixels.reset_index(drop=True)], axis=1)
    out.index = pixels.index
    return out


def matrix(
    h5: h5py.Group,
    i0: int,
    i1: int,
    j0: int,
    j1: int,
    field: str | None = None,
    balance: bool | str = True,
    sparse: bool = False,
    as_pixels: bool = False,
    join: bool = True,
    ignore_index: bool = True,
    divisive_weights: bool = False,
    chunksize: int = 10000000,
    fill_lower: bool = True,
) -> np.ndarray | coo_matrix | pd.DataFrame:
    """
    Two-dimensional range query on the Hi-C contact heatmap.
    Depending on the options, returns either a 2D NumPy array, a rectangular
    sparse ``coo_matrix``, or a data frame of pixels.

    Parameters
    ----------
    h5 : :py:class:`h5py.File` or :py:class:`h5py.Group`
        Open handle to cooler file.
    i0, i1 : int, optional
        Bin range along the 0th (row) axis of the heatmap.
    j0, j1 : int, optional
        Bin range along the 1st (col) axis of the heatmap.
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
        field = "count"

    if isinstance(balance, str):
        name = balance
    elif balance:
        name = "weight"

    if balance and name not in h5["bins"]:
        raise ValueError(
            f"No column 'bins/{name}'"
            + "found. Use ``cooler.balance_cooler`` to "
            + "calculate balancing weights or set balance=False."
        )

    reader = CSRReader(h5["pixels"], h5["indexes/bin1_offset"][:])

    if as_pixels:
        # The historical behavior for as_pixels is to return only explicitly stored
        # pixels so we ignore the ``fill_lower`` parameter in this case.
        engine = DirectRangeQuery2D(
            reader, field, (i0, i1, j0, j1), chunksize, return_index=not ignore_index
        )
        df = engine.to_frame()

        if balance:
            weights = Cooler(h5).bins()[[name]]
            df2 = annotate(df, weights, replace=False)
            if divisive_weights:
                df2[name + "1"] = 1 / df2[name + "1"]
                df2[name + "2"] = 1 / df2[name + "2"]
            df["balanced"] = df2[name + "1"] * df2[name + "2"] * df2[field]

        if join:
            bins = Cooler(h5).bins()[["chrom", "start", "end"]]
            df = annotate(df, bins, replace=True)

        return df

    elif sparse:
        if fill_lower:
            engine = FillLowerRangeQuery2D(reader, field, (i0, i1, j0, j1), chunksize)
        else:
            engine = DirectRangeQuery2D(reader, field, (i0, i1, j0, j1), chunksize)
        mat = engine.to_sparse_matrix()

        if balance:
            weights = h5["bins"][name]
            bias1 = weights[i0:i1]
            bias2 = bias1 if (i0, i1) == (j0, j1) else weights[j0:j1]
            if divisive_weights:
                bias1 = 1 / bias1
                bias2 = 1 / bias2
            mat.data = bias1[mat.row] * bias2[mat.col] * mat.data

        return mat

    else:
        if fill_lower:
            engine = FillLowerRangeQuery2D(reader, field, (i0, i1, j0, j1), chunksize)
        else:
            engine = DirectRangeQuery2D(reader, field, (i0, i1, j0, j1), chunksize)
        arr = engine.to_array()

        if balance:
            weights = h5["bins"][name]
            bias1 = weights[i0:i1]
            bias2 = bias1 if (i0, i1) == (j0, j1) else weights[j0:j1]
            if divisive_weights:
                bias1 = 1 / bias1
                bias2 = 1 / bias2
            arr = arr * np.outer(bias1, bias2)

        return arr
