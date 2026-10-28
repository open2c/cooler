from __future__ import annotations

import sys
import warnings
from collections.abc import Iterable
from typing import Any

import h5py
import numpy as np
import pandas as pd

from .._logging import get_logger
from .._typing import Tabular
from ..core import get
from ..util import (
    infer_meta,
)
from . import (
    CHROM_DTYPE,
    CHROMID_DTYPE,
    PIXEL_FIELDS,
)
from ._ingest import BadInputError
from ._writer import (
    _CoolerWriter,
    _pixel_schema,
    _ScoolWriter,
    _set_h5opts,
    _UnorderedCoolerWriter,
    _value_aggs,
)

logger = get_logger("cooler.create")


def _rename_chroms(grp: h5py.Group, rename_dict: dict[str, str], h5opts: dict) -> None:
    chroms = get(grp["chroms"]).set_index("name")
    n_chroms = len(chroms)
    new_names = np.array(
        chroms.rename(rename_dict).index.values, dtype=CHROM_DTYPE
    )  # auto-adjusts char length

    # Replace chroms/name
    del grp["chroms/name"]
    grp["chroms"].create_dataset(
        "name", shape=(n_chroms,), dtype=new_names.dtype, data=new_names, **h5opts
    )

    # Replace the bins/chroms enum mapping if applicable
    bins = get(grp["bins"])
    n_bins = len(bins)
    if isinstance(bins["chrom"].dtype, pd.CategoricalDtype):
        idmap = dict(zip(new_names, range(n_chroms)))
        chrom_ids = bins["chrom"].cat.codes
        chrom_dtype = h5py.special_dtype(enum=(CHROMID_DTYPE, idmap))
        del grp["bins/chrom"]
        try:
            grp["bins"].create_dataset(
                "chrom", shape=(n_bins,), dtype=chrom_dtype, data=chrom_ids, **h5opts
            )
        except ValueError:
            # If HDF5 enum header would be too large,
            # try storing chrom IDs as raw int instead
            chrom_dtype = CHROMID_DTYPE
            grp["bins"].create_dataset(
                "chrom", shape=(n_bins,), dtype=chrom_dtype, data=chrom_ids, **h5opts
            )


def rename_chroms(clr, rename_dict: dict[str, str], h5opts: dict | None = None) -> None:
    """
    Substitute existing chromosome/contig names for new ones. They will be
    written to the file and the Cooler object will be refreshed.

    Parameters
    ----------
    clr : Cooler
        Cooler object that can be opened with write permissions.
    rename_dict : dict
        Dictionary of old -> new chromosome names. Any names omitted from
        the dictionary will be kept as is.
    h5opts : dict, optional
        HDF5 filter options.

    """
    h5opts = _set_h5opts(h5opts)

    with clr.open("r+") as f:
        _rename_chroms(f, rename_dict, h5opts)
    clr._refresh()


def _get_dtypes_arg(
    dtypes: pd.Series | dict | None, kwargs: dict[str, Any]
) -> pd.Series | dict | None:
    if "dtype" in kwargs:
        if dtypes is None:
            dtypes = kwargs.pop("dtype")
            warnings.warn("Use dtypes= instead of dtype=", FutureWarning, stacklevel=2)
        else:
            raise ValueError(
                'Received both "dtypes" and "dtype" arguments. '
                'Please use "dtypes" to provide a column name -> dtype mapping. '
                '"dtype" remains as an alias but is deprecated.'
            )
    return dtypes


def create(
    cool_uri,
    bins,
    pixels,
    columns=None,
    dtypes=None,
    metadata=None,
    assembly=None,
    symmetric_upper=True,
    mode=None,
    h5opts=None,
    boundscheck=True,
    triucheck=None,
    dupcheck=True,
    ensure_sorted=False,
    lock=None,
    binsize=None,
    **kwargs,
):
    """
    Create a new Cooler.

    ``binsize`` is the bin size of the bin table, if known. By default, it is
    inferred from the bin table, which isn't possible when every chromosome is
    a single bin.

    Deprecated parameters
    ---------------------
    chromsizes : Series
        Chromsizes are now inferred from ``bins``.
    dtype : dict, optional
        Dictionary mapping column names in the pixel table to dtypes.
        Use the ``dtypes`` argument instead.

    """
    if mode is None:
        mode = "w"

    if not isinstance(bins, pd.DataFrame):
        raise ValueError(  # noqa: TRY004
            "Second positional argument must be a pandas DataFrame. "
            "Note that the `chromsizes` argument is now deprecated: "
            "see documentation for `create`."
        )
    dtypes = _get_dtypes_arg(dtypes, kwargs)
    if kwargs:
        raise TypeError(
            f"create() got unexpected keyword arguments: {', '.join(sorted(kwargs))}"
        )

    # Determine the appropriate iterable. A dask DataFrame means dask is
    # already imported, so don't import it for every other input.
    dd = sys.modules.get("dask.dataframe")
    if dd is not None and isinstance(pixels, dd.DataFrame):
        iterable = (x.compute() for x in pixels.to_delayed())
        input_columns = infer_meta(pixels).columns
    elif isinstance(pixels, pd.DataFrame):
        iterable = (pixels,)
        input_columns = infer_meta(pixels).columns
    elif isinstance(pixels, dict):
        iterable = (pixels,)
        input_columns = infer_meta([(k, v.dtype) for (k, v) in pixels.items()]).columns
    else:
        iterable = pixels
        input_columns = None

    # If possible, ensure all expected columns are available before writing
    if input_columns is not None:
        for col in _pixel_schema(columns, dtypes)[0]:
            if col not in input_columns:
                col_type = "Standard" if col in PIXEL_FIELDS else "User"
                raise ValueError(f"{col_type} column not found in input: '{col}'")

    writer = _CoolerWriter(
        cool_uri,
        bins,
        columns=columns,
        dtypes=dtypes,
        metadata=metadata,
        assembly=assembly,
        symmetric_upper=symmetric_upper,
        mode=mode,
        h5opts=h5opts,
        boundscheck=boundscheck,
        triucheck=triucheck,
        dupcheck=dupcheck,
        ensure_sorted=ensure_sorted,
        lock=lock,
        binsize=binsize,
    )
    with writer:
        logger.info("Writing pixels")
        for chunk in iterable:
            writer.write(chunk)


def create_from_unordered(
    cool_uri,
    bins,
    chunks,
    columns=None,
    dtypes=None,
    mode=None,
    mergebuf=20_000_000,
    delete_temp=True,
    temp_dir=None,
    max_merge=200,
    agg=None,
    duplicates="combine",
    **kwargs,
):
    """
    Create a Cooler in two passes via an external sort mechanism. In the first
    pass, a sequence of data chunks are processed and sorted in memory and saved
    to temporary Coolers. In the second pass, the temporary Coolers are merged
    into the output. This way the individual chunks do not need to be provided
    in any particular order.

    Pixels that occur more than once in the input are combined using ``agg``,
    a mapping of value column names to "sum" (default), "min" or "max". If
    ``duplicates`` is "error", a pixel occurring more than once anywhere in
    the input raises instead.

    """
    dtypes = _get_dtypes_arg(dtypes, kwargs)
    if columns is not None:
        columns = [col for col in columns if col not in {"bin1_id", "bin2_id"}]

    writer = _UnorderedCoolerWriter(
        cool_uri,
        bins,
        columns=columns,
        dtypes=dtypes,
        agg=agg,
        duplicates=duplicates,
        mergebuf=mergebuf,
        max_merge=max_merge,
        temp_dir=temp_dir,
        delete_temp=delete_temp,
        mode="w" if mode is None else mode,
        **kwargs,
    )
    with writer:
        for chunk in chunks:
            writer.write(chunk)


_DOC_OTHER_PARAMS = """
    columns : sequence of str, optional
        Customize which value columns from the input pixels to store in the
        cooler. Non-standard value columns will be given dtype ``float64``
        unless overriden using the ``dtypes`` argument. If ``None``, we only
        attempt to store a value column named ``"count"``.
    dtypes : dict, optional
        Dictionary mapping column names to dtypes. Can be used to override the
        default dtypes of ``bin1_id``, ``bin2_id`` or ``count`` or assign
        dtypes to custom value columns. Non-standard value columns given in
        ``dtypes`` must also be provided in the ``columns`` argument or they
        will be ignored.
    metadata : dict, optional
        Experiment metadata to store in the file. Must be JSON compatible.
    assembly : str, optional
        Name of genome assembly.
    ordered : bool, optional [default: False]
        If the input chunks of pixels are provided with correct triangularity
        and in ascending order of (``bin1_id``, ``bin2_id``), set this to
        ``True`` to write the cooler in one step.
        If ``False`` (default), we create the cooler in two steps using an
        external sort mechanism. See Notes for more details.
    symmetric_upper : bool, optional [default: True]
        If True, sets the file's storage-mode property to ``symmetric-upper``:
        use this only if the input data references the upper triangle of a
        symmetric matrix! For all other cases, set this option to False.
    mode : {'w' , 'a'}, optional [default: 'w']
        Write mode for the output file. 'a': if the output file exists, append
        the new cooler to it. 'w': if the output file exists, it will be
        truncated. Default is 'w'.

    Other parameters
    ----------------
    mergebuf : int, optional
        Maximum number of records to buffer in memory at any give time during
        the merge step.
    delete_temp : bool, optional
        Whether to delete temporary files when finished.
        Useful for debugging. Default is True.
    temp_dir : str, optional
        Create temporary files in a specified directory instead of the same
        directory as the output file. Pass ``-`` to use the system default.
    agg : dict, optional
        When ``ordered=False``, how to combine the values of a pixel that
        occurs more than once in the input: a mapping of value column names
        to "sum" (default), "min" or "max".
    duplicates : {"combine", "error"}, optional
        When ``ordered=False``, what a pixel occurring more than once in the
        input (within or across chunks) means: "combine" (default) combines
        its values using ``agg``; "error" raises. With ``ordered=True``,
        duplicates always raise unless ``dupcheck`` is False.
    max_merge : int, optional
        If merging more than ``max_merge`` chunks, do the merge recursively in
        two passes.
    h5opts : dict, optional
        HDF5 dataset filter options to use (compression, shuffling,
        checksumming, etc.). Default is to use autochunking and GZIP
        compression, level 6.
    lock : multiprocessing.Lock, optional
        Optional lock to control concurrent access to the output file.
    ensure_sorted : bool, optional
        Sort each input chunk. With ``ordered=True``, the chunks must still
        be in order relative to each other.
    boundscheck : bool, optional
        Input validation: Check that all bin IDs lie in the expected range.
    dupcheck : bool, optional
        Input validation: with ``ordered=True``, check that no pixel occurs
        more than once in the input. With ``ordered=False``, use
        ``duplicates`` instead.
    triucheck : bool, optional
        Input validation: Check that ``bin1_id`` <= ``bin2_id`` when creating
        coolers in symmetric-upper mode.
""".strip()

_DOC_NOTES = """
    Notes
    -----
    If the pixel chunks are provided in the correct order required for the
    output to be properly sorted, then the cooler can be created in a single
    step by setting ``ordered=True``.

    If not, the cooler is created in two steps via an external sort mechanism.
    In the first pass, the sequence of pixel chunks are processed and sorted in
    memory and saved to temporary coolers. In the second pass, the temporary
    coolers are merged into the output file. This way the individual chunks do
    not need to be provided in any particular order. When ``ordered=False``,
    the following options for the merge step are available: ``mergebuf``,
    ``delete_temp``, ``temp_dir``, ``max_merge``.

    Each chunk of pixels will go through a validation pipeline, which can be
    customized with the following options: ``boundscheck``, ``triucheck``,
    ``dupcheck``, ``ensure_sorted``.
""".strip()


def _format_docstring(**kwargs):
    def decorate(func):
        func.__doc__ = func.__doc__.format(**kwargs)
        return func

    return decorate


@_format_docstring(other_parameters=_DOC_OTHER_PARAMS, notes=_DOC_NOTES)
def create_cooler(
    cool_uri: str,
    bins: pd.DataFrame,
    pixels: Tabular | Iterable[Tabular],
    columns: list[str] | None = None,
    dtypes: dict[str, Any] | None = None,
    metadata: dict | None = None,
    assembly: str | None = None,
    ordered: bool = False,
    symmetric_upper: bool = True,
    mode: str = "w",
    mergebuf: int = 20_000_000,
    delete_temp: bool = True,
    temp_dir: str | None = None,
    max_merge: int = 200,
    agg: dict[str, str] | None = None,
    duplicates: str = "combine",
    boundscheck: bool = True,
    dupcheck: bool = True,
    triucheck: bool | None = None,
    ensure_sorted: bool = False,
    h5opts: dict | None = None,
    lock: Any | None = None,
) -> None:
    r"""
    Create a cooler from bins and pixels at the specified URI.

    Because the number of pixels is often very large, the input pixels are
    normally provided as an iterable (e.g., an iterator or generator) of
    DataFrame **chunks** that fit in memory.

    .. versionadded:: 0.8.0

    Parameters
    ----------
    cool_uri : str
        Path to cooler file or URI string. If the file does not exist,
        it will be created.
    bins : pandas.DataFrame
        Segmentation of the chromosomes into genomic bins as a BED-like
        DataFrame with columns ``chrom``, ``start`` and ``end``. May contain
        additional columns.
    pixels : DataFrame, dictionary, or iterable of either
        A table, given as a dataframe or a column-oriented dict, containing
        columns labeled ``bin1_id``, ``bin2_id`` and ``count``, sorted by
        (``bin1_id``, ``bin2_id``). If additional columns are included in the
        pixel table, their names and dtypes must be specified using the
        ``columns`` and ``dtypes`` arguments. For larger input data, an
        **iterable** can be provided that yields the pixel data as a sequence
        of chunks. If the input is a dask DataFrame, it will also be processed
        one chunk at a time.
    {other_parameters}

    See also
    --------
    cooler.create_scool
    cooler.create.sanitize_records
    cooler.create.sanitize_pixels

    {notes}

    """
    # dispatch to the approprate creation method
    if isinstance(pixels, (pd.DataFrame, dict)):
        pixels = pd.DataFrame(pixels).sort_values(["bin1_id", "bin2_id"])
        if not ordered:
            # Repeated pixels are combined, as for an iterable of chunks
            pixels = _combine_duplicates(pixels, columns, agg, duplicates)
        ordered = True

    if ordered:
        create(
            cool_uri,
            bins,
            pixels,
            columns=columns,
            dtypes=dtypes,
            metadata=metadata,
            assembly=assembly,
            symmetric_upper=symmetric_upper,
            mode=mode,
            boundscheck=boundscheck,
            dupcheck=dupcheck,
            triucheck=triucheck,
            ensure_sorted=ensure_sorted,
            h5opts=h5opts,
            lock=lock,
        )
    else:
        create_from_unordered(
            cool_uri,
            bins,
            pixels,
            columns=columns,
            dtypes=dtypes,
            metadata=metadata,
            assembly=assembly,
            symmetric_upper=symmetric_upper,
            mode=mode,
            boundscheck=boundscheck,
            dupcheck=dupcheck,
            triucheck=triucheck,
            ensure_sorted=ensure_sorted,
            h5opts=h5opts,
            lock=lock,
            mergebuf=mergebuf,
            delete_temp=delete_temp,
            temp_dir=temp_dir,
            max_merge=max_merge,
            agg=agg,
            duplicates=duplicates,
        )


def _combine_duplicates(
    pixels: pd.DataFrame,
    columns: list[str] | None,
    agg: dict[str, str] | None,
    duplicates: str,
) -> pd.DataFrame:
    """Combine the values of repeated pixels in sorted in-memory pixels."""
    if duplicates not in ("combine", "error"):
        raise ValueError(
            f"duplicates must be 'combine' or 'error', not '{duplicates}'."
        )
    keys = ["bin1_id", "bin2_id"]
    if columns is None:
        columns = ["count"]
    value_cols = [c for c in columns if c not in keys and c in pixels.columns]
    aggs = _value_aggs(value_cols, agg)
    repeated = pixels.duplicated(keys)
    if not repeated.any():
        return pixels
    if duplicates == "error":
        i = int(np.argmax(repeated.to_numpy()))
        pixel = tuple(int(x) for x in pixels[keys].iloc[i])
        raise BadInputError(f"Found a duplicate pixel {pixel} in the input.")
    return pixels.groupby(keys, sort=False).agg(aggs).reset_index()


@_format_docstring(other_parameters=_DOC_OTHER_PARAMS, notes=_DOC_NOTES)
def create_scool(
    cool_uri: str,
    bins: pd.DataFrame,
    cell_name_pixels_dict: dict[str, pd.DataFrame],
    columns: list[str] | None = None,
    dtypes: dict[str, Any] | None = None,
    metadata: dict | None = None,
    assembly: str | None = None,
    ordered: bool = False,
    symmetric_upper: bool = True,
    mode: str = "w",
    mergebuf: int = 20_000_000,
    delete_temp: bool = True,
    temp_dir: str | None = None,
    max_merge: int = 200,
    boundscheck: bool = True,
    dupcheck: bool = True,
    triucheck: bool | None = None,
    ensure_sorted: bool = False,
    h5opts: dict | None = None,
    lock: Any | None = None,
    **kwargs,
) -> None:
    r"""
    Create a single-cell (scool) file.

    For each cell store a cooler matrix under **/cells**, where all matrices
    have the same dimensions.

    Each cell is a regular cooler data collection, so the input must be a
    bin table and pixel table for each cell. The pixel tables are provided as
    a dictionary where the key is a unique cell name. The bin tables can be
    provided as a dict with the same keys or a single common bin table can be
    given.

    .. versionadded:: 0.8.9

    Parameters
    ----------
    cool_uri : str
        Path to scool file or URI string. If the file does not exist,
        it will be created.
    bins : :class:`pandas.DataFrame` or Dict[str, DataFrame]
        A single bin table or dictionary of cell names to bins tables. A bin
        table is a dataframe with columns ``chrom``, ``start`` and ``end``.
        May contain additional columns.
    cell_name_pixels_dict : Dict[str, DataFrame]
        Cell name as key and pixel table DataFrame as value.
        A table, given as a dataframe or a column-oriented dict, containing
        columns labeled ``bin1_id``, ``bin2_id`` and ``count``, sorted by
        (``bin1_id``, ``bin2_id``). If additional columns are included in the
        pixel table, their names and dtypes must be specified using the
        ``columns`` and ``dtypes`` arguments. For larger input data, an
        **iterable** can be provided that yields the pixel data as a sequence
        of chunks. If the input is a dask DataFrame, it will also be processed
        one chunk at a time.
    {other_parameters}

    See also
    --------
    cooler.create_cooler
    cooler.zoomify_cooler

    {notes}

    """

    if isinstance(bins, pd.DataFrame):
        bins_dict = None
        shared_bins = bins
    else:
        # A dict of cell name -> bin table, with matching keys
        bins_dict = bins
        if len(bins_dict) == 0:
            raise ValueError("At least one bin must be given.")
        if set(bins_dict) != set(cell_name_pixels_dict):
            raise ValueError("Bins and pixel dicts do not have matching keys")
        shared_bins = bins_dict[next(iter(bins_dict))][["chrom", "start", "end"]]

    dtypes = _get_dtypes_arg(dtypes, kwargs)
    if kwargs:
        raise TypeError(
            "create_scool() got unexpected keyword arguments: "
            f"{', '.join(sorted(kwargs))}"
        )

    # Cells are stored under the last component of their name
    cell_names = sorted(cell_name_pixels_dict)
    stored_names = pd.Series([key.split("/")[-1] for key in cell_names], dtype=object)
    if stored_names.duplicated().any():
        dups = sorted(set(stored_names[stored_names.duplicated()]))
        raise ValueError(f"Cell names would be stored under the same name: {dups}")

    writer = _ScoolWriter(
        cool_uri,
        shared_bins,
        columns=columns,
        dtypes=dtypes,
        metadata=metadata,
        assembly=assembly,
        symmetric_upper=symmetric_upper,
        mode=mode,
        h5opts=h5opts,
        boundscheck=boundscheck,
        triucheck=triucheck,
        dupcheck=dupcheck,
        ensure_sorted=ensure_sorted,
        lock=lock,
    )
    with writer:
        for key, name in zip(cell_names, stored_names):
            writer.add_cell(
                name,
                cell_name_pixels_dict[key],
                bins=None if bins_dict is None else bins_dict[key],
                ordered=ordered,
                mergebuf=mergebuf,
                max_merge=max_merge,
                temp_dir=temp_dir,
                delete_temp=delete_temp,
            )
