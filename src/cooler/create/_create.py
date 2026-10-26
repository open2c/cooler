from __future__ import annotations

import json
import warnings
from collections.abc import Iterable
from datetime import datetime
from typing import Any

import h5py
import numpy as np
import pandas as pd

from .._logging import get_logger
from .._typing import Tabular
from .._version import __format_version__, __format_version_scool__, __version__
from ..core import get, put
from ..util import (
    infer_meta,
    parse_cooler_uri,
    rlencode,
)
from . import (
    BIN1OFFSET_DTYPE,
    BIN_DTYPE,
    CHROM_DTYPE,
    CHROMID_DTYPE,
    CHROMOFFSET_DTYPE,
    CHROMSIZE_DTYPE,
    COORD_DTYPE,
    COUNT_DTYPE,
    MAGIC,
    MAGIC_SCOOL,
    PIXEL_FIELDS,
    URL,
)

logger = get_logger("cooler.create")


def write_chroms(grp: h5py.Group, chroms: pd.DataFrame, h5opts: dict) -> None:
    """
    Write the chromosome table.

    Parameters
    ----------
    grp : h5py.Group
        Group handle of an open HDF5 file with write permissions.
    chroms : DataFrame
        Chromosome table containing at least 'chrom' and 'length' columns
    h5opts : dict
        HDF5 dataset filter options.

    """
    n_chroms = len(chroms)
    names = np.array(chroms["name"], dtype=CHROM_DTYPE)  # auto-adjusts char length
    grp.create_dataset(
        "name", shape=(n_chroms,), dtype=names.dtype, data=names, **h5opts
    )
    grp.create_dataset(
        "length",
        shape=(n_chroms,),
        dtype=CHROMSIZE_DTYPE,
        data=chroms["length"],
        **h5opts,
    )

    # Extra columns
    columns = list(chroms.keys())
    for col in ["name", "length"]:
        columns.remove(col)
    if columns:
        put(grp, chroms[columns])


def write_bins(
    grp: h5py.Group,
    bins: pd.DataFrame,
    chromnames: list[str],
    h5opts: dict,
    chrom_as_enum: bool = True,
) -> None:
    """
    Write the genomic bin table.

    Parameters
    ----------
    grp : h5py.Group
        Group handle of an open HDF5 file with write permissions.
    bins : pandas.DataFrame
        BED-like data frame with at least three columns: ``chrom``, ``start``,
        ``end``, sorted by ``chrom`` then ``start``, and forming a complete
        genome segmentation. The ``chrom`` column must be sorted according to
        the ordering in ``chroms``.
    chromnames : sequence of str
        Contig names.
    h5opts : dict
        HDF5 dataset filter options.

    """
    n_chroms = len(chromnames)
    n_bins = len(bins)
    idmap = dict(zip(chromnames, range(n_chroms)))

    # Convert chrom names to enum
    chrom_ids = [idmap[chrom] for chrom in bins["chrom"]]
    if chrom_as_enum:
        chrom_dtype = h5py.special_dtype(enum=(CHROMID_DTYPE, idmap))
    else:
        chrom_dtype = CHROMID_DTYPE

    # Store bins
    try:
        chrom_dset = grp.create_dataset(
            "chrom", shape=(n_bins,), dtype=chrom_dtype, data=chrom_ids, **h5opts
        )
    except ValueError:
        # If too many scaffolds for HDF5 enum header,
        # try storing chrom IDs as raw int instead
        if chrom_as_enum:
            chrom_as_enum = False
            chrom_dtype = CHROMID_DTYPE
            chrom_dset = grp.create_dataset(
                "chrom", shape=(n_bins,), dtype=chrom_dtype, data=chrom_ids, **h5opts
            )
        else:
            raise
    if not chrom_as_enum:
        chrom_dset.attrs["enum_path"] = "/chroms/name"

    grp.create_dataset(
        "start", shape=(n_bins,), dtype=COORD_DTYPE, data=bins["start"], **h5opts
    )
    grp.create_dataset(
        "end", shape=(n_bins,), dtype=COORD_DTYPE, data=bins["end"], **h5opts
    )

    # Extra columns
    columns = list(bins.keys())
    for col in ["chrom", "start", "end"]:
        columns.remove(col)
    if columns:
        put(grp, bins[columns])


def prepare_pixels(
    grp: h5py.Group,
    n_bins: int,
    max_size: int,
    columns: pd.Index | list[str],
    dtypes: pd.Series | dict[str, Any],
    h5opts: dict,
) -> None:
    columns = list(columns)
    init_size = min(5 * n_bins, max_size)
    grp.create_dataset(
        "bin1_id",
        dtype=dtypes.get("bin1_id", BIN_DTYPE),
        shape=(init_size,),
        maxshape=(max_size,),
        **h5opts,
    )
    grp.create_dataset(
        "bin2_id",
        dtype=dtypes.get("bin2_id", BIN_DTYPE),
        shape=(init_size,),
        maxshape=(max_size,),
        **h5opts,
    )

    if "count" in columns:
        grp.create_dataset(
            "count",
            dtype=dtypes.get("count", COUNT_DTYPE),
            shape=(init_size,),
            maxshape=(max_size,),
            **h5opts,
        )

    for col in ["bin1_id", "bin2_id", "count"]:
        try:
            columns.remove(col)
        except ValueError:
            pass

    if columns:
        for col in columns:
            grp.create_dataset(
                col,
                dtype=dtypes.get(col, float),
                shape=(init_size,),
                maxshape=(max_size,),
                **h5opts,
            )


def _check_fits_dtype(col: str, values: np.ndarray, dtype: np.dtype) -> None:
    """Raise if integer values would be clipped when stored as ``dtype``."""
    dtype = np.dtype(dtype)
    values = np.asarray(values)
    if dtype.kind not in "iu" or values.dtype.kind not in "iuf" or not len(values):
        return
    info = np.iinfo(dtype)
    vmin, vmax = values.min(), values.max()
    if vmin < info.min or vmax > info.max:
        raise ValueError(
            f"Values in pixel column '{col}' range from {vmin} to {vmax}, which "
            f"doesn't fit in the column's data type ({dtype}). Retry storing the "
            f"column as int64, e.g. with dtypes={{'{col}': 'int64'}} in Python "
            f"or --field {col}:dtype=int64 on the command line."
        )


def write_pixels(
    filepath: str,
    grouppath: str,
    columns: pd.Index | list[str],
    iterable: Tabular | Iterable[Tabular],
    h5opts: dict,
    lock: Any | None,
) -> tuple[int, int]:
    """
    Write the non-zero pixel table.

    Parameters
    ----------
    filepath : str
        Path to HDF5 output file.
    grouppath : str
        Qualified path to destination HDF5 group.
    columns : sequence
        Sequence of column names
    iterable : an iterable object
        An object that processes and yields binned contacts from some input
        source as a stream of chunks. The chunks must be either pandas
        DataFrames or mappings of column names to arrays.
    h5opts : dict
        HDF5 filter options.
    lock : multiprocessing.Lock, optional
        Optional lock to synchronize concurrent HDF5 file access.

    """
    nnz = 0
    total = 0
    for i, chunk in enumerate(iterable):
        if isinstance(chunk, pd.DataFrame):
            chunk = {k: v.values for k, v in chunk.items()}

        try:
            if lock is not None:
                lock.acquire()

            logger.debug(f"writing chunk {i}")

            with h5py.File(filepath, "r+") as fw:
                grp = fw[grouppath]
                dsets = [grp[col] for col in columns]

                n = len(chunk[columns[0]])
                for col, dset in zip(columns, dsets):
                    _check_fits_dtype(col, chunk[col], dset.dtype)
                for col, dset in zip(columns, dsets):
                    dset.resize((nnz + n,))
                    dset[nnz : nnz + n] = chunk[col]
                nnz += n
                if "count" in chunk:
                    total += chunk["count"].sum()

                fw.flush()

        finally:
            if lock is not None:
                lock.release()

    return nnz, total


def index_pixels(grp: h5py.Group, n_bins: int, nnz: int) -> np.ndarray:
    bin1 = grp["bin1_id"]
    bin1_offset = np.zeros(n_bins + 1, dtype=BIN1OFFSET_DTYPE)
    curr_val = 0
    for start, _length, value in zip(*rlencode(bin1, 1000000)):
        bin1_offset[curr_val : value + 1] = start
        curr_val = value + 1
    bin1_offset[curr_val:] = nnz
    return bin1_offset


def index_bins(grp: h5py.Group, n_chroms: int, n_bins: int) -> np.ndarray:
    chrom_ids = grp["chrom"]
    chrom_offset = np.zeros(n_chroms + 1, dtype=CHROMOFFSET_DTYPE)
    curr_val = 0
    for start, _length, value in zip(*rlencode(chrom_ids)):
        chrom_offset[curr_val : value + 1] = start
        curr_val = value + 1
    chrom_offset[curr_val:] = n_bins
    return chrom_offset


def write_indexes(
    grp: h5py.Group,
    chrom_offset: np.ndarray,
    bin1_offset: np.ndarray,
    h5opts: dict
) -> None:
    """
    Write the indexes.

    Parameters
    ----------
    grp : h5py.Group
        Group handle of an open HDF5 file with write permissions.
    chrom_offset : sequence
        Lookup table: chromosome ID -> first row in bin table (bin ID)
        corresponding to that chromosome.
    bin1_offset : sequence
        Lookup table: genomic bin ID -> first row in pixel table (pixel ID)
        having that bin on the first axis.

    """
    grp.create_dataset(
        "chrom_offset",
        shape=(len(chrom_offset),),
        dtype=CHROMOFFSET_DTYPE,
        data=chrom_offset,
        **h5opts,
    )
    grp.create_dataset(
        "bin1_offset",
        shape=(len(bin1_offset),),
        dtype=BIN1OFFSET_DTYPE,
        data=bin1_offset,
        **h5opts,
    )


def write_info(grp: h5py.Group, info: dict, scool: bool = False) -> None:
    """
    Write the file description and metadata attributes.

    Parameters
    ----------
    grp : h5py.Group
        Group handle of an open HDF5 file with write permissions.
    info : dict
        Dictionary, unnested with the possible exception of the ``metadata``
        key. ``metadata``, if present, must be JSON-serializable.

    Required keys
    -------------
    nbins : int
        number of genomic bins
    nnz : int
        number of non-zero pixels

    """
    assert "nbins" in info
    if not scool:
        assert "nnz" in info
    info.setdefault("genome-assembly", "unknown")
    info["metadata"] = json.dumps(info.get("metadata", {}))
    info["creation-date"] = datetime.now().isoformat()
    info["generated-by"] = "cooler-" + __version__
    if scool:
        info["format"] = MAGIC_SCOOL
        info["format-version"] = __format_version_scool__
    else:
        info["format"] = MAGIC
        info["format-version"] = __format_version__
    info["format-url"] = URL
    grp.attrs.update(info)


def _rename_chroms(
    grp: h5py.Group, rename_dict: dict[str, str], h5opts: dict
) -> None:
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


def rename_chroms(
    clr,
    rename_dict: dict[str, str],
    h5opts: dict | None = None
) -> None:
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
    dtypes: pd.Series | dict | None,
    kwargs: dict[str, Any]
) -> pd.Series | dict | None:
    if "dtype" in kwargs:
        if dtypes is None:
            dtypes = kwargs.pop("dtype")
            warnings.warn(
                "Use dtypes= instead of dtype=", FutureWarning, stacklevel=2
            )
        else:
            raise ValueError(
                'Received both "dtypes" and "dtype" arguments. '
                'Please use "dtypes" to provide a column name -> dtype mapping. '
                '"dtype" remains as an alias but is deprecated.'
            )
    return dtypes


# Aggregations that can be applied to per-chunk results to combine them, and
# the aggregation to combine them with.
_MERGEABLE_AGGS = {"sum": "sum", "min": "min", "max": "max", "count": "sum"}


def _set_h5opts(h5opts: dict[str, Any]) -> dict[str, Any]:
    result = {}
    if h5opts is not None:
        result.update(h5opts)
    available_opts = {
        "chunks",
        "maxshape",
        "compression",
        "compression_opts",
        "scaleoffset",
        "shuffle",
        "fletcher32",
        "fillvalue",
        "track_times",
    }
    for key in result:
        if key not in available_opts:
            raise ValueError(f"Unknown storage option '{key}'.")
    result.setdefault("compression", "gzip")
    if result["compression"] == "gzip" and "compression_opts" not in result:
        result["compression_opts"] = 6
    result.setdefault("shuffle", True)
    return result


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
    triucheck=True,
    dupcheck=True,
    ensure_sorted=False,
    lock=None,
    append=False,
    **kwargs,
):
    """
    Create a new Cooler.

    Deprecated parameters
    ---------------------
    chromsizes : Series
        Chromsizes are now inferred from ``bins``.
    append : bool, optional
        Append new Cooler to the file if it exists. If False, an existing file
        with the same name will be truncated. Default is False.
        Use the ``mode`` argument instead.
    dtype : dict, optional
        Dictionary mapping column names in the pixel table to dtypes.
        Use the ``dtypes`` argument instead.

    """
    from ._writer import _CoolerWriter, _pixel_schema

    if mode is None:
        mode = "a" if append else "w"

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

    # Determine the appropriate iterable
    try:
        from dask.dataframe import DataFrame as dask_df
    except (ImportError, AttributeError):  # pragma: no cover
        dask_df = ()

    if isinstance(pixels, dask_df):
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
    **kwargs,
):
    """
    Create a Cooler in two passes via an external sort mechanism. In the first
    pass, a sequence of data chunks are processed and sorted in memory and saved
    to temporary Coolers. In the second pass, the temporary Coolers are merged
    into the output. This way the individual chunks do not need to be provided
    in any particular order.

    Pixels that occur in more than one chunk are combined using ``agg``, a
    mapping of value column names to aggregation names (default: "sum"). When
    there is more than one chunk, only "sum", "min", "max", "count" and "mean"
    are supported. For "mean", each occurrence of a pixel counts as one value,
    unless chunks provide the sum and number of values of the column (see
    ``aggregate_records(..., split_means=True)``).

    """
    from ._writer import _UnorderedCoolerWriter

    dtypes = _get_dtypes_arg(dtypes, kwargs)
    if columns is not None:
        columns = [col for col in columns if col not in {"bin1_id", "bin2_id"}]

    writer = _UnorderedCoolerWriter(
        cool_uri,
        bins,
        columns=columns,
        dtypes=dtypes,
        agg=agg,
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


def append(
    cool_uri, table, data, chunked=False, force=False, h5opts=None, lock=None
):  # pragma: no cover
    """
    Append one or more data columns to an existing table.

    Parameters
    ----------
    cool_uri : str
        Path to Cooler file or URI to Cooler group.
    table : str
        Name of table (HDF5 group).
    data : dict-like
        DataFrame, Series or mapping of column names to data. If the input is a
        dask DataFrame or Series, the data is written in chunks.
    chunked : bool, optional
        If True, the values of the data dict are treated as separate chunk
        iterators of column data.
    force : bool, optional
        If True, replace existing columns with the same name as the input.
    h5opts : dict, optional
        HDF5 dataset filter options to use (compression, shuffling,
        checksumming, etc.). Default is to use autochunking and GZIP
        compression, level 6.
    lock : multiprocessing.Lock, optional
        Optional lock to synchronize concurrent HDF5 file access.

    """
    h5opts = _set_h5opts(h5opts)

    file_path, group_path = parse_cooler_uri(cool_uri)

    try:
        from dask.dataframe import DataFrame as dask_df
        from dask.dataframe import Series as dask_series
    except (ImportError, AttributeError):
        dask_df = ()
        dask_series = ()

    if isinstance(data, dask_series):
        data = data.to_frame()

    try:
        names = data.keys()
    except AttributeError:
        names = data.columns

    with h5py.File(file_path, "r+") as f:
        h5 = f[group_path]
        for name in names:
            if name in h5[table]:
                if not force:
                    raise ValueError(
                        f"'{name}' column already exists. "
                        + "Use --force option to overwrite."
                    )
                else:
                    del h5[table][name]

        if isinstance(data, dask_df):
            # iterate over dataframe chunks
            for chunk in data.to_delayed():
                i = 0
                for chunk in data.to_delayed():
                    chunk = chunk.compute()
                    try:
                        if lock is not None:
                            lock.acquire()
                        put(h5[table], chunk, lo=i, h5opts=h5opts)
                    finally:
                        if lock is not None:
                            lock.release()
                    i += len(chunk)
        elif chunked:
            # iterate over chunks from each column
            for key in data:
                i = 0
                for chunk in data[key]:
                    try:
                        if lock is not None:
                            lock.acquire()
                        put(h5[table], {key: chunk}, lo=i, h5opts=h5opts)
                    finally:
                        if lock is not None:
                            lock.release()
                    i += len(chunk)
        else:
            # write all the data
            try:
                if lock is not None:
                    lock.acquire()
                put(h5[table], data, lo=0, h5opts=h5opts)
            finally:
                if lock is not None:
                    lock.release()


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
        occurs in more than one input chunk: a mapping of value column names
        to "sum" (default), "min", "max", "count" or "mean".
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
        Ensure that each input chunk is properly sorted.
    boundscheck : bool, optional
        Input validation: Check that all bin IDs lie in the expected range.
    dupcheck : bool, optional
        Input validation: Check that no duplicate pixels exist within any chunk.
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
    boundscheck: bool = True,
    dupcheck: bool = True,
    triucheck: bool = True,
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
        )


def _same_chroms_and_bins(
    grp: h5py.Group, chroms: pd.DataFrame, bins: pd.DataFrame
) -> bool:
    if "chroms" not in grp or "bins" not in grp:
        return False
    old_chroms = get(grp["chroms"], 0, None, ["name", "length"], convert_enum=False)
    old_bins = get(grp["bins"], 0, None, ["chrom", "start", "end"], convert_enum=True)
    new_bins = bins[["chrom", "start", "end"]]
    return (
        old_chroms["name"].astype(str).tolist() == chroms["name"].astype(str).tolist()
        and np.array_equal(old_chroms["length"], chroms["length"])
        and len(old_bins) == len(new_bins)
        and np.array_equal(old_bins["chrom"].astype(str), new_bins["chrom"].astype(str))
        and np.array_equal(old_bins["start"], new_bins["start"])
        and np.array_equal(old_bins["end"], new_bins["end"])
    )


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
    triucheck: bool = True,
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

    from ._writer import _ScoolWriter

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
    stored_names = pd.Series([key.split("/")[-1] for key in cell_names])
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
