from __future__ import annotations

import bisect
import inspect
import json
import os
import os.path as op
import posixpath
import tempfile
import warnings
from datetime import datetime
from typing import Any

import h5py
import numpy as np
import pandas as pd

from .._logging import get_logger
from .._version import __format_version__, __format_version_scool__, __version__
from ..core import get, put
from ..util import get_binsize, get_chromsizes, get_meta, parse_cooler_uri, rlencode
from ._constants import (
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
    PIXEL_DTYPES,
    URL,
)
from ._ingest import BadInputError, _check_pixel_ids

logger = get_logger("cooler.create")

# Suffix of the temporary group a cooler is written into before it is moved
# into place.
_WRITING_SUFFIX = ".__writing__"

# Standard tables of a cooler group.
_TABLES = ("chroms", "bins", "pixels", "indexes")


# Aggregations that can be applied to per-chunk results to combine them, and
# the aggregation to combine them with.
_MERGEABLE_AGGS = {"sum": "sum", "min": "min", "max": "max", "count": "sum"}

def _set_h5opts(h5opts: dict[str, Any]) -> dict[str, Any]:
    result = {}
    if h5opts is not None:
        result.update(h5opts)
    available_opts = {
        "chunks",
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

# Storage options that only apply to numeric datasets
_NUMERIC_H5OPTS = {"scaleoffset", "fillvalue"}


def _h5opts_for_dtype(h5opts: dict, dtype: Any) -> dict:
    """Drop storage options that don't apply to datasets of this dtype."""
    if np.dtype(dtype).kind in "biuf":
        return h5opts
    return {k: v for k, v in h5opts.items() if k not in _NUMERIC_H5OPTS}


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
        "name",
        shape=(n_chroms,),
        dtype=names.dtype,
        data=names,
        **_h5opts_for_dtype(h5opts, names.dtype),
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


def _check_bins(bins: pd.DataFrame) -> None:
    """
    Check that the bins of each chromosome are consecutive rows that tile the
    chromosome from 0 (#254). Chromosome lengths are taken from the last bins,
    so they can't be checked.
    """
    codes, names = pd.factorize(bins["chrom"])
    starts = bins["start"].to_numpy()
    ends = bins["end"].to_numpy()
    first = np.r_[True, codes[1:] != codes[:-1]] if len(codes) else codes == 0

    run_codes = codes[first]
    if len(np.unique(run_codes)) < len(run_codes):
        _, counts = np.unique(run_codes, return_counts=True)
        chrom = names[np.unique(run_codes)[counts > 1][0]]
        raise ValueError(
            f"The bins of '{chrom}' are not consecutive rows of the bin table."
        )
    problems = [
        (ends <= starts, "has a bin whose end is not after its start"),
        (first & (starts != 0), "does not start at 0"),
        (~first & (starts != np.r_[0, ends[:-1]]), "has gaps or overlapping bins"),
    ]
    for mask, message in problems:
        if mask.any():
            chrom = names[codes[np.argmax(mask)]]
            raise ValueError(f"The bin table for '{chrom}' {message}.")


def _check_binsize(bins: pd.DataFrame, binsize: int) -> None:
    sizes = (bins["end"] - bins["start"]).to_numpy()
    codes = pd.factorize(bins["chrom"])[0]
    last = np.r_[codes[1:] != codes[:-1], True]
    if (sizes[~last] != binsize).any() or (sizes[last] > binsize).any():
        raise ValueError(f"The bins don't have a uniform size of {binsize}.")


def _prepare_tables(
    bins: pd.DataFrame, binsize: int | None = None
) -> tuple[pd.DataFrame, pd.DataFrame, Any]:
    """
    Return the bin table to write, the chrom table, and the bin size.

    A known bin size can be given, since it can't always be inferred from the
    table, e.g. when every chromosome is a single bin.
    """
    for col in ["chrom", "start", "end"]:
        if col not in bins.columns:
            raise ValueError(f"Missing column from bin table: '{col}'.")
    _check_bins(bins)
    if binsize is None:
        binsize = get_binsize(bins)
    else:
        binsize = int(binsize)
        _check_binsize(bins, binsize)
    bins = bins.copy()
    bins["chrom"] = bins["chrom"].astype(object)
    chromsizes = get_chromsizes(bins)
    try:
        chromsizes = chromsizes.items()
    except AttributeError:
        pass
    names, lengths = zip(*chromsizes)
    chroms = pd.DataFrame(
        {"name": names, "length": lengths}, columns=["name", "length"]
    )
    return bins, chroms, binsize


def _pixel_schema(
    columns: list[str] | None, dtypes: dict | None
) -> tuple[list[str], dict[str, Any]]:
    """Return the pixel column names and dtypes to store."""
    if columns is None:
        columns = ["bin1_id", "bin2_id", "count"]
    else:
        columns = list(columns)
        for col in ["bin1_id", "bin2_id"]:  # don't include count!
            if col not in columns:
                columns.insert(0, col)
    all_dtypes = dict(PIXEL_DTYPES)
    if dtypes is not None:
        all_dtypes.update(dict(dtypes))
    meta = get_meta(columns, all_dtypes, default_dtype=float)
    return list(meta.columns), dict(meta.dtypes)


class _CoolerWriter:
    """
    Write a cooler from a stream of pixel chunks sorted by (bin1_id, bin2_id).

    The chrom and bin tables are written when the writer is created. Each call
    to :meth:`write` validates a chunk of pixels and appends it to the pixel
    table, and :meth:`finish` writes the indexes and metadata.

    The cooler is written into a temporary group and only moved to its
    destination when it is finished, so a failure never leaves a partial cooler
    at the destination. Used as a context manager, the temporary group is
    removed if an exception is raised.

    Parameters
    ----------
    cool_uri : str
        Destination file path, optionally followed by ``::`` and a group path.
    bins : pandas.DataFrame
        Bin table with at least the columns ``chrom``, ``start`` and ``end``.
    columns : list of str, optional
        Value columns to store. Default is ``["count"]``.
    dtypes : dict, optional
        Dtypes of the stored pixel columns, overriding the defaults.
    metadata : dict, optional
        JSON-serializable metadata to store.
    assembly : str, optional
        Name of the genome assembly.
    symmetric_upper : bool, optional
        Whether pixels are the upper triangle of a symmetric matrix.
    mode : {"w", "a"}, optional
        "w" truncates the destination file; "a" adds to or replaces a group in
        an existing file.
    h5opts : dict, optional
        HDF5 dataset storage options.
    boundscheck : bool, optional
        Check that bin IDs are within the bin table.
    triucheck : bool, optional
        Check that pixels are in the upper triangle (``bin1_id <= bin2_id``).
        Default: only for ``symmetric_upper`` storage.
    dupcheck : bool, optional
        Check that no pixel occurs more than once. Pixels must always be in
        sorted order, across all chunks; if ``dupcheck`` is False, repeated
        pixels are allowed.
    ensure_sorted : bool, optional
        Sort each chunk before writing it. Chunks must still be in order
        relative to each other.
    lock : multiprocessing.Lock, optional
        If given, the file is only opened while holding the lock, for each
        write. Use this when other processes read from the same file while it
        is being written. Otherwise, the file is kept open until the writer is
        finished.
    binsize : int, optional
        Bin size of the bin table, if known. By default, it is inferred from
        the bin table, which isn't possible when every chromosome is a single
        bin.

    """

    def __init__(
        self,
        cool_uri: str,
        bins: pd.DataFrame,
        *,
        columns: list[str] | None = None,
        dtypes: dict | None = None,
        metadata: dict | None = None,
        assembly: str | None = None,
        symmetric_upper: bool = True,
        mode: str = "w",
        h5opts: dict | None = None,
        boundscheck: bool = True,
        triucheck: bool | None = None,
        dupcheck: bool = True,
        ensure_sorted: bool = False,
        lock: Any | None = None,
        link_tables_from: str | None = None,
        binsize: int | None = None,
    ):
        file_path, group_path = parse_cooler_uri(cool_uri)
        self.file_path = op.realpath(file_path)
        self.group_path = group_path
        # In mode "w", write a new file next to the destination and replace the
        # destination with it when finished, so a failure leaves any existing
        # file untouched.
        self._replace_path = None
        if mode == "w":
            self._replace_path = self.file_path
            self.file_path = self.file_path + _WRITING_SUFFIX
        if group_path == "/":
            self._tmp_path = "/" + _WRITING_SUFFIX
        else:
            self._tmp_path = group_path.rstrip("/") + _WRITING_SUFFIX

        if not isinstance(bins, pd.DataFrame):
            raise ValueError(  # noqa: TRY004
                "The bin table must be a pandas DataFrame."
            )
        self.bins, self.chroms, self.binsize = _prepare_tables(bins, binsize)
        self.n_bins = len(self.bins)
        self.columns, self.dtypes = _pixel_schema(columns, dtypes)
        self.h5opts = _set_h5opts(h5opts)
        self.metadata = metadata
        self.assembly = assembly
        self.symmetric_upper = symmetric_upper

        if triucheck is None:
            triucheck = symmetric_upper
        elif triucheck and not symmetric_upper:
            warnings.warn(
                "Creating a non-symmetric matrix, but `triucheck` was set to "
                "True. Changing to False.",
                stacklevel=2,
            )
            triucheck = False
        self.boundscheck = boundscheck
        self.triucheck = triucheck
        self.dupcheck = dupcheck
        self.ensure_sorted = ensure_sorted
        self.lock = lock

        self.nnz = 0
        self.sum = 0
        self._last_key = -1
        self._file = None
        self._done = False

        with self._open(mode) as f:
            logger.info(f'Creating cooler at "{self.file_path}::{group_path}"')
            if self._tmp_path in f:  # left over from a failed run
                del f[self._tmp_path]
            grp = f.create_group(self._tmp_path)
            if link_tables_from is None:
                write_chroms(grp.create_group("chroms"), self.chroms, self.h5opts)
                write_bins(
                    grp.create_group("bins"), self.bins, self.chroms["name"],
                    self.h5opts,
                )
            else:
                self._link_tables(f, grp, link_tables_from)
            if symmetric_upper:
                max_size = self.n_bins * (self.n_bins - 1) // 2 + self.n_bins
            else:
                max_size = self.n_bins * self.n_bins
            prepare_pixels(
                grp.create_group("pixels"), self.n_bins, max_size, self.columns,
                self.dtypes, self.h5opts,
            )

        if lock is None:
            self._file = h5py.File(self.file_path, "r+")

    def _link_tables(self, f: h5py.File, grp: h5py.Group, source_uri: str) -> None:
        # Hard link the chroms table and the chrom, start and end columns of
        # the bins table from another cooler group in the same file (used for
        # scool cells). Any other bin columns are stored separately.
        _, src_group = parse_cooler_uri(source_uri)
        src = f[src_group]
        grp["chroms"] = src["chroms"]
        bins_grp = grp.create_group("bins")
        for col in ["chrom", "start", "end"]:
            bins_grp[col] = src["bins"][col]
        extra = [c for c in self.bins.columns if c not in ("chrom", "start", "end")]
        if extra:
            put(bins_grp, self.bins[extra])

    class _Handle:
        """Context manager yielding the open file, under the lock if any."""

        def __init__(self, writer: _CoolerWriter, mode: str):
            self.writer = writer
            self.mode = mode
            self.file = None

        def __enter__(self) -> h5py.File:
            w = self.writer
            if w._file is not None:
                return w._file
            if w.lock is not None:
                w.lock.acquire()
            try:
                self.file = h5py.File(w.file_path, self.mode)
            except BaseException:
                if w.lock is not None:
                    w.lock.release()
                raise
            return self.file

        def __exit__(self, *exc) -> None:
            if self.file is not None:
                try:
                    self.file.flush()
                    self.file.close()
                finally:
                    if self.writer.lock is not None:
                        self.writer.lock.release()

    def _open(self, mode: str = "r+") -> _CoolerWriter._Handle:
        return self._Handle(self, mode)

    def _as_columns(self, chunk: Any) -> dict[str, np.ndarray]:
        if isinstance(chunk, pd.DataFrame):
            chunk = {k: v.to_numpy() for k, v in chunk.items()}
        missing = [col for col in self.columns if col not in chunk]
        if missing:
            raise ValueError(f"Pixel columns not found in input: {missing}")
        return {col: np.asarray(chunk[col]) for col in self.columns}

    def _validate(self, chunk: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
        bin1 = chunk["bin1_id"]
        bin2 = chunk["bin2_id"]
        _check_pixel_ids(bin1, bin2, self.n_bins, self.boundscheck, self.triucheck)

        keys = bin1.astype(np.int64) * self.n_bins + bin2.astype(np.int64)
        if self.ensure_sorted:
            order = np.argsort(keys, kind="stable")
            chunk = {col: values[order] for col, values in chunk.items()}
            keys = keys[order]

        # Pixels must be in order within and across chunks
        steps = np.diff(keys, prepend=self._last_key)
        bad = (steps <= 0) if self.dupcheck else (steps < 0)
        if np.any(bad):
            i = int(np.argmax(bad))
            pixel = (int(chunk["bin1_id"][i]), int(chunk["bin2_id"][i]))
            problem = "duplicate" if steps[i] == 0 else "out-of-order"
            raise BadInputError(
                f"Found a {problem} pixel {pixel} after {self.nnz + i} pixels. "
                "Pixels must be sorted by (bin1_id, bin2_id), with no "
                "duplicates. Pass `ordered=False` to `create_cooler` for "
                "input that isn't sorted."
            )
        self._last_key = int(keys[-1])
        return chunk

    def write(self, chunk: Any) -> None:
        """Validate a chunk of pixels and append it to the pixel table."""
        if self._done:
            raise ValueError("Writer has already finished.")
        chunk = self._as_columns(chunk)
        n = len(chunk["bin1_id"])
        if n == 0:
            return
        chunk = self._validate(chunk)
        for col in self.columns:
            _check_fits_dtype(col, chunk[col], self.dtypes[col])

        with self._open() as f:
            grp = f[self._tmp_path]["pixels"]
            for col in self.columns:
                dset = grp[col]
                dset.resize((self.nnz + n,))
                dset[self.nnz : self.nnz + n] = chunk[col]
        self.nnz += n
        if "count" in chunk:
            self.sum += chunk["count"].sum()

    def finish(self) -> None:
        """Write the indexes and metadata, and move the cooler into place."""
        if self._done:
            return
        try:
            with self._open() as f:
                grp = f[self._tmp_path]
                for col in self.columns:
                    grp["pixels"][col].resize((self.nnz,))
                chrom_offset = index_bins(grp["bins"], len(self.chroms), self.n_bins)
                bin1_offset = index_pixels(grp["pixels"], self.n_bins, self.nnz)
                write_indexes(
                    grp.create_group("indexes"), chrom_offset, bin1_offset,
                    self.h5opts,
                )
                write_info(grp, self._info())
                self._move_into_place(f)
        except BaseException:
            self._done = False
            self.abort()
            raise
        self._close_file()
        self._done = True
        if self._replace_path is not None:
            os.replace(self.file_path, self._replace_path)
            self.file_path = self._replace_path

    def abort(self) -> None:
        """Remove the partially written cooler."""
        if self._done:
            return
        try:
            if self._replace_path is None:
                with self._open() as f:
                    if self._tmp_path in f:
                        del f[self._tmp_path]
        finally:
            self._close_file()
            self._done = True
            if self._replace_path is not None and op.exists(self.file_path):
                os.remove(self.file_path)

    def _close_file(self) -> None:
        if self._file is not None:
            self._file.close()
            self._file = None

    def _info(self) -> dict[str, Any]:
        info = {
            "bin-type": "fixed" if self.binsize is not None else "variable",
            "bin-size": self.binsize if self.binsize is not None else "null",
            "storage-mode": "symmetric-upper" if self.symmetric_upper else "square",
            "nchroms": len(self.chroms),
            "nbins": self.n_bins,
            "sum": self.sum,
            "nnz": self.nnz,
        }
        if self.assembly is not None:
            info["genome-assembly"] = self.assembly
        if self.metadata is not None:
            info["metadata"] = self.metadata
        return info

    def _move_into_place(self, f: h5py.File) -> None:
        tmp = f[self._tmp_path]
        if self.group_path == "/":
            for name in _TABLES:
                if name in f:
                    del f[name]
            for name in list(tmp):
                f.move(posixpath.join(self._tmp_path, name), "/" + name)
            f.attrs.update(tmp.attrs)
            del f[self._tmp_path]
        else:
            if f.get(self.group_path, getlink=True) is not None:
                del f[self.group_path]
            parent = posixpath.dirname(self.group_path.rstrip("/"))
            f.require_group(parent)
            f.move(self._tmp_path, self.group_path)

    def __enter__(self) -> _CoolerWriter:  # noqa: PYI034
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        if exc_type is None:
            self.finish()
        else:
            self.abort()


# Hidden column counting the occurrences of each pixel in unordered input
_OCCURRENCES = "__occurrences"


class _UnorderedCoolerWriter:
    """
    Write a cooler from pixel chunks in any order, using an external sort.

    Each chunk passed to :meth:`write` is sorted and saved as a "run" in a
    temporary file. :meth:`finish` merges the runs, combining the values of
    pixels that occur in more than one run according to ``agg``, and writes
    the result with a :class:`_CoolerWriter`.

    Aggregations ``"sum"`` (the default), ``"min"``, ``"max"`` and
    ``"count"`` combine directly across runs. ``"mean"`` is supported by
    carrying the sum and the number of values of the column through the runs,
    in the hidden columns named by :func:`~cooler.create._ingest.mean_parts`.
    Chunks can provide these columns (see ``aggregate_records(...,
    split_means=True)``); otherwise each occurrence of a pixel counts as one
    value. Other aggregations raise an error if there is more than one run.

    Parameters
    ----------
    cool_uri, bins, columns, dtypes
        As for :class:`_CoolerWriter`.
    agg : dict, optional
        Mapping of value column names to aggregations.
    duplicates : {"combine", "error"}, optional
        What a pixel that occurs more than once in the input, within a chunk
        or across chunks, means. "combine" (default): its values are combined
        with ``agg``, as for contributions from different chunks of records.
        "error": the input must list each pixel at most once.
    mergebuf : int, optional
        Maximum number of pixels to hold in memory while merging.
    max_merge : int, optional
        Maximum number of runs to merge at once. More runs are merged in two
        passes.
    temp_dir : str, optional
        Directory for temporary files. Defaults to the output's directory;
        ``"-"`` uses the system default.
    delete_temp : bool, optional
        Delete temporary files when done.
    **kwargs
        Other options for the output :class:`_CoolerWriter`. The validation
        options also apply to each chunk.

    """

    def __init__(
        self,
        cool_uri: str,
        bins: pd.DataFrame,
        *,
        columns: list[str] | None = None,
        dtypes: dict | None = None,
        agg: dict[str, str] | None = None,
        duplicates: str = "combine",
        mergebuf: int = 20_000_000,
        max_merge: int = 200,
        temp_dir: str | None = None,
        delete_temp: bool = True,
        **kwargs,
    ):
        from ._ingest import mean_parts

        if duplicates not in ("combine", "error"):
            raise ValueError(
                f"duplicates must be 'combine' or 'error', not '{duplicates}'."
            )
        self.duplicates = duplicates

        allowed = set(inspect.signature(_CoolerWriter).parameters) - {
            "cool_uri", "bins", "columns", "dtypes"
        }
        unknown = sorted(set(kwargs) - allowed)
        if unknown:
            raise TypeError(f"Unexpected keyword arguments: {', '.join(unknown)}")

        self.cool_uri = cool_uri
        self.bins = bins
        self.columns, self.dtypes = _pixel_schema(columns, dtypes)
        self.writer_kwargs = kwargs
        self.mergebuf = mergebuf
        self.max_merge = max_merge

        value_cols = [c for c in self.columns if c not in ("bin1_id", "bin2_id")]
        self.agg = {col: "sum" for col in value_cols}
        if agg is not None:
            unknown = sorted(set(agg) - set(value_cols))
            if unknown:
                raise ValueError(f"Aggregations given for unknown columns: {unknown}")
            self.agg.update(agg)
        self.mean_cols = [col for col, func in self.agg.items() if func == "mean"]
        self.uncombinable = [
            col
            for col, func in self.agg.items()
            if func != "mean" and func not in _MERGEABLE_AGGS
        ]

        # Schema of the runs: mean columns are replaced by their parts
        self.run_columns, self.run_dtypes = [], {}
        self.merge_agg = {}
        for col in self.columns:
            if col in self.mean_cols:
                sum_col, n_col = mean_parts(col)
                self.run_columns += [sum_col, n_col]
                self.run_dtypes.update({sum_col: np.float64, n_col: np.int64})
                self.merge_agg.update({sum_col: "sum", n_col: "sum"})
            else:
                self.run_columns.append(col)
                self.run_dtypes[col] = self.dtypes[col]
                if col in self.agg:
                    self.merge_agg[col] = _MERGEABLE_AGGS.get(self.agg[col], "sum")
        if duplicates == "error":
            # Count occurrences of each pixel, to detect duplicates across runs
            self.run_columns.append(_OCCURRENCES)
            self.run_dtypes[_OCCURRENCES] = np.int64
            self.merge_agg[_OCCURRENCES] = "sum"
        self.run_kwargs = {
            key: kwargs[key]
            for key in ("symmetric_upper", "h5opts", "boundscheck", "triucheck")
            if key in kwargs
        }
        self.run_kwargs["dupcheck"] = duplicates == "error"

        if temp_dir is None:
            temp_dir = op.dirname(parse_cooler_uri(cool_uri)[0])
        elif temp_dir == "-":
            temp_dir = None  # makes tempfile module use the system dir
        self.temp_dir = temp_dir
        self.delete_temp = delete_temp
        # Temporary files can't be deleted while open on Windows
        self._delete_on_close = delete_temp and os.name != "nt"
        self._temp_files = []
        self._runs = []
        self._n_chunks = 0
        self._done = False
        self._run_file = self._new_temp_file()

    def _new_temp_file(self) -> Any:
        tf = tempfile.NamedTemporaryFile(
            suffix=".multi.cool", delete=self._delete_on_close, dir=self.temp_dir
        )
        self._temp_files.append(tf)
        return tf

    def _as_run(self, chunk: Any) -> dict[str, np.ndarray]:
        from ._ingest import mean_parts

        if isinstance(chunk, pd.DataFrame):
            chunk = {k: v.to_numpy() for k, v in chunk.items()}
        chunk = dict(chunk)
        for col in self.mean_cols:
            sum_col, n_col = mean_parts(col)
            if sum_col in chunk and n_col in chunk:
                continue
            if col not in chunk:
                raise ValueError(f"Pixel columns not found in input: ['{col}']")
            values = np.asarray(chunk.pop(col), dtype=float)
            missing = np.isnan(values)
            chunk[sum_col] = np.where(missing, 0.0, values)
            chunk[n_col] = (~missing).astype(np.int64)
        if self.duplicates == "error":
            chunk[_OCCURRENCES] = np.ones(len(chunk["bin1_id"]), dtype=np.int64)
        return chunk

    def write(self, chunk: Any) -> None:
        """Sort a chunk of pixels and save it as a run."""
        if self._done:
            raise ValueError("Writer has already finished.")
        chunk = self._as_run(chunk)
        uri = f"{self._run_file.name}::{self._n_chunks}"
        logger.info(f"Writing chunk {self._n_chunks}: {uri}")
        self._n_chunks += 1
        with _CoolerWriter(
            uri,
            self.bins,
            columns=self.run_columns,
            dtypes=self.run_dtypes,
            mode="a",
            ensure_sorted=True,
            **self.run_kwargs,
        ) as w:
            w.write(chunk)
        if w.nnz > 0:
            self._runs.append(uri)

    def _merged(self, uris: list[str]) -> Any:
        from .._reduce import CoolerMerger
        from ..api import Cooler

        return CoolerMerger(
            [Cooler(uri) for uri in uris],
            self.mergebuf,
            columns=[c for c in self.run_columns if c not in ("bin1_id", "bin2_id")],
            agg=self.merge_agg,
        )

    def _finish_chunk(self, chunk: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
        from ._ingest import mean_parts

        chunk = dict(chunk)
        if self.duplicates == "error":
            occurrences = chunk.pop(_OCCURRENCES)
            if np.any(occurrences > 1):
                i = int(np.argmax(occurrences > 1))
                pixel = (int(chunk["bin1_id"][i]), int(chunk["bin2_id"][i]))
                raise BadInputError(
                    f"Found a duplicate pixel {pixel}, which occurs "
                    f"{int(occurrences[i])} times in the input."
                )
        for col in self.mean_cols:
            sum_col, n_col = mean_parts(col)
            total, n = chunk.pop(sum_col), chunk.pop(n_col)
            with np.errstate(invalid="ignore", divide="ignore"):
                chunk[col] = np.where(n > 0, total / n, np.nan)
        return chunk

    def finish(self) -> None:
        """Merge the runs and write the output cooler."""
        if self._done:
            return
        try:
            if len(self._runs) > 1 and self.uncombinable:
                details = ", ".join(
                    f"'{self.agg[c]}' for '{c}'" for c in self.uncombinable
                )
                raise ValueError(
                    f"Aggregation {details} can't be combined across chunks of "
                    "input. Use 'sum', 'min', 'max', 'count' or 'mean', or a "
                    "chunk size large enough for the whole input to fit in one "
                    "chunk."
                )

            runs = self._runs
            if len(runs) > self.max_merge > 0:
                # Merge in two passes: groups of ~sqrt(n) runs, keeping both
                # the group size and the number of groups within max_merge
                # where possible.
                n = len(runs)
                group_size = max(
                    int(np.ceil(np.sqrt(n))), int(np.ceil(n / self.max_merge))
                )
                group_size = max(2, min(group_size, self.max_merge))
                tf2 = self._new_temp_file()
                merged_runs = []
                for lo in range(0, n, group_size):
                    hi = min(lo + group_size, n)
                    uri = f"{tf2.name}::{lo}-{hi}"
                    logger.info(f"Merging chunks {lo}-{hi}: {uri}")
                    with _CoolerWriter(
                        uri, self.bins, columns=self.run_columns,
                        dtypes=self.run_dtypes, mode="a", **self.run_kwargs,
                    ) as w:
                        for chunk in self._merged(runs[lo:hi]):
                            w.write(chunk)
                    merged_runs.append(uri)
                runs = merged_runs

            logger.info(f"Merging into {self.cool_uri}")
            with _CoolerWriter(
                self.cool_uri,
                self.bins,
                columns=self.columns,
                dtypes=self.dtypes,
                **self.writer_kwargs,
            ) as w:
                if runs:
                    for chunk in self._merged(runs):
                        w.write(self._finish_chunk(chunk))
        finally:
            self._cleanup()

    def abort(self) -> None:
        """Remove the temporary files without writing the output."""
        if not self._done:
            self._cleanup()

    def _cleanup(self) -> None:
        self._done = True
        for tf in self._temp_files:
            if not tf.closed:
                tf.close()
            if self.delete_temp and not self._delete_on_close and op.exists(tf.name):
                os.remove(tf.name)
        self._temp_files = []

    def __enter__(self) -> _UnorderedCoolerWriter:  # noqa: PYI034
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        if exc_type is None:
            self.finish()
        else:
            self.abort()


class _ScoolWriter:
    """
    Write a single-cell (scool) file, one cell at a time.

    The shared chrom and bin tables are written when the writer is created.
    Each call to :meth:`add_cell` writes one cell as a cooler under
    ``cells/<name>`` that hard-links the shared tables, and :meth:`finish`
    writes the file's metadata. Each cell is written atomically, so a failure
    while adding a cell leaves the cells already added intact.

    In mode "a", cells can be added to an existing scool file with the same
    bin table, and adding a cell with an existing name replaces it.

    Parameters
    ----------
    scool_uri : str
        Destination file path, optionally followed by ``::`` and a group path.
    bins : pandas.DataFrame
        The shared bin table: ``chrom``, ``start`` and ``end`` columns. Any
        other columns are ignored here; per-cell bin columns can be passed to
        :meth:`add_cell`.
    columns, dtypes, symmetric_upper, h5opts, boundscheck, triucheck, dupcheck,
    ensure_sorted, lock
        As for :class:`_CoolerWriter`, applied to every cell.
    metadata : dict, optional
        JSON-serializable metadata, stored on the file and on every cell.
    assembly : str, optional
        Name of the genome assembly.
    mode : {"w", "a"}, optional
        "w" truncates the file; "a" adds cells to an existing scool.

    """

    def __init__(
        self,
        scool_uri: str,
        bins: pd.DataFrame,
        *,
        columns: list[str] | None = None,
        dtypes: dict | None = None,
        metadata: dict | None = None,
        assembly: str | None = None,
        symmetric_upper: bool = True,
        mode: str = "w",
        h5opts: dict | None = None,
        boundscheck: bool = True,
        triucheck: bool | None = None,
        dupcheck: bool = True,
        ensure_sorted: bool = False,
        lock: Any | None = None,
    ):
        file_path, group_path = parse_cooler_uri(scool_uri)
        self.file_path = op.realpath(file_path)
        self.group_path = group_path
        self.scool_uri = f"{self.file_path}::{group_path}"
        if not isinstance(bins, pd.DataFrame):
            raise ValueError(  # noqa: TRY004
                "The bin table must be a pandas DataFrame."
            )
        self.bins, self.chroms, self.binsize = _prepare_tables(
            bins[["chrom", "start", "end"]]
        )
        self.metadata = metadata
        self.assembly = assembly
        self.h5opts = _set_h5opts(h5opts)
        self.cell_kwargs = {
            "columns": columns,
            "dtypes": dtypes,
            "metadata": metadata,
            "assembly": assembly,
            "symmetric_upper": symmetric_upper,
            "h5opts": h5opts,
            "boundscheck": boundscheck,
            "triucheck": triucheck,
            "dupcheck": dupcheck,
            "ensure_sorted": ensure_sorted,
            "lock": lock,
        }
        self._added = set()
        self._done = False

        with h5py.File(self.file_path, mode) as f:
            logger.info(f'Creating scool at "{self.file_path}::{group_path}"')
            root = f.require_group(group_path)
            # Cells hard-link the root chroms and bins tables, so when
            # appending, keep the existing tables rather than rewriting them.
            if _same_chroms_and_bins(root, self.chroms, self.bins):
                return
            if "cells" in root and len(root["cells"]):
                raise ValueError(
                    "Cannot append cells with a different bin table to an "
                    f"existing scool file: '{self.file_path}'."
                )
            for name in ["chroms", "bins"]:
                if name in root:
                    del root[name]
            write_chroms(root.create_group("chroms"), self.chroms, self.h5opts)
            write_bins(
                root.create_group("bins"), self.bins, self.chroms["name"],
                self.h5opts,
            )

    def add_cell(
        self,
        name: str,
        pixels: Any,
        *,
        bins: pd.DataFrame | None = None,
        ordered: bool = True,
        agg: dict[str, str] | None = None,
        mergebuf: int = 20_000_000,
        max_merge: int = 200,
        temp_dir: str | None = None,
        delete_temp: bool = True,
    ) -> None:
        """
        Write a cell.

        Parameters
        ----------
        name : str
            Name of the cell. It can't contain "/".
        pixels : DataFrame, dict or iterable of chunks
            The cell's pixels. A DataFrame or dict of columns is sorted in
            memory. An iterable of chunks must be sorted, unless ``ordered``
            is False.
        bins : pandas.DataFrame, optional
            Bin table with additional per-cell columns (e.g. ``weight``). Its
            ``chrom``, ``start`` and ``end`` columns must match the shared
            bin table.
        ordered : bool, optional
            Whether the chunks of ``pixels`` are sorted. If False, they are
            sorted with an external sort, using ``agg``, ``mergebuf``,
            ``max_merge``, ``temp_dir`` and ``delete_temp`` as for
            :func:`cooler.create_cooler`.

        """
        if self._done:
            raise ValueError("Writer has already finished.")
        if not name or "/" in name:
            raise ValueError(f"Invalid cell name: '{name}'.")
        if name in self._added:
            raise ValueError(f"Cell '{name}' has already been added.")

        cell_bins = self.bins
        if bins is not None:
            shared = bins[["chrom", "start", "end"]].reset_index(drop=True)
            if not (
                len(shared) == len(self.bins)
                and np.array_equal(
                    shared["chrom"].astype(str), self.bins["chrom"].astype(str)
                )
                and np.array_equal(shared["start"], self.bins["start"])
                and np.array_equal(shared["end"], self.bins["end"])
            ):
                raise ValueError(
                    f"The bin table of cell '{name}' doesn't match the shared "
                    "bin table."
                )
            cell_bins = bins.reset_index(drop=True)

        if isinstance(pixels, (pd.DataFrame, dict)):
            pixels = [pd.DataFrame(pixels).sort_values(["bin1_id", "bin2_id"])]
            ordered = True

        cell_uri = f"{self.file_path}::" + posixpath.join(
            self.group_path, "cells", name
        )
        kwargs = dict(self.cell_kwargs, mode="a", link_tables_from=self.scool_uri)
        if ordered:
            writer = _CoolerWriter(cell_uri, cell_bins, **kwargs)
        else:
            writer = _UnorderedCoolerWriter(
                cell_uri,
                cell_bins,
                agg=agg,
                mergebuf=mergebuf,
                max_merge=max_merge,
                temp_dir=temp_dir,
                delete_temp=delete_temp,
                **kwargs,
            )
        with writer:
            for chunk in pixels:
                writer.write(chunk)
        self._added.add(name)

    def finish(self) -> None:
        """Write the scool file's metadata."""
        if self._done:
            return
        self._done = True
        info = {
            "bin-type": "fixed" if self.binsize is not None else "variable",
            "bin-size": self.binsize if self.binsize is not None else "null",
            "nchroms": len(self.chroms),
            "nbins": len(self.bins),
        }
        if self.assembly is not None:
            info["genome-assembly"] = self.assembly
        if self.metadata is not None:
            info["metadata"] = self.metadata
        with h5py.File(self.file_path, "r+") as f:
            root = f[self.group_path]
            # Count all cells, including those from earlier appends
            info["ncells"] = len(root["cells"]) if "cells" in root else 0
            write_info(root, info, scool=True)

    def abort(self) -> None:
        """Stop without writing the metadata. Cells already added are kept."""
        self._done = True

    def __enter__(self) -> _ScoolWriter:  # noqa: PYI034
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        if exc_type is None:
            self.finish()
        else:
            self.abort()


# Columns that define a table and can't be added or replaced as data columns
_CORE_COLUMNS = {
    "chroms": {"name", "length"},
    "bins": {"chrom", "start", "end"},
    "pixels": {"bin1_id", "bin2_id", "count"},
}


class _ColumnWriter:
    """
    Add columns to a table of an existing cooler, in batches of rows.

    Each call to :meth:`write` takes a batch with every new column, for a
    range of rows. Batches are written in order by default, or at an explicit
    row offset. :meth:`finish` checks that every row was written and only then
    adds the columns to the table, replacing existing ones if ``replace`` is
    True, so the table never has partially written columns and replaced
    columns survive a failure.

    Parameters
    ----------
    cool_uri : str
        URI of an existing cooler.
    table : {"chroms", "bins", "pixels"}
        Table to add columns to.
    columns : dict
        Mapping of new column names to dtypes. The dtype can be None, for
        float64.
    attrs : dict, optional
        Mapping of column names to HDF5 attributes to store on the column.
    replace : bool, optional
        Replace existing columns with the same names. Otherwise, adding an
        existing column is an error.
    h5opts : dict, optional
        HDF5 dataset storage options.
    lock : multiprocessing.Lock, optional
        If given, the file is only opened while holding the lock, for each
        write. Otherwise, it is kept open until the writer is finished.

    """

    def __init__(
        self,
        cool_uri: str,
        table: str,
        columns: dict[str, Any],
        *,
        attrs: dict[str, dict] | None = None,
        replace: bool = False,
        h5opts: dict | None = None,
        lock: Any | None = None,
    ):
        if table not in _CORE_COLUMNS:
            raise ValueError(
                f"Unknown table '{table}'. Use one of {sorted(_CORE_COLUMNS)}."
            )
        if not columns:
            raise ValueError("No columns given.")
        core = sorted(set(columns) & _CORE_COLUMNS[table])
        if core:
            raise ValueError(f"Can't add or replace core columns of '{table}': {core}")

        file_path, group_path = parse_cooler_uri(cool_uri)
        self.file_path = op.realpath(file_path)
        self.table_path = posixpath.join(group_path, table)
        self._tmp_path = posixpath.join(group_path, ".__writing__columns", table)
        self.columns = {
            col: np.dtype(np.float64 if dtype is None else dtype)
            for col, dtype in columns.items()
        }
        self.attrs = attrs or {}
        self.replace = replace
        self.h5opts = _set_h5opts(h5opts)
        self.lock = lock
        self._file = None
        self._done = False
        self._cursor = 0
        self._written = []  # sorted, non-overlapping (lo, hi) row ranges

        with self._open() as f:
            if self.table_path not in f:
                raise KeyError(f"No table '{self.table_path}' in '{self.file_path}'.")
            grp = f[self.table_path]
            existing = sorted(set(columns) & set(grp))
            if existing and not replace:
                raise ValueError(
                    f"Columns already exist in '{table}': {existing}. Use replace "
                    "to replace them."
                )
            first_core = min(_CORE_COLUMNS[table] & set(grp))
            self.nrows = len(grp[first_core])

            if self._tmp_path in f:  # left over from a failed run
                del f[self._tmp_path]
            tmp = f.require_group(self._tmp_path)
            h5opts = self.h5opts if self.nrows > 0 else {}
            for col, dtype in self.columns.items():
                # Tracking attribute order gives the dataset a newer object
                # header that can hold attributes larger than 64 KiB (#226)
                dset = tmp.create_dataset(
                    col, shape=(self.nrows,), dtype=dtype, track_order=True,
                    **h5opts
                )
                dset.attrs.update(self.attrs.get(col, {}))

        if lock is None:
            self._file = h5py.File(self.file_path, "r+")

    def _open(self) -> _CoolerWriter._Handle:
        return _CoolerWriter._Handle(self, "r+")

    def write(self, batch: Any, lo: int | None = None) -> None:
        """
        Write a batch of rows of the new columns.

        Parameters
        ----------
        batch : dict, pandas.DataFrame or polars.DataFrame
            Values of every new column, for consecutive rows.
        lo : int, optional
            Row offset of the batch. Default is right after the previous
            batch.

        """
        if self._done:
            raise ValueError("Writer has already finished.")
        missing = [col for col in self.columns if col not in batch]
        if missing:
            raise ValueError(f"Columns missing from batch: {missing}")
        values = {col: np.asarray(batch[col]) for col in self.columns}
        lengths = {len(v) for v in values.values()}
        if len(lengths) != 1:
            raise ValueError("Columns in a batch must have the same length.")
        n = lengths.pop()
        if n == 0:
            return

        lo = self._cursor if lo is None else int(lo)
        hi = lo + n
        if lo < 0 or hi > self.nrows:
            raise ValueError(
                f"Rows [{lo}, {hi}) are out of range for a table of {self.nrows} rows."
            )
        i = bisect.bisect_left(self._written, (lo, hi))
        if (i > 0 and self._written[i - 1][1] > lo) or (
            i < len(self._written) and self._written[i][0] < hi
        ):
            raise ValueError(f"Rows [{lo}, {hi}) overlap rows already written.")
        for col, dtype in self.columns.items():
            _check_fits_dtype(col, values[col], dtype)

        with self._open() as f:
            tmp = f[self._tmp_path]
            for col in self.columns:
                tmp[col][lo:hi] = values[col]
        self._written.insert(i, (lo, hi))
        self._cursor = hi

    def _missing_rows(self) -> list[tuple[int, int]]:
        missing, pos = [], 0
        for lo, hi in self._written:
            if lo > pos:
                missing.append((pos, lo))
            pos = hi
        if pos < self.nrows:
            missing.append((pos, self.nrows))
        return missing

    def finish(self) -> None:
        """Check that all rows were written, and add the columns to the table."""
        if self._done:
            return
        missing = self._missing_rows()
        if missing:
            self.abort()
            shown = ", ".join(f"[{lo}, {hi})" for lo, hi in missing[:3])
            more = f" and {len(missing) - 3} more" if len(missing) > 3 else ""
            raise ValueError(f"Rows were not written: {shown}{more}.")
        try:
            with self._open() as f:
                grp = f[self.table_path]
                for col in self.columns:
                    if col in grp:
                        del grp[col]
                    f.move(
                        posixpath.join(self._tmp_path, col), f"{self.table_path}/{col}"
                    )
                self._remove_tmp(f)
        finally:
            self._close_file()
            self._done = True

    def abort(self) -> None:
        """Discard the new columns."""
        if self._done:
            return
        try:
            with self._open() as f:
                self._remove_tmp(f)
        finally:
            self._close_file()
            self._done = True

    def _remove_tmp(self, f: h5py.File) -> None:
        if self._tmp_path in f:
            del f[self._tmp_path]
        parent = posixpath.dirname(self._tmp_path)
        if parent in f and len(f[parent]) == 0:
            del f[parent]

    def _close_file(self) -> None:
        if self._file is not None:
            self._file.close()
            self._file = None

    def __enter__(self) -> _ColumnWriter:  # noqa: PYI034
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        if exc_type is None:
            self.finish()
        else:
            self.abort()


def _add_columns(
    cool_uri: str,
    table: str,
    data: Any,
    *,
    dtypes: dict[str, Any] | None = None,
    attrs: dict[str, dict] | None = None,
    replace: bool = False,
    h5opts: dict | None = None,
) -> None:
    """
    Add columns held in memory to a table of an existing cooler.

    ``data`` is a dict of arrays, or a pandas or polars DataFrame, with one
    value per row of the table. Dtypes default to those of the data.
    """
    names = list(data.columns) if hasattr(data, "columns") else list(data.keys())
    dtypes = dtypes or {}
    columns = {col: dtypes.get(col, np.asarray(data[col]).dtype) for col in names}
    with _ColumnWriter(
        cool_uri, table, columns, attrs=attrs, replace=replace, h5opts=h5opts
    ) as writer:
        writer.write(data)
