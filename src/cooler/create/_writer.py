from __future__ import annotations

import os.path as op
import posixpath
import warnings
from typing import Any

import h5py
import numpy as np
import pandas as pd

from .._logging import get_logger
from ..core import put
from ..util import get_binsize, get_chromsizes, get_meta, parse_cooler_uri
from ._constants import PIXEL_DTYPES
from ._create import (
    _check_fits_dtype,
    _set_h5opts,
    index_bins,
    index_pixels,
    prepare_pixels,
    write_bins,
    write_chroms,
    write_indexes,
    write_info,
)
from ._ingest import BadInputError

logger = get_logger("cooler.create")

# Suffix of the temporary group a cooler is written into before it is moved
# into place.
_WRITING_SUFFIX = ".__writing__"

# Standard tables of a cooler group.
_TABLES = ("chroms", "bins", "pixels", "indexes")


def _prepare_tables(bins: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, Any]:
    """Return the bin table to write, the chrom table, and the bin size."""
    for col in ["chrom", "start", "end"]:
        if col not in bins.columns:
            raise ValueError(f"Missing column from bin table: '{col}'.")
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
    return bins, chroms, get_binsize(bins)


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
        Ignored if ``symmetric_upper`` is False.
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
        triucheck: bool = True,
        dupcheck: bool = True,
        ensure_sorted: bool = False,
        lock: Any | None = None,
        link_tables_from: str | None = None,
    ):
        file_path, group_path = parse_cooler_uri(cool_uri)
        self.file_path = op.realpath(file_path)
        self.group_path = group_path
        if group_path == "/":
            self._tmp_path = "/" + _WRITING_SUFFIX
        else:
            self._tmp_path = group_path.rstrip("/") + _WRITING_SUFFIX

        if not isinstance(bins, pd.DataFrame):
            raise ValueError(  # noqa: TRY004
                "The bin table must be a pandas DataFrame."
            )
        self.bins, self.chroms, self.binsize = _prepare_tables(bins)
        self.n_bins = len(self.bins)
        self.columns, self.dtypes = _pixel_schema(columns, dtypes)
        self.h5opts = _set_h5opts(h5opts)
        self.metadata = metadata
        self.assembly = assembly
        self.symmetric_upper = symmetric_upper

        if not symmetric_upper and triucheck:
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
        if self.boundscheck:
            if np.any(bin1 < 0) or np.any(bin2 < 0):
                raise BadInputError("Found bin ID < 0")
            if np.any(bin1 >= self.n_bins) or np.any(bin2 >= self.n_bins):
                raise BadInputError(
                    "Found a bin ID that exceeds the declared number of bins. "
                    "Check whether your bin table is correct."
                )
        if self.triucheck and np.any(bin1 > bin2):
            raise BadInputError("Found bin1_id greater than bin2_id")

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
        finally:
            self._close_file()
            self._done = True

    def abort(self) -> None:
        """Remove the partially written cooler."""
        if self._done:
            return
        try:
            with self._open() as f:
                if self._tmp_path in f:
                    del f[self._tmp_path]
        finally:
            self._close_file()
            self._done = True

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
