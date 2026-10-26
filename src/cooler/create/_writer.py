from __future__ import annotations

import inspect
import os
import os.path as op
import posixpath
import tempfile
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
        mergebuf: int = 20_000_000,
        max_merge: int = 200,
        temp_dir: str | None = None,
        delete_temp: bool = True,
        **kwargs,
    ):
        from ._create import _MERGEABLE_AGGS
        from ._ingest import mean_parts

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
        self.run_kwargs = {
            key: kwargs[key]
            for key in ("symmetric_upper", "h5opts", "boundscheck", "triucheck",
                        "dupcheck")
            if key in kwargs
        }

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

    def _finish_means(self, chunk: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
        from ._ingest import mean_parts

        chunk = dict(chunk)
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
                        w.write(self._finish_means(chunk))
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
        triucheck: bool = True,
        dupcheck: bool = True,
        ensure_sorted: bool = False,
        lock: Any | None = None,
    ):
        from ._create import _same_chroms_and_bins

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
