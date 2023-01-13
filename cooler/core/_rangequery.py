# import h5py
import numpy as np
import pandas as pd
from cytoolz import compose


def _region_to_extent(h5, chrom_ids, region, binsize):
    chrom, start, end = region
    cid = chrom_ids[chrom]
    if binsize is not None:
        chrom_offset = h5["indexes"]["chrom_offset"][cid]
        yield chrom_offset + int(np.floor(start / binsize))
        yield chrom_offset + int(np.ceil(end / binsize))
    else:
        chrom_lo = h5["indexes"]["chrom_offset"][cid]
        chrom_hi = h5["indexes"]["chrom_offset"][cid + 1]
        chrom_bins = h5["bins"]["start"][chrom_lo:chrom_hi]
        yield chrom_lo + chrom_lo.dtype.type(
            np.searchsorted(chrom_bins, start, "right") - 1
        )
        yield chrom_lo + chrom_lo.dtype.type(np.searchsorted(chrom_bins, end, "left"))


def region_to_offset(h5, chrom_ids, region, binsize=None):
    return next(_region_to_extent(h5, chrom_ids, region, binsize))


def region_to_extent(h5, chrom_ids, region, binsize=None):
    return tuple(_region_to_extent(h5, chrom_ids, region, binsize))


def _comes_before(a0, a1, b0, b1, strict=False):
    if a0 < b0:
        return a1 <= b0 if strict else a1 <= b1
    return False


def _contains(a0, a1, b0, b1, strict=False):
    if a0 > b0 or a1 < b1:
        return False
    if strict and (a0 == b0 or a1 == b1):
        return False
    return a0 <= b0 and a1 >= b1


def concat(*dcts):
    if not dcts:
        return {}
    return {key: np.concatenate([dct[key] for dct in dcts]) for key in dcts[0]}


def transpose(dct):
    x, y = dct["bin1_id"], dct["bin2_id"]
    dct["bin1_id"], dct["bin2_id"] = y, x
    return dct


def frame_slice_from_dict(dct, field):
    index = dct.get("__index")
    return pd.DataFrame(dct, columns=["bin1_id", "bin2_id", field], index=index)


def sparray_slice_from_dict(dct, row_start, row_stop, col_start, col_stop, field):
    from sparse import COO

    shape = (row_stop - row_start, col_stop - col_start)
    return COO(
        (dct["bin1_id"] - row_start, dct["bin2_id"] - col_start),
        dct[field],
        shape=shape,
    )


def spmatrix_slice_from_dict(dct, row_start, row_stop, col_start, col_stop, field):
    from scipy.sparse import coo_matrix

    shape = (row_stop - row_start, col_stop - col_start)
    return coo_matrix(
        (dct[field], (dct["bin1_id"] - row_start, dct["bin2_id"] - col_start)),
        shape=shape,
    )


def array_slice_from_dict(dct, row_start, row_stop, col_start, col_stop, field):
    mat = spmatrix_slice_from_dict(dct, row_start, row_stop, col_start, col_stop, field)
    return mat.toarray()


def arg_prune_partition(seq, step):
    """
    Take a monotonic sequence of integers and downsample it such that they
    are at least ``step`` apart (roughly), preserving the first and last
    elements. Returns indices, not values.

    """
    lo, hi = seq[0], seq[-1]
    num = 2 + (hi - lo) // step
    cuts = np.linspace(lo, hi, num, dtype=int)
    return np.unique(np.searchsorted(seq, cuts))


class CSRReader:
    """
    Process full or partial 2D range queries from a CSR matrix stored as a group
    of columns.

    Parameters
    ----------
    pixel_grp : h5py.Group or dict-like
        Pixel group with keys {'bin1_id', 'bin2_id'}.
    bin1_offsets : 1D array-like
        The offsets of each bin1 in the pixel table (aka indptr).

    """

    def __init__(
        self,
        pixel_grp,
        bin1_offsets,
    ):
        # TODO: replace with file_path/handle, pixel_table_path
        self.pixel_grp = pixel_grp
        self.dtypes = {col: pixel_grp[col].dtype for col in pixel_grp}
        self.bin1_offsets = bin1_offsets

    def get_spans(self, bbox, chunksize):
        # Prune away (downsample) some bin1 offsets so that we extract big
        # enough chunks of matrix rows at a time.
        i0, i1, j0, j1 = bbox
        if (i1 - i0 < 1) or (j1 - j0 < 1):
            edges = np.array([], dtype=int)
        else:
            edges = i0 + arg_prune_partition(self.bin1_offsets[i0 : i1 + 1], chunksize)
        return list(zip(edges[:-1], edges[1:]))

    def get_dict_meta(self, field, return_index=False):
        dct = {
            "bin1_id": np.empty((0,), dtype=self.dtypes["bin1_id"]),
            "bin2_id": np.empty((0,), dtype=self.dtypes["bin2_id"]),
            field: np.empty((0,), dtype=self.dtypes[field]),
        }
        if return_index:
            dct["__index"] = np.empty((0,), dtype=np.int64)
        return dct

    def get_frame_meta(self, field):
        return pd.DataFrame(self.get_dict_meta(field))

    def __call__(self, field, bbox, row_span=None, reflect=False, return_index=False):
        """
        Materialize a sparse 2D range query as a dict.

        Parameters
        ----------
        field : str
            Name of value column to fetch from.

        bbox : 4-tuple
            Bounding box of the range query
            (row_start, row_stop, col_start, col_stop)

        row_span : 2-tuple, optional
            A subinterval of the bbox row span to process. If not provided, use
            all of (bbox[0], bbox[1]).

        reflect : bool, optional
            If the query bounding box covers parts of both upper and lower
            triangles of the parent matrix, reflect (copy) the pixels in the
            upper triangle part to the lower triangle part. Note that this only
            applies to the data within the bounding box. [Default: False]

        return_index : bool, optional
            Return the index values from the pixel table. Reflected elements
            carry the same index as the pixels they were reflected from. Stored
            using extra dictionary key "__index".

        Returns
        -------
        dict of columns with keys {'bin_id', 'bin2_id', field}

        """
        i0, i1, j0, j1 = bbox
        if row_span is None:
            s0, s1 = i0, i1
        else:
            s0, s1 = row_span

        # Initialize output dictionary
        result = {"bin1_id": [], "bin2_id": [], field: []}
        if return_index:
            result["__index"] = []

        # Find the offsets of our row limits in the pixel table.
        offset_lo, offset_hi = self.bin1_offsets[s0], self.bin1_offsets[s1]
        slc = slice(offset_lo, offset_hi)

        # TODO: open file in context manager in here
        bin1_selector = self.pixel_grp["bin1_id"]
        bin2_selector = self.pixel_grp["bin2_id"]
        data_selector = self.pixel_grp[field]
        # Extract the j coordinates and values of the pixels
        bin2_extracted = bin2_selector[slc]
        data_extracted = data_selector[slc]

        # Optionally, include the true index values from the pixel table.
        if return_index:
            ind_extracted = np.arange(slc.start, slc.stop)

        # Now, go row by row, filter out unwanted columns, and accumulate
        # the results.
        for i in range(s0, s1):
            # Shift the global offsets to relative ones.
            lo = self.bin1_offsets[i] - offset_lo
            hi = self.bin1_offsets[i + 1] - offset_lo

            # Get the j coordinates for this row and filter for the range
            # of j values we want.
            bin2 = bin2_extracted[lo:hi]
            mask = (bin2 >= j0) & (bin2 < j1)
            cols = bin2[mask]

            # Apply same mask to the pixel values.
            data = data_extracted[lo:hi][mask]

            # Shortcut to get i coordinates.
            rows = np.full(len(cols), i, dtype=bin1_selector.dtype)

            result["bin1_id"].append(rows)
            result["bin2_id"].append(cols)
            result[field].append(data)
            if return_index:
                result["__index"].append(ind_extracted[lo:hi][mask])

        # Concatenate outputs
        if len(result["bin1_id"]):
            for key in result.keys():
                result[key] = np.concatenate(result[key], axis=0)

            if reflect:
                to_duplex = (result["bin1_id"] != result["bin2_id"]) & (
                    result["bin2_id"] < i1
                )
                x = np.r_[result["bin1_id"], result["bin2_id"][to_duplex]]
                y = np.r_[result["bin2_id"], result["bin1_id"][to_duplex]]
                result["bin1_id"] = x
                result["bin2_id"] = y
                result[field] = np.r_[result[field], result[field][to_duplex]]

                if return_index:
                    result["__index"] = np.r_[
                        result["__index"], result["__index"][to_duplex]
                    ]
        else:
            result = self.get_dict_meta(field, return_index)

        return result


class BaseRangeQuery2D:
    """
    Mixin class for materializing 2D range queries from a sequence of
    pre-assembled tasks.

    """

    def __iter__(self):
        for task in self.tasks:
            yield task[0](*task[1:])

    @property
    def n_chunks(self):
        return len(self.tasks)

    def get(self):
        dct = concat(*self.__iter__())
        if not dct:
            return self.reader.get_dict_meta(self.field, self.return_index)
        return dct

    def get_chunk(self, i):
        if not (0 <= i < len(self.tasks)):
            raise IndexError
        task = self.tasks[i]
        return task[0](*task[1:])

    def to_delayed(self):
        from dask import delayed

        out = []
        for task in self.tasks:
            fetcher_delayed = delayed(task[0])
            out.append(fetcher_delayed(*task[1:]))
        return out

    def to_sparse_matrix(self):
        return spmatrix_slice_from_dict(self.get(), *self.bbox, self.field)

    def to_sparse_array(self):
        return sparray_slice_from_dict(self.get(), *self.bbox, self.field)

    def to_array(self):
        return array_slice_from_dict(self.get(), *self.bbox, self.field)

    def to_frame(self):
        return frame_slice_from_dict(self.get(), self.field)

    def to_dask_frame(self):
        from dask.base import tokenize
        from dask.dataframe import DataFrame

        meta = self.reader.get_frame_meta(self.field)
        tasks = self.tasks
        # spans = [task[3] for task in tasks]

        name = (
            "cooler-"
            + self.__class__.__name__
            + tokenize(self.bbox, self.field, self.return_index)
        )
        df_tasks = [(frame_slice_from_dict, task, self.field) for task in tasks]
        divisions = [None] * (len(tasks) + 1)

        keys = [(name, i) for i in range(len(df_tasks))]
        dsk = dict(zip(keys, df_tasks))

        return DataFrame(dsk, name, meta, divisions)


class DirectRangeQuery2D(BaseRangeQuery2D):
    """
    Query engine for a matrix that interprets data exactly as stored.

    Parameters
    ----------
    reader : CSRReader
        Reads pixel records from a CSR matrix.

    field : str
        Name of value column to select from.

    bbox : 4-tuple
        Query bounding box (row_start, row_stop, col_start, col_stop).
        Interpreted as half-open intervals:
        [row_start, row_stop), [col_start, col_stop).

    chunksize : int
        Rough number of pixel records to read from disk at a time (before
        filtering).

    return_index : bool, optional
        Extract the pixel table index values along with the pixels.

    Attributes
    ----------
    bbox
    chunksize
    field
    return_index
    tasks : list of tuple
        Partial query tasks represented as tuples of (callable, *args) as in Dask.
        See https://docs.dask.org/en/latest/spec.html for more details.

    """

    def __init__(self, reader, field, bbox, chunksize, return_index=False):
        self.reader = reader
        self.field = field
        self.bbox = bbox
        self.chunksize = chunksize
        self.return_index = return_index
        self.tasks = [
            (reader, field, bbox, span, False, return_index)
            for span in reader.get_spans(bbox, chunksize)
        ]


class FillLowerRangeQuery2D(BaseRangeQuery2D):
    """
    Query engine for a symmetric-upper matrix that generates the additional
    tasks required to fill in any implicit lower triangle area inside the query
    bounding box.

    Parameters
    ----------
    reader : CSRReader
        Reads pixel records from a symm-upper matrix.

    field : str
        Name of value column to select from.

    bbox : 4-tuple
        Query bounding box (row_start, row_stop, col_start, col_stop).
        Interpreted as half-open intervals:
        [row_start, row_stop), [col_start, col_stop).

    chunksize : int
        Rough number of pixel records to read from disk at a time (before
        filtering).

    return_index : bool, optional
        Extract the pixel table index values along with the pixels.

    Attributes
    ----------
    bbox
    chunksize
    field
    return_index
    tasks : list of tuple
        Partial query tasks represented as tuples of (callable, *args) as in Dask.
        See https://docs.dask.org/en/latest/spec.html for more details.

    Notes
    -----
    Each task generates a dict containing a subset of pixels from a horizontal
    or L-shaped strip of the full query bounding box. Tasks can be executed
    eagerly or lazily.

    """

    def __init__(self, reader, field, bbox, chunksize, return_index=False):
        self.reader = reader
        self.field = field
        self.bbox = bbox
        self.chunksize = chunksize
        self.return_index = return_index

        _fetch = self.reader
        _fetch_then_transpose = compose(transpose, self.reader)

        # If the lower limit of the query exceeds the right limit, we transpose
        # the query bbox to fetch data, then we transpose the result.
        i0, i1, j0, j1 = bbox
        use_transpose = i1 > j1
        if use_transpose:
            i0, i1, j0, j1 = j0, j1, i0, i1
            fetcher = _fetch_then_transpose
        else:
            fetcher = _fetch

        # Base cases:
        # Bounding box is anchored on the main diagonal or is completely off
        # the main diagonal.
        if i0 == j0 or _comes_before(i0, i1, j0, j1, strict=True):
            self._bboxes = [(i0, i1, j0, j1)]
            fetchers = [fetcher]

        # Mixed case I: partial overlap between i- and j-interval, but not
        # anchored on the main diagonal. Split the query bounding box into two
        # vertically stacked boxes.
        elif _comes_before(i0, i1, j0, j1):
            self._bboxes = [(i0, j0, j0, j1), (j0, i1, j0, j1)]
            fetchers = [fetcher, fetcher]

        # Mixed case II: i-interval nested in j-interval
        # Split the query bounding box into two horizontally stacked boxes.
        elif _contains(j0, j1, i0, i1):
            # The first block is completely in the lower triangle of the parent
            # matrix, so we query the transpose and transpose the result.
            # However, if we are already transposing, we can remove the
            # operation instead of doing it twice.
            self._bboxes = [(j0, i0, i0, i1), (i0, i1, i0, j1)]
            fetchers = [_fetch if use_transpose else _fetch_then_transpose, fetcher]

        else:
            raise ValueError("This shouldn't happen.")

        self.tasks = []
        for fetcher, bbox in zip(fetchers, self._bboxes):
            spans = self.reader.get_spans(bbox, chunksize)
            self.tasks += [
                (fetcher, field, bbox, span, True, return_index) for span in spans
            ]
