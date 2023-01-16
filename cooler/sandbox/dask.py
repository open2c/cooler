from itertools import product
from math import ceil

import h5py
import numpy as np
import pandas as pd

try:
    from sparse import COO
except ImportError:
    raise ImportError("The 'sparse' package is required to use dask") from None

import dask.array as da
import dask.dataframe as dd
from dask.base import tokenize

from ..core import CSRReader, query_rect
from ..util import parse_cooler_uri, partition


def _get_group_info(path, grouppath, keys):
    with h5py.File(path, "r") as f:
        grp = f[grouppath]

        if keys is None:
            keys = list(grp.keys())

        nrows = len(grp[keys[0]])

        categoricals = {}
        for key in keys:
            dt = h5py.check_dtype(enum=grp[key].dtype)
            if dt is not None:
                categoricals[key] = sorted(dt, key=dt.__getitem__)

        # Meta is an empty dataframe that serves as a compound "dtype"
        meta = pd.DataFrame(
            {key: np.array([], dtype=grp[key].dtype) for key in keys}, columns=keys
        )

        for key in categoricals:
            meta[key] = pd.Categorical([], categories=categoricals[key], ordered=True)

    return nrows, keys, meta, categoricals


def _slice_dataset(filepath, grouppath, key, slc, lock=None):
    try:
        if lock is not None:
            lock.acquire()
        with h5py.File(filepath, "r") as f:
            return f[grouppath][key][slc]
    finally:
        if lock is not None:
            lock.release()


def _slice_group(filepath, grouppath, keys, slc, lock=None):
    try:
        if lock is not None:
            lock.acquire()
        with h5py.File(filepath, "r") as f:
            return {key: f[grouppath][key][slc] for key in keys}
    finally:
        if lock is not None:
            lock.release()


def _restore_categories(data, categorical_columns):
    for key, category_dict in categorical_columns.items():
        data[key] = pd.Categorical.from_codes(data[key], category_dict, ordered=True)
    return data


def read_table(group_uri, keys=None, chunksize=10_000_000, index=None, lock=None):
    """
    Create a dask dataframe around a column-oriented table in HDF5.

    A table is a group containing equal-length 1D datasets.

    Parameters
    ----------
    group_uri : str
        URI to the HDF5 group storing the table.
    keys : list, optional
        list of HDF5 Dataset keys, default is to use all keys in the group
    chunksize : int, optional
        Chunk size
    index : str, optional
        Sorted column to use as index
    lock : multiprocessing.Lock, optional
        Lock to serialize HDF5 read/write access. Default is no lock.

    Returns
    -------
    :class:`dask.dataframe.DataFrame`

    Notes
    -----
    Learn more about the `dask <https://docs.dask.org/en/latest/>`_ project.

    """
    filepath, grouppath = parse_cooler_uri(group_uri)
    nrows, keys, meta, categoricals = _get_group_info(filepath, grouppath, keys)

    # Make a unique task name
    token = tokenize(filepath, grouppath, chunksize, keys)
    task_name = "daskify-h5py-table-" + token

    # Partition the table
    divisions = (0,) + tuple(range(-1, nrows, chunksize))[1:]
    if divisions[-1] != nrows - 1:
        divisions = divisions + (nrows - 1,)

    # Build the task graph
    dsk = {}
    for i in range(0, int(ceil(nrows / chunksize))):
        slc = slice(i * chunksize, (i + 1) * chunksize)
        data_dict = (_slice_group, filepath, grouppath, keys, slc, lock)
        if categoricals:
            data_dict = (_restore_categories, data_dict, categoricals)
        dsk[task_name, i] = (pd.DataFrame, data_dict, None, meta.columns)

    # Generate ddf from dask graph
    df = dd.DataFrame(dsk, task_name, meta, divisions)
    if index is not None:
        df = df.set_index(index, sorted=True, drop=False)
    return df


def _array_select(clr, i0, i1, j0, j1, field, sparse_array):
    is_upper = clr._is_symm_upper
    with clr.open("r") as h5:
        dtype = h5['pixels'][field].dtype
        reader = CSRReader(h5, field, max_chunk=500000000)
        if is_upper:
            i, j, v = query_rect(reader.query, i0, i1, j0, j1, duplex=True)
        else:
            i, j, v = reader.query(i0, i1, j0, j1)
        if not len(v):
            v = v.astype(dtype)
    arr = COO((i - i0, j - j0), v, shape=(i1 - i0, j1 - j0))
    if not sparse_array:
        arr = arr.todense()
    return arr


def load_dask_array(
    clr, i0, i1, j0, j1, field="count", sparse_array=False, chunksize=256
):
    """
    Create a parallel Dask array around the matrix representation of a cooler.

    Parameters
    ----------
    clr : :class:`cooler.Cooler`
        Cooler object
    i0, i1 : int
        Row query range
    j0, j1 : int
        Column query range
    field : str
        Value column to query
    sparse_array : bool, optional
        Create a dask array backed by :class:`sparse.COO` sparse arrays
        instead of dense numpy arrays (default).
    chunksize : int, optional
        Length of the rowwise chunks to partition the underlying data into.

    Returns
    -------
    :class:`dask.array.Array`

    """
    token = tokenize(clr.uri, i0, i1, j0, i1, field, chunksize)
    task_name = "cooler-array-slice-" + token

    shape = (i1 - i0, j1 - j0)
    meta = _array_select(clr, 0, 0, 0, 0, field, sparse_array)
    slices = [(lo, hi, j0, j1) for lo, hi in partition(0, shape[0], chunksize)]
    chunks = (tuple(s[1] - s[0] for s in slices), (shape[1],))
    keys = list(product([task_name], *[range(len(dim)) for dim in chunks]))
    values = [(_array_select, clr, *slc, field, sparse_array) for slc in slices]
    dsk = dict(zip(keys, values))

    return da.Array(dsk, task_name, chunks, meta=meta, shape=shape)
