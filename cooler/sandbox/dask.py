from __future__ import division, print_function, absolute_import
from math import ceil

import numpy as np
import pandas as pd
import h5py
import cooler

import dask
from dask.base import tokenize
import dask.dataframe as dd
import dask.array as da


def get_group_info(path, grouppath, keys):
    with h5py.File(path, 'r') as f:
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
            {key: np.array([], dtype=grp[key].dtype) for key in keys},
            columns=keys)

        for key in categoricals:
            meta[key] = pd.Categorical([],
                categories=categoricals[key], ordered=True)

    return nrows, keys, meta, categoricals


def slice_dataset(filepath, grouppath, key, slc, lock=None):
    try:
        if lock is not None:
            lock.acquire()
        with h5py.File(filepath, 'r') as f:
            return f[grouppath][key][slc]
    finally:
        if lock is not None:
            lock.release()


def slice_group(filepath, grouppath, keys, slc, lock=None):
    try:
        if lock is not None:
            lock.acquire()
        with h5py.File(filepath, 'r') as f:
            return {key: f[grouppath][key][slc] for key in keys}
    finally:
        if lock is not None:
            lock.release()


def restore_categories(data, categorical_columns):
    for key, category_dict in categorical_columns.items():
        data[key] = pd.Categorical.from_codes(
                data[key],
                category_dict,
                ordered=True)
    return data


def read_table(group_uri, keys=None, chunksize=int(10e6), index=None,
            lock=None):
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
    :py:class:`dask.dataframe.DataFrame`

    Notes
    -----
    Learn more about the `dask <https://docs.dask.org/en/latest/>`_ project.

    """
    filepath, grouppath = cooler.util.parse_cooler_uri(group_uri)
    nrows, keys, meta, categoricals = get_group_info(filepath, grouppath, keys)

    # Make a unique task name
    token = tokenize(filepath, grouppath, chunksize, keys)
    task_name = 'daskify-h5py-table-' + token

    # Partition the table
    divisions = (0,) + tuple(range(-1, nrows, chunksize))[1:]
    if divisions[-1] != nrows - 1:
        divisions = divisions + (nrows - 1,)

    # Build the task graph
    dsk = {}
    for i in range(0, int(ceil(nrows / chunksize))):
        slc = slice(i * chunksize, (i + 1) * chunksize)
        data_dict = (slice_group, filepath, grouppath, keys, slc, lock)
        if categoricals:
            data_dict = (restore_categories, data_dict, categoricals)
        dsk[task_name, i] = (pd.DataFrame, data_dict, None, meta.columns)

    # Generate ddf from dask graph
    df = dd.DataFrame(dsk, task_name, meta, divisions)
    if index is not None:
        df = df.set_index(index, sorted=True, drop=False)
    return df
