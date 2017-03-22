from __future__ import division, print_function
from operator import getitem, setitem
from math import ceil

import numpy as np
import pandas as pd
import h5py
import cooler

from dask.dataframe.core import new_dd_object
from dask.base import tokenize
import dask.dataframe as dd
import dask.array as da


def _coldict_accessor(clr, group_path, slc, columns, lock=None):
    columns = list(columns)
    if lock:
        lock.acquire()
    try:
        with clr.open() as h5:
            grp = h5[group_path]
            chunks = [grp[name][slc] for name in columns]
        data = dict(zip(columns, chunks))
    finally:
        if lock:
            lock.release()
    return data


def _restore_categories(df, col2cat):
    for col, category_dict in col2cat.items():
        df[col] = pd.Categorical.from_codes(
                df[col], sorted(category_dict, key=dt.__getitem__), ordered=True)
    return df

    
def daskify(clr, group_path, chunksize=int(10e6), columns=None, index=None, lock=None):     
    with clr.open() as h5:
        grp = h5[group_path]
        if columns is None:
            columns = tuple(grp.keys())

        shape = (len(grp[columns[0]]), len(columns))
        dtype = np.dtype([(col, grp[col].dtype) for col in columns])
        ncols = len(columns)
        nrows = shape[0]
    
        # Meta is an empty dataframe that serves as a compound "dtype"    
        meta = pd.DataFrame(
            {col: np.array([], dtype=grp[col].dtype) for col in columns}, 
            columns=columns)
    
        col2cat = {}
        for col in columns:
            dt = h5py.check_dtype(enum=grp[col].dtype)
            if dt is not None:
                meta[col] = pd.Categorical.from_codes(
                    meta[col], sorted(dt, key=dt.__getitem__), ordered=True)
                col2cat[col] = dt

    # Partition the table
    divisions = (0,) + tuple(range(-1, nrows, chunksize))[1:]             
    if divisions[-1] != nrows - 1:
        divisions = divisions + (nrows - 1,)
    
    # Make a unique task name
    token = tokenize(id(clr), group_path, chunksize, columns)
    task_name = 'cooler-daskify-' + token

    # Create the task graph
    dsk = {}
    for i in range(0, int(ceil(nrows / chunksize))):
        slc = slice(i * chunksize, (i + 1) * chunksize)
        data = (_coldict_accessor, clr, group_path, slc, columns, lock)
        if col2cat:
            data = (_restore_categories, data, col2cat)
        dsk[task_name, i] = (pd.DataFrame, data, None, meta.columns)
    
    return new_dd_object(dsk, task_name, meta, divisions)
