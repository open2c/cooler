import h5py
import numpy as np
import pandas as pd
from pandas.api.types import is_categorical_dtype


def get(grp, lo=0, hi=None, fields=None, convert_enum=True, as_dict=False):
    """
    Query a range of rows from a table as a dataframe.

    A table is an HDF5 group containing equal-length 1D datasets serving as
    columns.

    Parameters
    ----------
    grp : ``h5py.Group`` or any dict-like of array-likes
        Handle to an HDF5 group containing only 1D datasets or any similar
        collection of 1D datasets or arrays
    lo, hi : int, optional
        Range of rows to select from the table.
    fields : str or sequence of str, optional
        Column or list of columns to query. Defaults to all available columns.
        A single string returns a Series instead of a DataFrame.
    convert_enum : bool, optional
        Whether to convert HDF5 enum datasets into ``pandas.Categorical``
        columns instead of plain integer columns. Default is True.
    kwargs : optional
        Options to pass to ``pandas.DataFrame`` or ``pandas.Series``.

    Returns
    -------
    DataFrame or Series

    Notes
    -----
    HDF5 ASCII datasets are converted to Unicode.

    """
    series = False
    if fields is None:
        fields = list(grp.keys())
    elif isinstance(fields, str):
        fields = [fields]
        series = True

    data = {}
    for field in fields:
        dset = grp[field]

        if convert_enum:
            dt = h5py.check_dtype(enum=dset.dtype)
        else:
            dt = None

        if dt is not None:
            data[field] = pd.Categorical.from_codes(
                dset[lo:hi], sorted(dt, key=dt.__getitem__), ordered=True
            )
        elif dset.dtype.type == np.string_:
            data[field] = dset[lo:hi].astype("U")
        else:
            data[field] = dset[lo:hi]

    if as_dict:
        return data

    if data and lo is not None:
        index = np.arange(lo, lo + len(next(iter(data.values()))))
    else:
        index = None

    if series:
        return pd.Series(data[fields[0]], index=index, name=field)
    else:
        return pd.DataFrame(data, columns=fields, index=index)


def put(grp, df, lo=0, store_categories=True, h5opts=None):
    """
    Store a dataframe into a column-oriented table store.

    A table is an HDF5 group containing equal-length 1D datasets serving as
    columns.

    Parameters
    ----------
    h5 : ``h5py.Group``
        Handle to an HDF5 group containing only 1D datasets or any similar
        collection of 1D datasets or arrays
    df : DataFrame or Series
        Data columns to write to the HDF5 group
    lo : int, optional
        Row offset for data to be stored.
    store_categories : bool, optional
        Whether to convert ``pandas.Categorical`` columns into HDF5 enum
        datasets instead of plain integer datasets. Default is True.
    h5opts : dict, optional
        HDF5 dataset filter options to use (compression, shuffling,
        checksumming, etc.). Default is to use autochunking and GZIP
        compression, level 6.

    Notes
    -----
    Categorical data must be ASCII compatible.

    """
    if h5opts is None:
        h5opts = {"compression": "gzip", "compression_opts": 6}

    if isinstance(df, pd.Series):
        df = df.to_frame()

    # fields = df.keys()
    for field, data in df.items():

        if np.isscalar(data):
            data = np.array([data])
            dtype = data.dtype
            fillvalue = None
        elif is_categorical_dtype(data):
            if store_categories:
                cats = data.cat.categories
                enum = (data.cat.codes.dtype, dict(zip(cats, range(len(cats)))))
                data = data.cat.codes
                dtype = h5py.special_dtype(enum=enum)
                fillvalue = -1
            else:
                data = data.cat.codes
                dtype = data.dtype
                fillvalue = -1
        else:
            data = np.asarray(data)
            if data.dtype in (object, str, bytes):
                dtype = np.dtype("S")
                data = np.array(data, dtype=dtype)
                fillvalue = None
            else:
                dtype = data.dtype
                fillvalue = None

        hi = lo + len(data)
        try:
            dset = grp[field]
        except KeyError:
            dset = grp.create_dataset(
                field,
                shape=(hi,),
                dtype=dtype,
                maxshape=(None,),
                fillvalue=fillvalue,
                **h5opts
            )
        if hi > len(dset):
            dset.resize((hi,))

        dset[lo:hi] = data


def delete(grp, fields=None):
    """
    Delete columns from a table.

    A table is an HDF5 group containing equal-length 1D datasets serving as
    columns.

    Parameters
    ----------
    grp : ``h5py.Group``
        Handle to an HDF5 group containing only 1D datasets or any similar
        collection of 1D datasets or arrays
    fields : str or sequence of str, optional
        Column or list of columns to query. Defaults to all available columns.
        A single string returns a Series instead of a DataFrame.

    Notes
    -----
    Deleting objects leaves "holes" in HDF5 files and doesn't shrink the file.
    You will need to repack or copy the file contents to reclaim space.
    See the h5repack tool.

    """
    if fields is None:
        fields = list(grp.keys())
    elif isinstance(fields, str):
        fields = [fields]
    for field in fields:
        if field in grp.keys():
            del grp[field]
