# -*- coding: utf-8 -*-
from __future__ import division, print_function
import posixpath
import tempfile
import warnings
import six

import numpy as np
import pandas
import h5py

from .core import put
from .util import get_binsize, get_chromsizes, infer_meta, get_meta
from . import get_logger


logger = get_logger()


def parse_cooler_uri(s):
    """
    Parse a Cooler URI string

    e.g. /path/to/mycoolers.cool::/path/to/cooler

    """
    parts = s.split('::')
    if len(parts) == 1:
        file_path, group_path = parts[0], '/'
    elif len(parts) == 2:
        file_path, group_path = parts
        if not group_path.startswith('/'):
            group_path = '/' + group_path
    else:
        raise ValueError("Invalid Cooler URI string")
    return file_path, group_path


def ls(filepath):
    """
    Traverse a file's data hierarchy and list all Cooler nodes.

    Parameters
    ----------
    filepath : str

    Returns
    -------
    list of Cooler group paths in the file
    
    """
    listing = []
    keys = ['chroms', 'bins', 'pixels', 'indexes']
    def _check_group(pth, grp):
        fmt = grp.attrs.get('format', None)
        url = grp.attrs.get('format-url', None)
        if fmt == MAGIC or url == URL:
            if not all(name in grp.keys() for name in keys):
                warnings.warn(
                    'Cooler path /{} appears to be corrupt'.format(pth))
            listing.append('/' + pth if not pth.startswith('/') else pth)

    with h5py.File(filepath, 'r') as f:
        _check_group('/', f)
        f.visititems(_check_group)

    return listing


def is_cooler(filepath, group=None):
    """
    Determine if a file contains a valid Cooler data hierarchy.

    Parameters
    ----------
    filepath : str
    group : str, optional
        Path to specific group to check. Otherwise returns True
        if any Cooler paths are detected.

    Returns
    -------
    bool
    
    """
    if not h5py.is_hdf5(filepath):
        return False
    if group is None:
        return len(ls(filepath)) > 0
    if not group.startswith('/'):
        group = '/' + group
    return group in ls(filepath)


def create(cool_uri, bins, pixels, metadata=None, assembly=None, dtypes='<deprecated>',
           h5opts=None, append=False, lock=None, columns=None, dtype=None,
           chromsizes='<deprecated>'):
    """
    Create a new Cooler.

    Parameters
    ----------
    cool_uri : str
        Path to Cooler file or URI to Cooler group. If the file does not exist, 
        it will be created.
    bins : pandas.DataFrame
        Segmentation of the chromosomes into genomic bins as a BED-like
        DataFrame with columns ``chrom``, ``start`` and ``end``. May contain 
        additional columns.
    pixels : DataFrame, dictionary, or iterable of either
        A table (represented as a data frame or a column-oriented dict) 
        containing columns labeled 'bin1_id', 'bin2_id' and 'count', sorted 
        by ('bin1_id', 'bin2_id'). If additional columns are included in the 
        pixel table, their names and dtypes must be specified in the ``dtypes`` 
        argument. Alternatively, for larger data, an iterable can be provided 
        that yields the pixel data as a sequence of "chunks". If the input is a
        dask DataFrame, it will also be processed one chunk at a time.
    dtypes : dict, optional
        Deprecated in favor of `dtype` to mirror the pandas constructor.
    metadata : dict, optional
        Experiment metadata to store in the file. Must be JSON compatible.
    assembly : str, optional
        Name of genome assembly.
    h5opts : dict, optional
        HDF5 dataset filter options to use (compression, shuffling,
        checksumming, etc.). Default is to use autochunking and GZIP
        compression, level 6.
    append : bool, optional
        Append new Cooler to the file if it exists. If False, an existing file
        with the same name will be truncated. Default is False.
    lock : multiprocessing.Lock, optional
        Optional lock to control concurrent access to the output file.
    columns : sequence of str, optional
        Specify here the names of any additional value columns from the input 
        besides 'count' to store in the Cooler. The standard columns ['bin1_id', 
        'bin2_id', 'count'] can be provided, but are already assumed and don't 
        need to be given explicitly. Additional value columns provided here will 
        be stored as np.float64 unless otherwised specified using `dtype`.
    dtype : dict, optional
        Dictionary mapping column names in the pixel table to dtypes. Can be 
        used to override the default dtypes of 'bin1_id', 'bin2_id' or 'count'. 
        Any additional value column dtypes must also be provided in the
        `columns` argument, or will be ignored.

    Result
    ------
    Cooler data hierarchy stored in at the specified Cooler URI (HDF5 file and 
    group path).

    """
    file_path, group_path = parse_cooler_uri(cool_uri)
    mode = 'a' if append else 'w'
    if h5opts is None:
        h5opts = dict(compression='gzip', compression_opts=6)

    if not isinstance(bins, pandas.DataFrame):
        raise ValueError(
            "Second positional argument must be a pandas DataFrame. "
            "Note that the `chromsizes` argument is now deprecated: "
            "see documentation for `create`.")

    if dtypes != '<deprecated>':
        warnings.warn("Use dtype= instead of dtypes=", FutureWarning)
        dtype = dtypes

    for col in ['chrom', 'start', 'end']:
        if col not in bins.columns:
            raise ValueError("Missing column from bin table: '{}'.".format(col))

    # Populate expected pixel column names. Include user-provided value columns.
    if columns is None:
        columns = ['bin1_id', 'bin2_id', 'count']
    else:
        columns = list(columns)
        for col in ['bin1_id', 'bin2_id', 'count']:
            if col not in columns:
                columns.insert(0, col)

    # Populate dtypes for expected pixel columns, and apply user overrides.
    if dtype is None:
        dtype = dict(PIXEL_DTYPES)
    else:
        dtype_ = dict(dtype)
        dtype = dict(PIXEL_DTYPES)
        dtype.update(dtype_)

    # Get empty "meta" header frame (assigns the undeclared dtypes).
    # Any columns from the input not in meta will be ignored.
    meta = get_meta(columns, dtype, default_dtype=float)

    # Determine the appropriate iterable
    try:
        from dask.dataframe import DataFrame as dask_df
    except (ImportError, AttributeError):
        dask_df = ()

    if isinstance(pixels, dask_df):
        iterable = map(lambda x: x.compute(), pixels.to_delayed())
        input_columns = infer_meta(pixels).columns
    elif isinstance(pixels, pandas.DataFrame):
        iterable = (pixels,)
        input_columns = infer_meta(pixels).columns
    elif isinstance(pixels, dict):
        iterable = (pixels,)
        input_columns = infer_meta(
            [(k, v.dtype) for (k, v) in pixels.items()]).columns
    else:
        iterable = pixels
        input_columns = None

    # If possible, ensure all expected columns are available
    if input_columns is not None:
        for col in columns:
            if col not in input_columns:
                col_type = 'Standard' if col in PIXEL_FIELDS else 'User'
                raise ValueError(
                    "{} column not found in input: '{}'".format(col_type, col))

    # Prepare chroms and bins
    bins = bins.copy()
    bins['chrom'] = bins['chrom'].astype(object)
    chromsizes = get_chromsizes(bins)
    try:
        chromsizes = six.iteritems(chromsizes)
    except AttributeError:
        pass
    chromnames, lengths = zip(*chromsizes)
    chroms = pandas.DataFrame(
        {'name': chromnames, 'length': lengths}, 
        columns=['name', 'length'])
    binsize = get_binsize(bins)
    n_chroms = len(chroms)
    n_bins = len(bins)

    # Create root group
    with h5py.File(file_path, mode) as f:
        logger.info('Creating cooler at "{}::{}"'.format(file_path, group_path))
        if group_path is '/':
            for name in ['chroms', 'bins', 'pixels', 'indexes']:
                if name in f:
                    del f[name]
        else:
            try:
                f.create_group(group_path)
            except ValueError:
                del f[group_path]
                f.create_group(group_path)

    # Write chroms, bins and pixels
    with h5py.File(file_path, 'r+') as f:
        h5 = f[group_path]

        logger.info('Writing chroms')
        grp = h5.create_group('chroms')
        write_chroms(grp, chroms, h5opts)

        logger.info('Writing bins')
        grp = h5.create_group('bins')
        write_bins(grp, bins, chroms['name'], h5opts)
        
        grp = h5.create_group('pixels')
        prepare_pixels(grp, n_bins, meta.columns, dict(meta.dtypes), h5opts)

    # Multiprocess HDF5 reading is supported only if the same HDF5 file is not
    # open in write mode anywhere. To read and write to the same file, pass a
    # lock shared with the HDF5 reading processes. `write_pixels` will acquire
    # it and open the file for writing for the duration of each write step
    # only. After it closes the file and releases the lock, the reading
    # processes will have to re-acquire the lock and re-open the file to obtain
    # the updated file state for reading.
    logger.info('Writing pixels')
    target = posixpath.join(group_path, 'pixels')
    nnz, ncontacts = write_pixels(
        file_path, target, meta.columns, iterable, h5opts, lock)

    # Write indexes
    with h5py.File(file_path, 'r+') as f:
        h5 = f[group_path]

        logger.info('Writing indexes')
        grp = h5.create_group('indexes')
        
        chrom_offset = index_bins(h5['bins'], n_chroms, n_bins)
        bin1_offset = index_pixels(h5['pixels'], n_bins, nnz)
        write_indexes(grp, chrom_offset, bin1_offset, h5opts)

        logger.info('Writing info')
        info = {}
        info['bin-type'] = 'fixed' if binsize is not None else 'variable'
        info['bin-size'] = binsize if binsize is not None else 'null'
        info['nchroms'] = n_chroms
        info['nbins'] = n_bins
        info['sum'] = ncontacts
        info['nnz'] = nnz
        if assembly is not None:
            info['genome-assembly'] = assembly
        if metadata is not None:
            info['metadata'] = metadata
        write_info(h5, info)

    logger.info('Done')


def append(cool_uri, table, data, chunked=False, force=False, h5opts=None, 
           lock=None):
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
    if h5opts is None:
        h5opts = dict(compression='gzip', compression_opts=6)

    file_path, group_path = parse_cooler_uri(cool_uri)

    try:
        from dask.dataframe import (
            DataFrame as dask_df,
            Series as dask_series
        )
    except (ImportError, AttributeError):
        dask_df = ()
        dask_series = ()

    if isinstance(data, dask_series):
        data = data.to_frame()

    try:
        names = data.keys()
    except AttributeError:
        names = data.columns

    with h5py.File(file_path, 'r+') as f:
        h5 = f[group_path]
        for name in names:
            if name in h5[table]:
                if not force:
                    raise ValueError(
                        "'{}' column already exists. ".format(name) +
                        "Use --force option to overwrite.")
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
            for key in data.keys():
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


def create_from_unsorted(cool_uri, bins, chunks, columns=None, dtype=None, 
                         mergebuf=int(20e6), delete_temp=True, temp_dir=None, 
                         multifile_merge=False, **kwargs):
    """
    Create a Cooler in two passes via an external sort mechanism. In the first 
    pass, a sequence of data chunks are processed and sorted in memory and saved
    to temporary Coolers. In the second pass, the temporary Coolers are merged 
    into the output.
    
    Parameters
    ----------
    cool_uri : str
        Path to Cooler file or URI to Cooler group. If the file does not exist,
        it will be created.
    bins : DataFrame
        Segmentation of the chromosomes into genomic bins. May contain 
        additional columns.
    chunks : iterable of DataFrames
        Sequence of chunks that get processed and written to separate Coolers 
        and then subsequently merged.
    columns : sequence of str, optional
        Specify here the names of any additional value columns from the input 
        besides 'count' to store in the Cooler. The standard columns ['bin1_id', 
        'bin2_id', 'count'] can be provided, but are already assumed and don't 
        need to be given explicitly. Additional value columns provided here will 
        be stored as np.float64 unless otherwised specified using `dtype`.
    dtype : dict, optional
        Dictionary mapping column names in the pixel table to dtypes. Can be 
        used to override the default dtypes of 'bin1_id', 'bin2_id' or 'count'. 
        Any additional value column dtypes must also be provided in the
        `columns` argument, or will be ignored.
    mergebuf : int, optional
        Maximum number of records to buffer in memory at any give time during 
        the merge step.
    delete_temp : bool, optional
        Whether to delete temporary files when finished. 
        Useful for debugging. Default is False.
    temp_dir : str, optional
        Create temporary files in this directory.
    multifile_merge : bool, optional
        Store temporary merge chunks as separate .cool files rather than in 
        a single multi-cooler file. Default is False.
    metadata : dict, optional
        Experiment metadata to store in the file. Must be JSON compatible.
    assembly : str, optional
        Name of genome assembly.
    h5opts : dict, optional
        HDF5 dataset filter options to use (compression, shuffling,
        checksumming, etc.). Default is to use autochunking and GZIP
        compression, level 6.
    append : bool, optional
        Append new Cooler to the file if it exists. If False, an existing file
        with the same name will be truncated. Default is False.
    lock : multiprocessing.Lock, optional
        Optional lock to control concurrent access to the output file.

    See also
    --------
    sanitize_records
    sanitize_pixels

    """
    from .api import Cooler
    chromsizes = get_chromsizes(bins)
    bins = bins.copy()
    bins['chrom'] = bins['chrom'].astype(object)
    
    temp_files = []

    if multifile_merge:
        for i, chunk in enumerate(chunks):
            tf = tempfile.NamedTemporaryFile(
                suffix='.cool', 
                delete=delete_temp,
                dir=temp_dir)
            temp_files.append(tf)
            logger.info('Writing chunk {}: {}'.format(i, tf.name))
            create(tf.name, bins, chunk, columns=columns, dtype=dtype)
        chunks = CoolerMerger([Cooler(tf.name) for tf in temp_files], mergebuf)
    else:
        tf = tempfile.NamedTemporaryFile(
                suffix='.multi.cool', 
                delete=delete_temp,
                dir=temp_dir)
        temp_files.append(tf)
        
        uris = []
        for i, chunk in enumerate(chunks):
            uri = tf.name + '::' + str(i)
            uris.append(uri)
            logger.info('Writing chunk {}: {}'.format(i, uri))
            create(uri, bins, chunk, columns=columns, dtype=dtype)
        chunks = CoolerMerger([Cooler(uri) for uri in uris], mergebuf)

    logger.info('Merging into {}'.format(cool_uri))
    create(cool_uri, bins, chunks, columns=columns, dtype=dtype, **kwargs)

    del temp_files


# Exports
from ._binning import (ContactBinner, HDF5Aggregator, TabixAggregator,
                       PairixAggregator, CoolerAggregator, CoolerMerger,
                       SparseLoader, BedGraph2DLoader, ArrayLoader,
                       sanitize_pixels, sanitize_records)

from ._writer import (write_chroms, write_bins, prepare_pixels, write_pixels, 
                      write_indexes, write_info, index_bins, index_pixels, 
                      MAGIC, URL, PIXEL_FIELDS, PIXEL_DTYPES)
