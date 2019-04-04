# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
from datetime import datetime
from six.moves import map
from pandas.api.types import is_categorical, is_integer
import os.path as op
import pandas as pd
import numpy as np
import posixpath
import tempfile
import warnings
import h5py
import json
import six

from .._version import __version__, __format_version__
from .._logging import get_logger
from ..core import put, get
from ..util import (
    parse_cooler_uri,
    get_chromsizes,
    get_binsize,
    infer_meta,
    get_meta,
    rlencode,
)
from ._ingest import validate_pixels
from . import (
    MAGIC,
    URL,
    CHROM_DTYPE,
    CHROMID_DTYPE,
    CHROMSIZE_DTYPE,
    COORD_DTYPE,
    BIN_DTYPE,
    COUNT_DTYPE,
    CHROMOFFSET_DTYPE,
    BIN1OFFSET_DTYPE,
    PIXEL_FIELDS,
    PIXEL_DTYPES,
)

logger = get_logger('cooler.create')


def write_chroms(grp, chroms, h5opts):
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
    names = np.array(chroms['name'], dtype=CHROM_DTYPE)  # auto-adjusts char length
    grp.create_dataset('name',
                       shape=(n_chroms,),
                       dtype=names.dtype,
                       data=names,
                       **h5opts)
    grp.create_dataset('length',
                       shape=(n_chroms,),
                       dtype=CHROMSIZE_DTYPE,
                       data=chroms['length'],
                       **h5opts)

    # Extra columns
    columns = list(chroms.keys())
    for col in ['name', 'length']:
        columns.remove(col)
    if columns:
        put(grp, chroms[columns])


def write_bins(grp, bins, chromnames, h5opts, chrom_as_enum=True):
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
    chrom_ids = [idmap[chrom] for chrom in bins['chrom']]
    if chrom_as_enum:
        chrom_dtype = h5py.special_dtype(enum=(CHROMID_DTYPE, idmap))
    else:
        chrom_dtype = CHROMID_DTYPE

    # Store bins
    try:
        chrom_dset = grp.create_dataset('chrom',
                           shape=(n_bins,),
                           dtype=chrom_dtype,
                           data=chrom_ids,
                           **h5opts)
    except ValueError:
        # If too many scaffolds for HDF5 enum header,
        # try storing chrom IDs as raw int instead
        if chrom_as_enum:
            chrom_as_enum = False
            chrom_dtype = CHROMID_DTYPE
            chrom_dset = grp.create_dataset('chrom',
                               shape=(n_bins,),
                               dtype=chrom_dtype,
                               data=chrom_ids,
                               **h5opts)
        else:
            raise
    if not chrom_as_enum:
        chrom_dset.attrs['enum_path'] = u'/chroms/name'

    grp.create_dataset('start',
                       shape=(n_bins,),
                       dtype=COORD_DTYPE,
                       data=bins['start'],
                       **h5opts)
    grp.create_dataset('end',
                       shape=(n_bins,),
                       dtype=COORD_DTYPE,
                       data=bins['end'],
                       **h5opts)

    # Extra columns
    columns = list(bins.keys())
    for col in ['chrom', 'start', 'end']:
        columns.remove(col)
    if columns:
        put(grp, bins[columns])


def prepare_pixels(grp, n_bins, max_size, columns, dtypes, h5opts):
    columns = list(columns)
    init_size = min(5 * n_bins, max_size)
    grp.create_dataset('bin1_id',
                       dtype=dtypes.get('bin1_id', BIN_DTYPE),
                       shape=(init_size,),
                       maxshape=(max_size,),
                       **h5opts)
    grp.create_dataset('bin2_id',
                       dtype=dtypes.get('bin2_id', BIN_DTYPE),
                       shape=(init_size,),
                       maxshape=(max_size,),
                       **h5opts)

    if 'count'in columns:
        grp.create_dataset('count',
                           dtype=dtypes.get('count', COUNT_DTYPE),
                           shape=(init_size,),
                           maxshape=(max_size,),
                           **h5opts)

    for col in ['bin1_id', 'bin2_id', 'count']:
        try:
            columns.remove(col)
        except ValueError:
            pass

    if columns:
        for col in columns:
            grp.create_dataset(col,
                               dtype=dtypes.get(col, float),
                               shape=(init_size,),
                               maxshape=(max_size,),
                               **h5opts)


def write_pixels(filepath, grouppath, columns, iterable, h5opts, lock):
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
            chunk = {k: v.values for k, v in six.iteritems(chunk)}

        try:
            if lock is not None:
                lock.acquire()

            logger.debug("writing chunk {}".format(i))

            with h5py.File(filepath, 'r+') as fw:
                grp = fw[grouppath]
                dsets = [grp[col] for col in columns]

                n = len(chunk[columns[0]])
                for col, dset in zip(columns, dsets):
                    dset.resize((nnz + n,))
                    dset[nnz:nnz+n] = chunk[col]
                nnz += n
                if 'count' in chunk:
                    total += chunk['count'].sum()

                fw.flush()

        finally:
            if lock is not None:
                lock.release()

    return nnz, total


def index_pixels(grp, n_bins, nnz):
    bin1 = grp['bin1_id']
    bin1_offset = np.zeros(n_bins + 1, dtype=BIN1OFFSET_DTYPE)
    curr_val = 0
    for start, length, value in zip(*rlencode(bin1, 1000000)):
        bin1_offset[curr_val:value + 1] = start
        curr_val = value + 1
    bin1_offset[curr_val:] = nnz
    return bin1_offset


def index_bins(grp, n_chroms, n_bins):
    chrom_ids = grp['chrom']
    chrom_offset = np.zeros(n_chroms + 1, dtype=CHROMOFFSET_DTYPE)
    curr_val = 0
    for start, length, value in zip(*rlencode(chrom_ids)):
        chrom_offset[curr_val:value + 1] = start
        curr_val = value + 1
    chrom_offset[curr_val:] = n_bins
    return chrom_offset


def write_indexes(grp, chrom_offset, bin1_offset, h5opts):
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
        **h5opts)
    grp.create_dataset(
        "bin1_offset",
        shape=(len(bin1_offset),),
        dtype=BIN1OFFSET_DTYPE,
        data=bin1_offset,
        **h5opts)


def write_info(grp, info):
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
    assert 'nbins' in info
    assert 'nnz' in info
    info.setdefault('genome-assembly', 'unknown')
    info['metadata'] = json.dumps(info.get('metadata', {}))
    info['creation-date'] = datetime.now().isoformat()
    info['generated-by'] = six.text_type('cooler-' + __version__)
    info['format'] = MAGIC
    info['format-version'] = six.text_type(__format_version__)
    info['format-url'] = URL
    grp.attrs.update(info)


def _rename_chroms(grp, rename_dict, h5opts):
    chroms = get(grp['chroms']).set_index('name')
    n_chroms = len(chroms)
    new_names = np.array(chroms.rename(rename_dict).index.values,
                         dtype=CHROM_DTYPE)  # auto-adjusts char length

    del grp['chroms/name']
    grp['chroms'].create_dataset('name',
                       shape=(n_chroms,),
                       dtype=new_names.dtype,
                       data=new_names,
                       **h5opts)

    bins = get(grp['bins'])
    n_bins = len(bins)
    idmap = dict(zip(new_names, range(n_chroms)))
    if is_categorical(bins['chrom']) or is_integer(bins['chrom']):
        chrom_ids = bins['chrom'].cat.codes
        chrom_dtype = h5py.special_dtype(enum=(CHROMID_DTYPE, idmap))
        del grp['bins/chrom']
        try:
            chrom_dset = grp['bins'].create_dataset('chrom',
                               shape=(n_bins,),
                               dtype=chrom_dtype,
                               data=chrom_ids,
                               **h5opts)
        except ValueError:
            # If HDF5 enum header would be too large,
            # try storing chrom IDs as raw int instead
            chrom_dtype = CHROMID_DTYPE
            chrom_dset = grp['bins'].create_dataset('chrom',
                               shape=(n_bins,),
                               dtype=chrom_dtype,
                               data=chrom_ids,
                               **h5opts)


def rename_chroms(clr, rename_dict, h5opts=None):
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

    with clr.open('r+') as f:
        _rename_chroms(f, rename_dict, h5opts)
    clr._refresh()


def _get_dtypes_arg(dtypes, kwargs):
    if 'dtype' in kwargs:
        if dtypes is None:
            dtypes = kwargs.pop('dtype')
            warnings.warn("Use dtypes= instead of dtype=", FutureWarning)
        else:
            raise ValueError(
                'Received both "dtypes" and "dtype" arguments. '
                'Please use "dtypes" to provide a column name -> dtype mapping. '
                '"dtype" remains as an alias but is deprecated.')
    return dtypes


def _set_h5opts(h5opts):
    result = {}
    if h5opts is not None:
        result.update(h5opts)
    available_opts = {'chunks', 'maxshape', 'compression', 'compression_opts',
        'scaleoffset', 'shuffle', 'fletcher32' , 'fillvalue', 'track_times'}
    for key in result.keys():
        if key not in available_opts:
            raise ValueError("Unknown storage option '{}'.".format(key))
    result.setdefault('compression', 'gzip')
    if result['compression'] == 'gzip' and 'compression_opts' not in result:
        result['compression_opts'] = 6
    result.setdefault('shuffle', True)
    return result


def create(cool_uri, bins, pixels, columns=None, dtypes=None, metadata=None,
           assembly=None, symmetric_upper=True, mode=None, h5opts=None,
           boundscheck=True, triucheck=True, dupcheck=True,
           ensure_sorted=False, lock=None, append=False, **kwargs):
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
    file_path, group_path = parse_cooler_uri(cool_uri)

    if mode is None:
        mode = 'a' if append else 'w'

    h5opts = _set_h5opts(h5opts)

    if not isinstance(bins, pd.DataFrame):
        raise ValueError(
            "Second positional argument must be a pandas DataFrame. "
            "Note that the `chromsizes` argument is now deprecated: "
            "see documentation for `create`.")

    dtypes = _get_dtypes_arg(dtypes, kwargs)

    for col in ['chrom', 'start', 'end']:
        if col not in bins.columns:
            raise ValueError("Missing column from bin table: '{}'.".format(col))

    # Populate expected pixel column names. Include user-provided value columns.
    if columns is None:
        columns = ['bin1_id', 'bin2_id', 'count']
    else:
        columns = list(columns)
        for col in ['bin1_id', 'bin2_id']:  # don't include count!
            if col not in columns:
                columns.insert(0, col)

    # Populate dtypes for expected pixel columns, and apply user overrides.
    if dtypes is None:
        dtypes = dict(PIXEL_DTYPES)
    else:
        dtypes_ = dict(dtypes)
        dtypes = dict(PIXEL_DTYPES)
        dtypes.update(dtypes_)

    # Get empty "meta" header frame (assigns the undeclared dtypes).
    # Any columns from the input not in meta will be ignored.
    meta = get_meta(columns, dtypes, default_dtype=float)

    # Determine the appropriate iterable
    try:
        from dask.dataframe import DataFrame as dask_df
    except (ImportError, AttributeError):
        dask_df = ()

    if isinstance(pixels, dask_df):
        iterable = map(lambda x: x.compute(), pixels.to_delayed())
        input_columns = infer_meta(pixels).columns
    elif isinstance(pixels, pd.DataFrame):
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
    chroms = pd.DataFrame(
        {'name': chromnames, 'length': lengths},
        columns=['name', 'length'])
    binsize = get_binsize(bins)
    n_chroms = len(chroms)
    n_bins = len(bins)

    if not symmetric_upper and triucheck:
        warnings.warn(
            "Creating a non-symmetric matrix, but `triucheck` was set to True. "
            "Changing to False.")
        triucheck = False

    # Chain input validation to the end of the pipeline
    if boundscheck or triucheck or dupcheck or ensure_sorted:
        validator = validate_pixels(
            n_bins, boundscheck, triucheck, dupcheck, ensure_sorted)
        iterable = map(validator, iterable)

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
        if symmetric_upper:
            max_size = n_bins * (n_bins - 1) // 2 + n_bins
        else:
            max_size = n_bins * n_bins
        prepare_pixels(grp, n_bins, max_size, meta.columns, dict(meta.dtypes), h5opts)

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
        info['bin-type'] = u"fixed" if binsize is not None else u"variable"
        info['bin-size'] = binsize if binsize is not None else u"null"
        info['storage-mode'] = u"symmetric-upper" if symmetric_upper else u"square"
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


def create_from_unordered(cool_uri, bins, chunks, columns=None, dtypes=None,
                          mode=None, mergebuf=int(20e6), delete_temp=True,
                          temp_dir=None, max_merge=200, **kwargs):
    """
    Create a Cooler in two passes via an external sort mechanism. In the first
    pass, a sequence of data chunks are processed and sorted in memory and saved
    to temporary Coolers. In the second pass, the temporary Coolers are merged
    into the output. This way the individual chunks do not need to be provided
    in any particular order.

    """
    from ..api import Cooler
    from ..reduce import CoolerMerger
    chromsizes = get_chromsizes(bins)
    bins = bins.copy()
    bins['chrom'] = bins['chrom'].astype(object)

    if columns is not None:
        columns = [col for col in columns if col not in {'bin1_id', 'bin2_id'}]

    if temp_dir is None:
        temp_dir = op.dirname(parse_cooler_uri(cool_uri)[0])
    elif temp_dir == '-':
        temp_dir = None  # makes tempfile module use the system dir

    dtypes = _get_dtypes_arg(dtypes, kwargs)

    temp_files = []

    # Sort pass
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
        create(uri, bins, chunk,
               columns=columns, dtypes=dtypes, mode='a', **kwargs)

    # Merge passes
    n = len(uris)
    if n > max_merge > 0:
        # Recursive merge: do the first of two merge passes.
        # Divide into ~sqrt(n) merges
        edges = np.linspace(0, n, int(np.sqrt(n)), dtype=int)

        tf2 = tempfile.NamedTemporaryFile(
            suffix='.multi.cool',
            delete=delete_temp,
            dir=temp_dir)
        temp_files.append(tf2)
        uris2 = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            chunk_subset = CoolerMerger(
                [Cooler(uri) for uri in uris[lo:hi]],
                mergebuf,
                columns=columns)
            uri = tf2.name + '::' + "{}-{}".format(lo, hi)
            uris2.append(uri)
            logger.info('Merging chunks {}-{}: {}'.format(lo, hi, uri))
            create(uri, bins, chunk_subset,
                   columns=columns, dtypes=dtypes, mode='a', **kwargs)

        final_uris = uris2
    else:
        # Do a single merge pass
        final_uris = uris

    # Do the final merge pass
    chunks = CoolerMerger(
        [Cooler(uri) for uri in final_uris],
        mergebuf,
        columns=columns)
    logger.info('Merging into {}'.format(cool_uri))
    create(cool_uri, bins, chunks, columns=columns, dtypes=dtypes, mode=mode,
           **kwargs)

    del temp_files


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
    h5opts = _set_h5opts(h5opts)

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


def create_cooler(cool_uri, bins, pixels, columns=None, dtypes=None,
                  metadata=None, assembly=None, ordered=False,
                  symmetric_upper=True, mode=None, mergebuf=int(20e6),
                  delete_temp=True, temp_dir=None, max_merge=200,
                  boundscheck=True, dupcheck=True, triucheck=True,
                  ensure_sorted=False, h5opts=None, lock=None):
    """
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
        Useful for debugging. Default is False.
    temp_dir : str, optional
        Create temporary files in a specified directory instead of the same
        directory as the output file. Pass ``-`` to use the system default.
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

    See also
    --------
    cooler.create.sanitize_records
    cooler.create.sanitize_pixels

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

    """
    # dispatch to the approprate creation method
    if isinstance(pixels, (pd.DataFrame, dict)):
        pixels = pd.DataFrame(pixels).sort_values(['bin1_id', 'bin2_id'])
        ordered = True

    if ordered:
        create(
            cool_uri, bins, pixels,
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
            lock=lock)
    else:
        create_from_unordered(
            cool_uri, bins, pixels,
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
            max_merge=max_merge)
