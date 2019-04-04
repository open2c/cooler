# -*- coding: utf-8 -*-
from __future__ import division, print_function
from collections import OrderedDict, defaultdict
from contextlib import contextmanager
import six
import re
import os

from pandas.api.types import is_scalar, is_integer
import numpy as np
import pandas as pd
import h5py


def partition(start, stop, step):
    """Partition an integer interval into equally-sized subintervals.
    Like builtin :py:func:`range`, but yields pairs of end points.

    Examples
    --------
    >>> for lo, hi in partition(0, 9, 2):
           print(lo, hi)
    0 2
    2 4
    4 6
    6 8
    8 9

    """
    return ((i, min(i+step, stop))
                for i in range(start, stop, step))


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


def atoi(s):
    return int(s.replace(',', ''))


def parse_humanized(s):
    _NUMERIC_RE = re.compile('([0-9,.]+)')
    _, value, unit = _NUMERIC_RE.split(s.replace(',', ''))
    if not len(unit):
        return int(value)

    value = float(value)
    unit = unit.upper().strip()
    if unit in ('K', 'KB'):
        value *= 1000
    elif unit in ('M', 'MB'):
        value *= 1000000
    elif unit in ('G', 'GB'):
        value *= 1000000000
    else:
        raise ValueError("Unknown unit '{}'".format(unit))
    return int(value)


def parse_region_string(s):
    """
    Parse a UCSC-style genomic region string into a triple.

    Parameters
    ----------
    s : str
        UCSC-style string, e.g. "chr5:10,100,000-30,000,000". Ensembl and FASTA
        style sequence names are allowed. End coordinate must be greater than or
        equal to start.

    Returns
    -------
    (str, int or None, int or None)

    """
    def _tokenize(s):
        token_spec = [
            ('HYPHEN', r'-'),
            ('COORD',  r'[0-9,]+(\.[0-9]*)?(?:[a-z]+)?'),
            ('OTHER',  r'.+')
        ]
        tok_regex = r'\s*' + r'|\s*'.join(
            r'(?P<%s>%s)' % pair for pair in token_spec)
        tok_regex = re.compile(tok_regex, re.IGNORECASE)
        for match in tok_regex.finditer(s):
            typ = match.lastgroup
            yield typ, match.group(typ)


    def _check_token(typ, token, expected):
        if typ is None:
            raise ValueError('Expected {} token missing'.format(' or '.join(expected)))
        else:
            if typ not in expected:
                raise ValueError('Unexpected token "{}"'.format(token))


    def _expect(tokens):
        typ, token = next(tokens, (None, None))
        _check_token(typ, token, ['COORD'])
        start = parse_humanized(token)

        typ, token = next(tokens, (None, None))
        _check_token(typ, token, ['HYPHEN'])

        typ, token = next(tokens, (None, None))
        if typ is None:
            return start, None

        _check_token(typ, token, ['COORD'])
        end = parse_humanized(token)
        if end < start:
            raise ValueError('End coordinate less than start')

        return start, end

    parts = s.split(':')
    chrom = parts[0].strip()
    if not len(chrom):
        raise ValueError("Chromosome name cannot be empty")
    if len(parts) < 2:
        return (chrom, None, None)
    start, end = _expect(_tokenize(parts[1]))
    return (chrom, start, end)


def parse_region(reg, chromsizes=None):
    """
    Genomic regions are represented as half-open intervals (0-based starts,
    1-based ends) along the length coordinate of a contig/scaffold/chromosome.

    Parameters
    ----------
    reg : str or tuple
        UCSC-style genomic region string, or
        Triple (chrom, start, end), where ``start`` or ``end`` may be ``None``.
    chromsizes : mapping, optional
        Lookup table of scaffold lengths to check against ``chrom`` and the
        ``end`` coordinate. Required if ``end`` is not supplied.

    Returns
    -------
    A well-formed genomic region triple (str, int, int)

    """
    if isinstance(reg, six.string_types):
        chrom, start, end = parse_region_string(reg)
    else:
        chrom, start, end = reg
        start = int(start) if start is not None else start
        end = int(end) if end is not None else end

    try:
        clen = chromsizes[chrom] if chromsizes is not None else None
    except KeyError:
        raise ValueError("Unknown sequence label: {}".format(chrom))

    start = 0 if start is None else start
    if end is None:
        if clen is None:  # TODO --- remove?
            raise ValueError("Cannot determine end coordinate.")
        end = clen

    if end < start:
        raise ValueError("End cannot be less than start")

    if start < 0 or (clen is not None and end > clen):
        raise ValueError(
            "Genomic region out of bounds: [{}, {})".format(start, end))

    return chrom, start, end


def natsort_key(s, _NS_REGEX=re.compile(r'(\d+)', re.U)):
    return tuple([int(x) if x.isdigit() else x
                 for x in _NS_REGEX.split(s) if x])


def natsorted(iterable):
    return sorted(iterable, key=natsort_key)


def argnatsort(array):
    array = np.asarray(array)
    if not len(array): return np.array([], dtype=int)
    cols = tuple(zip(*(natsort_key(x) for x in array)))
    return np.lexsort(cols[::-1])


def read_chromsizes(filepath_or,
                   name_patterns=(r'^chr[0-9]+$', r'^chr[XY]$', r'^chrM$'),
                   all_names=False,
                   **kwargs):
    """
    Parse a ``<db>.chrom.sizes`` or ``<db>.chromInfo.txt`` file from the UCSC
    database, where ``db`` is a genome assembly name.

    Parameters
    ----------
    filepath_or : str or file-like
        Path or url to text file, or buffer.
    name_patterns : sequence, optional
        Sequence of regular expressions to capture desired sequence names.
        Each corresponding set of records will be sorted in natural order.
    all_names : bool, optional
        Whether to return all contigs listed in the file. Default is
        ``False``.

    Returns
    -------
    :py:class:`pandas.Series`
        Series of integer bp lengths indexed by sequence name.

    References
    ----------
    * `UCSC assembly terminology <http://genome.ucsc.edu/FAQ/FAQdownloads.html#download9>`_
    * `GRC assembly terminology <https://www.ncbi.nlm.nih.gov/grc/help/definitions>`_

    """
    if isinstance(filepath_or, six.string_types) and filepath_or.endswith('.gz'):
        kwargs.setdefault('compression', 'gzip')
    chromtable = pd.read_csv(
        filepath_or, sep='\t', usecols=[0, 1],
        names=['name', 'length'], dtype={'name':str}, **kwargs)
    if not all_names:
        parts = []
        for pattern in name_patterns:
            part = chromtable[chromtable['name'].str.contains(pattern)]
            part = part.iloc[argnatsort(part['name'])]
            parts.append(part)
        chromtable = pd.concat(parts, axis=0)
    chromtable.index = chromtable['name'].values
    return chromtable['length']


def fetch_chromsizes(db, **kwargs):
    """
    Download chromosome sizes from UCSC as a :py:class:`pandas.Series`, indexed
    by chromosome label.

    """
    return read_chromsizes(
        'http://hgdownload.cse.ucsc.edu/goldenPath/{}/database/chromInfo.txt.gz'.format(db),
        **kwargs)


def load_fasta(names, *filepaths):
    """
    Load lazy FASTA records from one or multiple files without reading them into
    memory.

    Parameters
    ----------
    names : sequence of str
        Names of sequence records in FASTA file or files.
    filepaths : str
        Paths to one or more FASTA files to gather records from.

    Returns
    -------
    OrderedDict of sequence name -> sequence record

    """
    import pyfaidx
    if len(filepaths) == 0:
        raise ValueError("Need at least one file")

    if len(filepaths) == 1:
        fa = pyfaidx.Fasta(filepaths[0], as_raw=True)

    else:
        fa = {}
        for filepath in filepaths:
            fa.update(pyfaidx.Fasta(filepath, as_raw=True).records)

    records = OrderedDict((chrom, fa[chrom]) for chrom in names)
    return records


def binnify(chromsizes, binsize):
    """
    Divide a genome into evenly sized bins.

    Parameters
    ----------
    chromsizes : Series
        pandas Series indexed by chromosome name with chromosome lengths in bp.
    binsize : int
        size of bins in bp

    Returns
    -------
    bins : :py:class:`pandas.DataFrame`
        Dataframe with columns: ``chrom``, ``start``, ``end``.

    """
    def _each(chrom):
        clen = chromsizes[chrom]
        n_bins = int(np.ceil(clen / binsize))
        binedges = np.arange(0, (n_bins+1)) * binsize
        binedges[-1] = clen
        return pd.DataFrame({
                'chrom': [chrom]*n_bins,
                'start': binedges[:-1],
                'end': binedges[1:],
            }, columns=['chrom', 'start', 'end'])

    bintable = pd.concat(
        map(_each, chromsizes.keys()),
        axis=0,
        ignore_index=True)

    bintable['chrom'] = pd.Categorical(
        bintable['chrom'],
        categories=list(chromsizes.index),
        ordered=True)

    return bintable

make_bintable = binnify


def digest(fasta_records, enzyme):
    """
    Divide a genome into restriction fragments.

    Parameters
    ----------
    fasta_records : OrderedDict
        Dictionary of chromosome names to sequence records.
    enzyme: str
        Name of restriction enzyme (e.g., 'DpnII').

    Returns
    -------
    frags : :py:class:`pandas.DataFrame`
        Dataframe with columns: ``chrom``, ``start``, ``end``.

    """
    import Bio.Restriction as biorst
    import Bio.Seq as bioseq
    # http://biopython.org/DIST/docs/cookbook/Restriction.html#mozTocId447698
    chroms = fasta_records.keys()
    try:
        cut_finder = getattr(biorst, enzyme).search
    except AttributeError:
        raise ValueError('Unknown enzyme name: {}'.format(enzyme))

    def _each(chrom):
        seq = bioseq.Seq(str(fasta_records[chrom]))
        cuts = np.r_[0, np.array(cut_finder(seq)) + 1, len(seq)].astype(int)
        n_frags = len(cuts) - 1

        frags = pd.DataFrame({
            'chrom': [chrom] * n_frags,
            'start': cuts[:-1],
            'end': cuts[1:]},
            columns=['chrom', 'start', 'end'])
        return frags

    return pd.concat(map(_each, chroms), axis=0, ignore_index=True)


def get_binsize(bins):
    """
    Infer bin size from a bin DataFrame. Assumes that the last bin of each
    contig is allowed to differ in size from the rest.

    Returns
    -------
    int or None if bins are non-uniform

    """
    sizes = set()
    for chrom, group in bins.groupby('chrom'):
        sizes.update((group['end'] - group['start']).iloc[:-1].unique())
        if len(sizes) > 1:
            return None
    if len(sizes) == 1:
        return next(iter(sizes))
    else:
        return None


def get_chromsizes(bins):
    """
    Infer chromsizes Series from a bin DataFrame. Assumes that the last bin of
    each contig is allowed to differ in size from the rest.

    Returns
    -------
    int or None if bins are non-uniform

    """
    chromtable = (
        bins.drop_duplicates(['chrom'], keep='last')[['chrom', 'end']]
            .reset_index(drop=True)
            .rename(columns={'chrom': 'name', 'end': 'length'})
    )
    chroms, lengths = list(chromtable['name']), list(chromtable['length'])
    return pd.Series(index=chroms, data=lengths)


def bedslice(grouped, chromsizes, region):
    """
    Range query on a BED-like dataframe with non-overlapping intervals.

    """
    chrom, start, end = parse_region(region, chromsizes)
    result = grouped.get_group(chrom)
    if start > 0 or end < chromsizes[chrom]:
        lo = result['end'].values.searchsorted(start, side='right')
        hi = lo + result['start'].values[lo:].searchsorted(end, side='left')
        result = result.iloc[lo:hi]
    return result


def lexbisect(arrays, values, side='left', lo=0, hi=None):
    """
    Bisection search on lexically sorted arrays.

    Parameters
    ----------
    arrays : sequence of k 1-D array-like
        Each "array" can be any sequence that supports scalar integer indexing.
        The arrays are assumed to have the same length and values lexsorted
        from left to right.
    values : sequence of k values
        Values that would be inserted into the arrays.
    side : {'left', 'right'}, optional
        If ‘left’, the index of the first suitable location found is given.
        If ‘right’, return the last such index. If there is no suitable index,
        return either 0 or N (where N is the length of each array).
    lo, hi : int, optional
        Bound the slice to be searched (default 'lo' is 0 and 'hi' is N).

    Returns
    -------
    i : int
        Insertion index.

    Examples
    --------
    >>> h5 = h5py.File('mytable.h5', 'r')  # doctest: +SKIP
    >>> lexbisect([h5['chrom'], h5['start']], [1, 100000], side='right')  # doctest: +SKIP
    2151688

    """
    arrays = tuple(arrays)
    values = tuple(values)

    if lo < 0:
        raise ValueError('lo must be non-negative')
    if hi is None:
        hi = len(arrays[0])

    if side == 'left':
        while lo < hi:
            mid = (lo + hi) // 2
            if tuple(arr[mid] for arr in arrays) < values:
                lo = mid + 1
            else:
                hi = mid
    elif side == 'right':
        while lo < hi:
            mid = (lo + hi) // 2
            if values < tuple(arr[mid] for arr in arrays):
                hi = mid
            else:
                lo = mid + 1
    else:
        raise ValueError("side must be 'left' or 'right'")

    return lo


def asarray_or_dataset(x):
    return x if isinstance(x, h5py.Dataset) else np.asarray(x)


def rlencode(array, chunksize=None):
    """
    Run length encoding.
    Based on http://stackoverflow.com/a/32681075, which is based on the rle
    function from R.

    Parameters
    ----------
    x : 1D array_like
        Input array to encode
    dropna: bool, optional
        Drop all runs of NaNs.

    Returns
    -------
    start positions, run lengths, run values

    """
    where = np.flatnonzero
    array = asarray_or_dataset(array)
    n = len(array)
    if n == 0:
        return (np.array([], dtype=int),
                np.array([], dtype=int),
                np.array([], dtype=array.dtype))

    if chunksize is None:
        chunksize = n

    starts, values = [], []
    last_val = np.nan
    for i in range(0, n, chunksize):
        x = array[i:i+chunksize]
        locs = where(x[1:] != x[:-1]) + 1
        if x[0] != last_val:
            locs = np.r_[0, locs]
        starts.append(i + locs)
        values.append(x[locs])
        last_val = x[-1]
    starts = np.concatenate(starts)
    lengths = np.diff(np.r_[starts, n])
    values = np.concatenate(values)

    return starts, lengths, values


def cmd_exists(cmd):
    return any(os.access(os.path.join(path, cmd), os.X_OK)
                for path in os.environ['PATH'].split(os.pathsep))


def mad(data, axis=None):
    return np.median(np.abs(data - np.median(data, axis)), axis)


@contextmanager
def open_hdf5(fp, mode='r', *args, **kwargs):
    """
    Context manager like ``h5py.File`` but accepts already open HDF5 file
    handles which do not get closed on teardown.

    Parameters
    ----------
    fp : str or ``h5py.File`` object
        If an open file object is provided, it passes through unchanged,
        provided that the requested mode is compatible.
        If a filepath is passed, the context manager will close the file on
        tear down.

    mode : str
        * r        Readonly, file must exist
        * r+       Read/write, file must exist
        * a        Read/write if exists, create otherwise
        * w        Truncate if exists, create otherwise
        * w- or x  Fail if exists, create otherwise

    """
    if isinstance(fp, six.string_types):
        own_fh = True
        fh = h5py.File(fp, mode, *args, **kwargs)
    else:
        own_fh = False
        if mode == 'r' and fp.file.mode == 'r+':
            #warnings.warn("File object provided is writeable but intent is read-only")
            pass
        elif mode in ('r+', 'a') and fp.file.mode == 'r':
            raise ValueError("File object provided is not writeable")
        elif mode == 'w':
            raise ValueError("Cannot truncate open file")
        elif mode in ('w-', 'x'):
            raise ValueError("File exists")
        fh = fp
    try:
        yield fh
    finally:
        if own_fh:
            fh.close()


class closing_hdf5(h5py.Group):
    def __init__(self, grp):
        super(closing_hdf5, self).__init__(grp.id)
    def __enter__(self):
        return self
    def __exit__(self, *exc_info):
        return self.file.close()
    def close(self):
        self.file.close()


def attrs_to_jsonable(attrs):
    out = dict(attrs)
    for k, v in attrs.items():
        try:
            out[k] = np.asscalar(v)
        except ValueError:
            out[k] = v.tolist()
        except AttributeError:
            out[k] = v
    return out


def unstar(func):
    def unstarred(args):
        return func(*args)
    return unstarred


def infer_meta(x, index=None):
    """
    Extracted and modified from dask/dataframe/utils.py :
        make_meta (BSD licensed)

    Create an empty pandas object containing the desired metadata.

    Parameters
    ----------
    x : dict, tuple, list, pd.Series, pd.DataFrame, pd.Index, dtype, scalar
        To create a DataFrame, provide a `dict` mapping of `{name: dtype}`, or
        an iterable of `(name, dtype)` tuples. To create a `Series`, provide a
        tuple of `(name, dtype)`. If a pandas object, names, dtypes, and index
        should match the desired output. If a dtype or scalar, a scalar of the
        same dtype is returned.
    index :  pd.Index, optional
        Any pandas index to use in the metadata. If none provided, a
        `RangeIndex` will be used.

    Examples
    --------
    >>> make_meta([('a', 'i8'), ('b', 'O')])
    Empty DataFrame
    Columns: [a, b]
    Index: []
    >>> make_meta(('a', 'f8'))
    Series([], Name: a, dtype: float64)
    >>> make_meta('i8')
    1

    """

    UNKNOWN_CATEGORIES = '__UNKNOWN_CATEGORIES__'

    def _scalar_from_dtype(dtype):
        if dtype.kind in ('i', 'f', 'u'):
            return dtype.type(1)
        elif dtype.kind == 'c':
            return dtype.type(complex(1, 0))
        elif dtype.kind in _simple_fake_mapping:
            o = _simple_fake_mapping[dtype.kind]
            return o.astype(dtype) if dtype.kind in ('m', 'M') else o
        else:
            raise TypeError("Can't handle dtype: {0}".format(dtype))

    def _nonempty_scalar(x):
        if isinstance(x, (pd.Timestamp, pd.Timedelta, pd.Period)):
            return x
        elif np.isscalar(x):
            dtype = x.dtype if hasattr(x, 'dtype') else np.dtype(type(x))
            return _scalar_from_dtype(dtype)
        else:
            raise TypeError("Can't handle meta of type "
                            "'{0}'".format(type(x).__name__))

    def _empty_series(name, dtype, index=None):
        if isinstance(dtype, str) and dtype == 'category':
            return pd.Series(pd.Categorical([UNKNOWN_CATEGORIES]),
                             name=name, index=index).iloc[:0]
        return pd.Series([], dtype=dtype, name=name, index=index)

    if hasattr(x, '_meta'):
        return x._meta
    if isinstance(x, (pd.Series, pd.DataFrame)):
        return x.iloc[0:0]
    elif isinstance(x, pd.Index):
        return x[0:0]
    index = index if index is None else index[0:0]

    if isinstance(x, dict):
        return pd.DataFrame({c: _empty_series(c, d, index=index)
                             for (c, d) in x.items()}, index=index)
    if isinstance(x, tuple) and len(x) == 2:
        return _empty_series(x[0], x[1], index=index)
    elif isinstance(x, (list, tuple)):
        if not all(isinstance(i, tuple) and len(i) == 2 for i in x):
            raise ValueError("Expected iterable of tuples of (name, dtype), "
                             "got {0}".format(x))
        return pd.DataFrame({c: _empty_series(c, d, index=index) for (c, d) in x},
                            columns=[c for c, d in x], index=index)
    elif not hasattr(x, 'dtype') and x is not None:
        # could be a string, a dtype object, or a python type. Skip `None`,
        # because it is implictly converted to `dtype('f8')`, which we don't
        # want here.
        try:
            dtype = np.dtype(x)
            return _scalar_from_dtype(dtype)
        except:
            # Continue on to next check
            pass

    if is_scalar(x):
        return _nonempty_scalar(x)

    raise TypeError("Don't know how to create metadata from {0}".format(x))


def get_meta(columns, dtype=None, index_columns=None, index_names=None,
             default_dtype=np.object):
    """
    Extracted and modified from pandas/io/parsers.py :
        _get_empty_meta (BSD licensed).

    """
    columns = list(columns)

    # Convert `dtype` to a defaultdict of some kind.
    # This will enable us to write `dtype[col_name]`
    # without worrying about KeyError issues later on.
    if not isinstance(dtype, dict):
        # if dtype == None, default will be default_dtype.
        dtype = defaultdict(lambda: dtype or default_dtype)
    else:
        # Save a copy of the dictionary.
        _dtype = dtype.copy()
        dtype = defaultdict(lambda: default_dtype)

        # Convert column indexes to column names.
        for k, v in six.iteritems(_dtype):
            col = columns[k] if is_integer(k) else k
            dtype[col] = v

    if index_columns is None or index_columns is False:
        index = pd.Index([])
    else:
        data = [pd.Series([], dtype=dtype[name]) for name in index_names]
        if len(data) == 1:
            index = pd.Index(data[0], name=index_names[0])
        else:
            index = pd.MultiIndex.from_arrays(data, names=index_names)
        index_columns.sort()
        for i, n in enumerate(index_columns):
            columns.pop(n - i)

    col_dict = {col_name: pd.Series([], dtype=dtype[col_name])
                for col_name in columns}

    return pd.DataFrame(col_dict, columns=columns, index=index)


def check_bins(bins, chromsizes):
    is_cat = pd.api.types.is_categorical(bins['chrom'])
    bins = bins.copy()
    if not is_cat:
        bins['chrom'] = pd.Categorical(
            bins.chrom,
            categories=list(chromsizes.index),
            ordered=True)
    else:
        assert (bins['chrom'].cat.categories == chromsizes.index).all()

    return bins


def balanced_partition(gs, n_chunk_max, file_contigs, loadings=None):
    n_bins = len(gs.bins)
    grouped = gs._bins_grouped

    chrom_nbins = grouped.size()
    if loadings is None:
        loadings = chrom_nbins
    chrmax = loadings.idxmax()
    loadings = loadings / loadings.loc[chrmax]
    const = chrom_nbins.loc[chrmax] / n_chunk_max

    granges = []
    for chrom, group in grouped:
        if chrom not in file_contigs:
            continue
        clen = gs.chromsizes[chrom]
        step = int(np.ceil(const / loadings.loc[chrom]))
        anchors = group.start.values[::step]
        if anchors[-1] != clen:
            anchors = np.r_[anchors, clen]
        granges.extend( (chrom, start, end)
            for start, end in zip(anchors[:-1], anchors[1:]))
    return granges


class GenomeSegmentation(object):
    def __init__(self, chromsizes, bins):
        bins = check_bins(bins, chromsizes)
        self._bins_grouped = bins.groupby('chrom', sort=False)
        nbins_per_chrom = self._bins_grouped.size().values

        self.chromsizes = chromsizes
        self.binsize = get_binsize(bins)
        self.contigs = list(chromsizes.keys())
        self.bins = bins
        self.idmap = pd.Series(
            index=chromsizes.keys(),
            data=range(len(chromsizes)))
        self.chrom_binoffset = np.r_[0, np.cumsum(nbins_per_chrom)]
        self.chrom_abspos = np.r_[0, np.cumsum(chromsizes.values)]
        self.start_abspos = (self.chrom_abspos[bins['chrom'].cat.codes] +
                             bins['start'].values)

    def fetch(self, region):
        chrom, start, end = parse_region(region, self.chromsizes)
        result = self._bins_grouped.get_group(chrom)
        if start > 0 or end < self.chromsizes[chrom]:
            lo = result['end'].values.searchsorted(start, side='right')
            hi = lo + result['start'].values[lo:].searchsorted(end, side='left')
            result = result.iloc[lo:hi]
        return result


def buffered(chunks, size=10000000):
    """
    Take an incoming iterator of small data frame chunks and buffer them into
    an outgoing iterator of larger chunks.

    Parameters
    ----------
    chunks : iterator of :py:class:`pandas.DataFrame`
        Each chunk should have the same column names.
    size : int
        Minimum length of output chunks.

    Yields
    ------
    Larger outgoing :py:class:`pandas.DataFrame` chunks made from concatenating
    the incoming ones.

    """
    buf = []
    n = 0
    for chunk in chunks:
        n += len(chunk)
        buf.append(chunk)
        if n > size:
            yield pd.concat(buf, axis=0)
            buf = []
            n = 0
    if len(buf):
        yield pd.concat(buf, axis=0)
