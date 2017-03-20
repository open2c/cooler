# -*- coding: utf-8 -*-
from __future__ import division, print_function
from contextlib import contextmanager
from collections import OrderedDict
import six
import re
import os

import numpy as np
import pandas
import h5py


def atoi(s):
    return int(s.replace(',', ''))


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
            ('INT',    r'[0-9,]+'),
            ('ALNUM',  r'[a-zA-z0-9_|]+'),
            ('COLON',  r':'),
            ('HYPHEN', r'-'),
        ]
        tok_regex = r'\s*' + r'|\s*'.join(
            r'(?P<%s>%s)' % pair for pair in token_spec)
        tok_regex = re.compile(tok_regex)
        for match in tok_regex.finditer(s):
            typ = match.lastgroup
            yield typ, match.group(typ)

    def _check_next(tokens, expected):
        try:
            token = next(tokens)
        except StopIteration:
            raise ValueError('Expected {} token missing'.format(expected))
        else:
            if token[0] not in expected:
                raise ValueError('Unexpected token "{}"'.format(token[1]))
        return token[1]

    def _expect(tokens):
        chrom = _check_next(tokens, ['ALNUM', 'INT'])
        try:
            token = next(tokens)
        except StopIteration:
            return (chrom, None, None)
        if token[0] != 'COLON':
            raise ValueError('Got "{}" after chromosome label'.format(token[1]))

        start = atoi(_check_next(tokens, ['INT']))
        _check_next(tokens, ['HYPHEN'])
        end = atoi(_check_next(tokens, ['INT']))
        if end < start:
            raise ValueError('End coordinate less than start')

        return chrom, start, end

    return _expect(_tokenize(s))


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
    Series of integer bp lengths indexed by sequence name.

    """
    if isinstance(filepath_or, six.string_types) and filepath_or.endswith('.gz'):
        kwargs.setdefault('compression', 'gzip')
    chromtable = pandas.read_csv(
        filepath_or, sep='\t', usecols=[0, 1],
        names=['name', 'length'], dtype={'name':str}, **kwargs)
    if not all_names:
        parts = []
        for pattern in name_patterns:
            part = chromtable[chromtable['name'].str.contains(pattern)]
            part = part.iloc[argnatsort(part['name'])]
            parts.append(part)
        chromtable = pandas.concat(parts, axis=0)
    chromtable.index = chromtable['name'].values
    return chromtable['length']


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
    Dataframe with columns: 'chrom', 'start', 'end'.

    """
    def _each(chrom):
        clen = chromsizes[chrom]
        n_bins = int(np.ceil(clen / binsize))
        binedges = np.arange(0, (n_bins+1)) * binsize
        binedges[-1] = clen
        return pandas.DataFrame({
                'chrom': [chrom]*n_bins,
                'start': binedges[:-1],
                'end': binedges[1:],
            }, columns=['chrom', 'start', 'end'])
    
    bintable = pandas.concat(
        map(_each, chromsizes.keys()),
        axis=0, 
        ignore_index=True)
    
    bintable['chrom'] = pandas.Categorical(
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
        Name of restriction enzyme.

    Returns
    -------
    Dataframe with columns: 'chrom', 'start', 'end'.

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

        frags = pandas.DataFrame({
            'chrom': [chrom] * n_frags,
            'start': cuts[:-1],
            'end': cuts[1:]},
            columns=['chrom', 'start', 'end'])
        return frags

    return pandas.concat(map(_each, chroms), axis=0, ignore_index=True)


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
    return pandas.Series(index=chroms, data=lengths)


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
