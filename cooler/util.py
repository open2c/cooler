# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals
import six
import re

import numpy as np
import pandas


def atoi(s):
    return int(s.replace(',', ''))


def parse_region_string(s):

    def _tokenize(s):
        token_spec = [
            ('COORD',   r'[0-9,]+'),
            ('CHROM',   r'\w+'),
            ('COLON',   r':'),
            ('HYPHEN',  r'-'),
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
            raise ValueError
        else:
            if token[0] != expected:
                raise ValueError
        return token[1]

    def _expect(tokens):
        chrom = _check_next(tokens, 'CHROM')

        try:
            token = next(tokens)
        except StopIteration:
            return (chrom, None, None)
        if token[0] != 'COLON':
            raise ValueError

        start = atoi(_check_next(tokens, 'COORD'))
        _check_next(tokens, 'HYPHEN')
        end = atoi(_check_next(tokens, 'COORD'))
        if end < start:
            raise ValueError

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
        start, end = map(int, (start, end))
    try:
        clen = chromsizes[chrom] if chromsizes is not None else None
    except KeyError:
        raise ValueError("Unknown scaffold {}".format(chrom))
    start = 0 if start is None else start
    if end is None:
        if clen is None:  # TODO --- remove?
            raise ValueError("Cannot determine end coordinate.")
        end = clen
    if end < start:
        raise ValueError("End cannot be less than start")
    if start < 0 or (clen is not None and end > clen):
        raise ValueError(
            "Genomic region out of bounds: [0, {})".format(clen))
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


def read_chrominfo(filepath_or,
                   name_patterns=(r'^chr[0-9]+$', r'^chr[XY]$', r'^chrM$'),
                   name_index=True,
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
    name_index : bool, optional
        Index table by chromosome name.
    all_names : bool, optional
        Whether to return all scaffolds listed in the file. Default is
        ``False``.

    Returns
    -------
    Data frame indexed by sequence name, with columns 'name' and 'length'.

    """
    if isinstance(filepath_or, six.string_types) and filepath_or.endswith('.gz'):
        kwargs.setdefault('compression', 'gzip')
    chromtable = pandas.read_csv(
        filepath_or, sep='\t', usecols=[0, 1],
        names=['name', 'length'], **kwargs)
    if not all_names:
        parts = []
        for pattern in name_patterns:
            part = chromtable[chromtable['name'].str.contains(pattern)]
            part = part.iloc[argnatsort(part['name'])]
            parts.append(part)
        chromtable = pandas.concat(parts, axis=0)
    chromtable.insert(0, 'id', np.arange(len(chromtable)))
    if name_index:
        chromtable.index = chromtable['name'].values
    return chromtable


def make_bintable(chromsizes, binsize):
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
    Data frame with columns: 'chrom', 'start', 'end'.

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
    bintable = pandas.concat(map(_each, chromsizes.index),
                             axis=0, ignore_index=True)
    return bintable


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
