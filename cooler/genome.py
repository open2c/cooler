# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals
from collections import namedtuple
from functools import partial
import six
import re

import numpy as np
import pandas

from .exception import RegionParseError
from .util import natsorted, read_tsv, atoi


__all__ = ['read_chromsizes', 'binnify', 'Region']


def read_chromsizes(filepath_or_fp, all_seqs=False,
                    name_patterns=(r'^chr[0-9]+$', r'^chr[XY]$', r'^chrM$')):
    """
    Parse a ``<db>.chrom.sizes`` or ``<db>.chromInfo.txt`` file from the UCSC 
    database, where ``db`` is a genome assembly name.

    Input
    -----
    filepath_or_fp : str or file-like
        ``<db>.chrom.sizes`` text file
    all_seqs : bool, optional
        Whether to return all scaffolds listed in the file. Default is 
        ``False``.
    name_patterns: sequence, optional
        Sequence of regular expressions to capture desired sequence names.
        Each corresponding set of records will be sorted in natural order.

    Returns
    -------
    Data frame indexed by sequence name, with columns 'name' and 'length'.

    """
    chrom_table = read_tsv(filepath_or_fp, 
                  usecols=[0, 1], names=['name', 'length'])
    if all_seqs:
        chrom_table = chrom_table.reindex(index=chrom_table['name'])
    else:
        parts = []
        for pattern in name_patterns:
            part = chrom_table[chrom_table['name'].str.contains(pattern)]
            part.index = part['name']
            part = part.reindex(index=natsorted(part.index))
            parts.append(part)
        chrom_table = pandas.concat(parts, axis=0)
    return chrom_table


def binnify(chrom_table, binsize):
    """
    Divide a genome into uniformly sized bins.

    Input
    -----
    chrom_table : DataFrame
        data frame indexed by chromosome name with chromosome lengths in bp.
    binsize : int
        size of bins in bp

    Returns
    -------
    Data frame with columns: 'chrom', 'start', 'end'.

    """
    def _binnify_each(chrom):
        clen = chrom_table['length'].at[chrom]
        n_bins = int(np.ceil(clen / binsize))
        binedges = np.arange(0, (n_bins+1)) * binsize
        binedges[-1] = clen
        return pandas.DataFrame({
                'chrom': [chrom]*n_bins,
                'start': binedges[:-1],
                'end': binedges[1:],
            }, columns=['chrom', 'start', 'end'])
    return pandas.concat(
        map(_binnify_each, chrom_table.index),
        axis=0, ignore_index=True)


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


def _expect_region(tokens):
    try:
        token = next(tokens)
    except StopIteration:
        raise RegionParseError

    if token[0] != 'CHROM':
        raise RegionParseError
    chrom = token[1]
    try:
        token = next(tokens)
    except StopIteration:
        return (chrom, None, None)

    if token[0] != 'COLON':
        raise RegionParseError
    try:
        token = next(tokens)
        if token[0] != 'COORD':
            raise RegionParseError
        start = atoi(token[1])

        token = next(tokens)
        if token[0] != 'HYPHEN':
            raise RegionParseError

        token = next(tokens)
        if token[0] != 'COORD':
            raise RegionParseError
        end = atoi(token[1])
    except StopIteration:
        raise RegionParseError
    
    if end < start:
        raise RegionParseError
    return chrom, start, end


def parse_region_string(s):
    return _expect_region(_tokenize(s))


class Region(namedtuple('Region', 'chrom start end')):
    """
    Genomic regions are represented as half-open intervals (0-based starts,
    1-based ends) along the length coordinate of an assembled sequence. 

    Parameters
    ----------
    reg : str or tuple
        Genomic region string, or 
        Triple (chrom, start, end), where ``start`` or ``end`` may be ``None``.
    chromsizes : mapping, optional
        Lookup table of scaffold lengths to check against ``chrom`` and the 
        ``end`` coordinate. Required if ``end`` is not supplied.
    
    Returns
    -------
    A well-formed genomic region triple (str, int, int)
    
    """
    def __new__(cls, reg, chromsizes=None):
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
            if clen is None:
                raise ValueError("Cannot determine end coordinate.")
            end = clen
        if end < start:
            raise ValueError("End cannot be less than start")
        if start < 0 or (clen is not None and end > clen):
            raise ValueError(
                "Genomic region out of bounds: [0, {})".format(clen))
        return super(Region, cls).__new__(cls, chrom, start, end)


    ##########################
    ### Binary logical ops ###
    ##########################
    def comes_before(self, other, strict=False):
        if self.chrom != other.chrom or self.start >= other.start: return False
        if self.start < other.start:
            if strict:
                return self.end <= other.start
            else:
                return self.end <= other.end

    def comes_after(self, other, strict=False):
        if self.chrom != other.chrom or self.end <= other.end: return False
        if self.end > other.end:
            if strict:
                return self.start >= other.end
            else:
                return self.start >= other.start

    def contains(self, other, strict=False):
        if (self.chrom != other.chrom 
                or self.start > other.start 
                or self.end < other.end): return False
        if (strict and 
                (self.start == other.start 
                    or self.end == other.end)): return False
        return self.start <= other.start and self.end >= other.end

    def overlaps(self, other):
        if self.chrom != other.chrom: return False
        return (not self.comes_before(other, strict=True) 
                and not self.comes_after(other, strict=True))

    comes_strictly_before = lambda self, other: self.comes_before(other, strict=True)
    comes_strictly_after = lambda self, other: self.comes_after(other, strict=True)
    strictly_contains = lambda self, other: self.contains(other, strict=True)


    ######################
    ### Binary set ops ###
    ######################
    def intersection(self, other):
        if self.chrom != other.chrom: 
            raise ValueError("Regions are on different chromosomes")
        start, end = max(self.start, other.start), min(self.end, other.end)
        if start > end: raise ValueError("Empty intersection.")
        return Region((self.chrom, start, end))

    def combined(self, other):
        if self.chrom != other.chrom: 
            raise ValueError("Regions are on different chromosomes")
        if not self.overlaps(other): raise ValueError("No overlap")
        start, end = min(self.start, other.start), max(self.end, other.end)
        return Region((self.chrom, start, end))

    def hull(self, other):
        if self.chrom != other.chrom: 
            raise ValueError("Regions are on different chromosomes")
        start, end = min(self.start, other.start), max(self.end, other.end)
        return Region((self.chrom, start, end))

    def diffleft(self, other):
        if self.chrom != other.chrom: 
            raise ValueError("Regions are on different chromosomes")
        if not self.overlaps(other) or self.start > other.start:
            raise ValueError("No difference")
        return Region((self.chrom, self.start, other.start))

    def diffright(self, other):
        if self.chrom != other.chrom: 
            raise ValueError("Regions are on different chromosomes")
        if not self.overlaps(other) or self.end < other.end:
            raise ValueError("No difference")
        return Region((self.chrom, other.end, self.end))
