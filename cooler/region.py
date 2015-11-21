from __future__ import division, print_function, unicode_literals
from collections import namedtuple
import math
import six
import re

import numpy as np
import pandas


Region = namedtuple('Region', 'chrom start end')


class ParseError(Exception):
    pass


def atoi(s):
    return int(s.replace(',', ''))


def tokenize(s):
    token_spec = [
        ('COORD',   r'[0-9,]+'),
        ('CHROM',   r'\w+'), 
        ('COLON',   r':'),   
        ('HYPHEN',  r'-'),
    ]
    tok_regex = r'\s*' + r'|\s*'.join( r'(?P<%s>%s)' % pair for pair in token_spec )
    tok_regex = re.compile(tok_regex)
    for match in tok_regex.finditer(s):
        typ = match.lastgroup
        yield typ, match.group(typ)


def expect_region(tokens):
    try:
        token = next(tokens)
    except StopIteration:
        raise ParseError
    if token[0] != 'CHROM':
        raise ParseError

    chrom = token[1]
    try:
        token = next(tokens)
    except StopIteration:
        return (chrom, None, None)
    if token[0] != 'COLON':
        raise ParseError

    try:
        token = next(tokens)
        if token[0] != 'COORD':
            raise ParseError
        start = atoi(token[1])

        token = next(tokens)
        if token[0] != 'HYPHEN':
            raise ParseError

        token = next(tokens)
        if token[0] != 'COORD':
            raise ParseError
        end = atoi(token[1])
    except StopIteration:
        raise ParseError
    
    if end < start:
        raise ParseError
    return chrom, start, end


def parse_region_string(s):
    return expect_region(tokenize(s))
    

def parse_region(reg, chromsizes=None):
    """
    Input
    -----
    reg : str or tuple
        Genomic region string or (chrom, start, end) triple.
        ``start`` or ``end`` may be ```None``.
    chromsizes : mapping, optional
        Lookup table for contig lengths to check against the
        ``end`` coordinate. Required if ``end`` is not supplied.
    
    Returns
    -------
    Region
        Well-formed genomic region triple (str, int, int)
    
    """
    if isinstance(reg, six.string_types):
        chrom, start, end = parse_region_string(reg)
    else:
        chrom, start, end = reg
        start, end = map(int, (start, end))
    clen = chromsizes[chrom] if chromsizes is not None else None
    start = 0 if start is None else start
    if end is None:
        if clen is None:
            raise ValueError("Cannot determine end coordinate.")
        end = clen
    if start < 0 or (clen is not None and end > clen):
        raise ValueError("Genomic region out of bounds")
    return Region(chrom, start, end)


def comes_before(a, b, strict=False):
    if a.chrom != b.chrom or a.start >= b.start: return False
    if a.start < b.start:
        if strict:
            return a.end <= b.start
        else:
            return a.end <= b.end


def comes_after(a, b, strict=False):
    if a.chrom != b.chrom or a.end <= b.end: return False
    if a.end > b.end:
        if strict:
            return a.start >= b.end
        else:
            return a.start >= b.start


def contains(a, b, strict=False):
    if a.chrom != b.chrom or a.start > b.start or a.end < b.end: return False
    if strict and a.start == b.start or a.end == b.end: return False
    return a.start <= b.start and a.end >= b.end


def overlaps(a, b):
    if a.chrom != b.chrom: return False
    return not comes_before(a, b, strict=True) and not comes_after(a, b, strict=True)


def intersection(a, b):
    if a.chrom != b.chrom: raise ValueError
    start, end = max(a.start, b.start), min(a.end, b.end)
    if start > end: raise ValueError("Empty intersection.")
    return Region(a.chrom, start, end)


def union(a, b):
    if a.chrom != b.chrom: raise ValueError
    if not overlaps(a, b): raise ValueError("Intervals not overlapping")
    start, end = min(a.start, b.start), max(a.end, b.end)
    return Region(a.chrom, start, end)


def hull(a, b):
    if a.chrom != b.chrom: raise ValueError
    start, end = min(a.start, b.start), max(a.end, b.end)
    return Region(a.chrom, min(a.start, b.start), max(a.end, b.end))


def diff(a, b):
    if a.chrom != b.chrom: raise ValueError
    if not overlaps(a, b) or contains(b, a): raise ValueError("Difference is empty")
    return Region(a.chrom, a.start, b.start) if comes_before(a, b) else Region(a.chrom, b.end, a.end)


def partition(a, b):
    if a.chrom != b.chrom: raise ValueError
    return diff(a, b), intersection(a, b), diff(b, a)

