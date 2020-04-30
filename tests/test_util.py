import os.path as op
import numpy as np
import pandas as pd

from cooler import util
import pytest

testdir = op.realpath(op.dirname(__file__))
datadir = op.join(testdir, "data")


def test_partition():
    p = list(util.partition(0, 9, 2))
    assert p == [(0, 2), (2, 4), (4, 6), (6, 8), (8, 9)]


def test_buffered():
    a = pd.DataFrame(np.random.zeros(4, 3), ['a', 'b', 'c'])
    b = pd.DataFrame(np.random.zeros(3, 3), ['a', 'b', 'c'])
    c = pd.DataFrame(np.random.zeros(3, 3), ['a', 'b', 'c'])
    it = util.buffered([a, b, c], size=6)
    assert len(next(it)) == 7
    assert len(next(it)) == 3


def test_rlencode():
    s, l, v = util.rlencode([1, 1, 1, 1, 5, 5, 5, 5, 3, 3, 8, 9, 9])
    assert list(s) == [0, 4, 8, 10, 11]
    assert list(l) == [4, 4, 2, 1, 2]
    assert list(v) == [1, 5, 3, 8, 9]


def test_parse_region_string():
    # UCSC-style names
    assert util.parse_region_string("chr21") == ("chr21", None, None)
    assert util.parse_region_string("chr21:1000-2000") == ("chr21", 1000, 2000)
    assert util.parse_region_string("chr21:1,000-2,000") == ("chr21", 1000, 2000)

    # Ensembl style names
    assert util.parse_region_string("6") == ("6", None, None)
    assert util.parse_region_string("6:1000-2000") == ("6", 1000, 2000)
    assert util.parse_region_string("6:1,000-2,000") == ("6", 1000, 2000)

    # FASTA style names
    assert util.parse_region_string("gb|accession|locus") == (
        "gb|accession|locus",
        None,
        None,
    )
    assert util.parse_region_string("gb|accession|locus:1000-2000") == (
        "gb|accession|locus",
        1000,
        2000,
    )
    assert util.parse_region_string("gb|accession|locus:1,000-2,000") == (
        "gb|accession|locus",
        1000,
        2000,
    )

    # Punctuation in names (aside from :)
    assert util.parse_region_string("name-with-hyphens-") == (
        "name-with-hyphens-",
        None,
        None,
    )
    assert util.parse_region_string("GL000207.1") == ("GL000207.1", None, None)
    assert util.parse_region_string("GL000207.1:1000-2000") == ("GL000207.1", 1000, 2000)

    # Trailing dash
    assert util.parse_region_string("chr21:1000-") == ("chr21", 1000, None)

    # Humanized units
    assert util.parse_region_string("6:1kb-2kb") == ("6", 1000, 2000)
    assert util.parse_region_string("6:1k-2000") == ("6", 1000, 2000)
    assert util.parse_region_string("6:1kb-2M") == ("6", 1000, 2000000)
    assert util.parse_region_string("6:1Gb-") == ("6", 1000000000, None)

    with pytest.raises(ValueError):
        util.parse_region_string("chr1:2,000-1,000")  # reverse selection
        util.parse_region_string("chr1::1000-2000")  # more than one colon


def test_parse_region():
    chromsizes = util.read_chromsizes(op.join(datadir, 'toy.chrom.sizes'))
    assert util.parse_region('chr1:0-10') == ('chr1', 0, 10)
    assert util.parse_region('chr1:0-') == ('chr1', 0, chromsizes['chr1'])


def test_read_chromsizes():
    util.read_chromsizes(op.join(datadir, 'toy.chrom.sizes'))


def test_fetch_chromsizes():
    util.fetch_chromsizes('hg19')


def test_load_fasta():
    fa = util.load_fasta(['chr1', 'chr2'], op.join(datadir, 'toy.fasta'))
    assert len(fa['chr1']) == 32
    assert len(fa['chr2']) == 32


def test_binnify():
    chromsizes = util.read_chromsizes(op.join(datadir, 'toy.chrom.sizes'))
    bins = util.binnify(chromsizes, 10)
    assert len(bins) == 8


def test_digest():
    fa = util.load_fasta(['chr1', 'chr2'], op.join(datadir, 'toy.fasta'))
    bins = util.digest(fa, 'HindIII')
    assert len(bins) == 2


def test_hdf5_contextmanagers():
    pass


def test_hdf5_attrs_to_jsonable_dict():
    pass


def test_dataframe_meta_getter():
    pass


def test_get_binsize():
    pass


def test_get_chromsizes():
    pass


def test_check_bins():
    pass


def test_genome_segmentation():
    pass
