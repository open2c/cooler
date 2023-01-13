import os.path as op
from io import BytesIO

import h5py
import numpy as np
import pandas as pd
import pytest

from cooler import util

testdir = op.realpath(op.dirname(__file__))
datadir = op.join(testdir, "data")


def test_partition():
    p = list(util.partition(0, 9, 2))
    assert p == [(0, 2), (2, 4), (4, 6), (6, 8), (8, 9)]


def test_buffered():
    a = pd.DataFrame(np.random.random((4, 3)), columns=['a', 'b', 'c'])
    b = pd.DataFrame(np.random.random((3, 3)), columns=['a', 'b', 'c'])
    c = pd.DataFrame(np.random.random((3, 3)), columns=['a', 'b', 'c'])
    it = util.buffered([a, b, c], size=6)
    assert len(next(it)) == 7
    assert len(next(it)) == 3


def test_rlencode():
    s, l, v = util.rlencode([1, 1, 1, 1, 5, 5, 5, 5, 3, 3, 8, 9, 9])
    assert list(s) == [0, 4, 8, 10, 11]
    assert list(l) == [4, 4, 2, 1, 2]
    assert list(v) == [1, 5, 3, 8, 9]

    s, l, v = util.rlencode([])
    assert list(s) == []
    assert list(l) == []
    assert list(v) == []


def test_parse_cooler_uri():
    for uri in [
        '/foo/bar/baz.mcool::resolutions/1000',
        '/foo/bar/baz.mcool::/resolutions/1000'
    ]:
        fp, gp = util.parse_cooler_uri(uri)
        assert fp == '/foo/bar/baz.mcool'
        assert gp == '/resolutions/1000'

    for uri in [
        '/foo/bar/baz.cool',
        '/foo/bar/baz.cool::/'
    ]:
        fp, gp = util.parse_cooler_uri(uri)
        assert fp == '/foo/bar/baz.cool'
        assert gp == '/'

    for uri in [
        '/foo/bar/baz.cool::/a/b::c.cool',
    ]:
        with pytest.raises(ValueError):
            util.parse_cooler_uri(uri)


def test_atoi():
    assert util.atoi('1,000') == 1000
    with pytest.raises(ValueError):
        assert util.atoi('1,000.05')  # not an integer


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

    # Bad inputs
    for region in [
        "chr1:2,000-1,000",  # reverse selection
        "chr1::1000-2000",  # more than one colon
        "chr1:1kb-2kDa",  # unknown unit kDa
        "chr1:1000",  # missing end
        "chr1:-2000",  # missing start
        ":1000-2000",  # missing chromosome name
        'chr1:$100-300',  # invalid token
    ]:
        with pytest.raises(ValueError):
            util.parse_region_string(region)


def test_parse_region():
    chromsizes = util.read_chromsizes(op.join(datadir, 'toy.chrom.sizes'))
    assert util.parse_region(('chr1', 0, 10)) == ('chr1', 0, 10)
    assert util.parse_region('chr1:0-10') == ('chr1', 0, 10)
    assert util.parse_region('chr1:0-', chromsizes) == ('chr1', 0, chromsizes['chr1'])

    # Don't accept undefined end unless chromsizes exists
    # NOTE: parse_region_string works here
    with pytest.raises(ValueError):
        util.parse_region('chr1:0-')

    # catch end < start in non-string case
    with pytest.raises(ValueError):
        util.parse_region(('chr1', 10, 0))

    # catch errors when chromsizes is given
    for region in [
        ('chr1', 0, 1000),
        ('chr1', -5, 10),
        ('DoesNotExist', 0, 10),
        'DoesNotExist',
    ]:
        with pytest.raises(ValueError):
            util.parse_region(region, chromsizes)


def test_natsort():
    chroms_alpha = ['chr1', 'chr10', 'chr2', 'chr3']
    chroms_nat = ['chr1', 'chr2', 'chr3', 'chr10']
    assert util.natsorted(chroms_alpha) == chroms_nat
    assert list(util.argnatsort(chroms_alpha)) == [0, 2, 3, 1]


def test_read_chromsizes():
    util.read_chromsizes(op.join(datadir, 'toy.chrom.sizes'))


def test_fetch_chromsizes():
    util.fetch_chromsizes('hg19')


def test_load_fasta():
    fa = util.load_fasta(['chr1', 'chr2'], op.join(datadir, 'toy.fasta'))
    assert len(fa['chr1']) == 32
    assert len(fa['chr2']) == 32

    with pytest.raises(ValueError):
        util.load_fasta(['chr1', 'chr2'])

    # s1 = StringIO(">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    # s2 = StringIO(">chr2\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")
    # fa = util.load_fasta(['chr1', 'chr2'], s1, s2)
    # assert len(fa['chr1']) == 32
    # assert len(fa['chr2']) == 31


def test_binnify():
    chromsizes = util.read_chromsizes(op.join(datadir, 'toy.chrom.sizes'))
    bins = util.binnify(chromsizes, 10)
    assert len(bins) == 8


def test_digest():
    fa = util.load_fasta(['chr1', 'chr2'], op.join(datadir, 'toy.fasta'))
    bins = util.digest(fa, 'HindIII')
    assert len(bins) == 2

    with pytest.raises(ValueError):
        util.digest(fa, 'HindMCMXCIX')


def test_get_binsize():
    chromsizes = util.read_chromsizes(op.join(datadir, 'toy.chrom.sizes'))
    bins = util.binnify(chromsizes, 10)
    assert util.get_binsize(bins) == 10

    # variable-sized bins
    bins = pd.read_csv(
        op.join(datadir, 'toy.bins.var.bed'),
        names=['chrom', 'start', 'end'],
        sep='\t'
    )
    assert util.get_binsize(bins) is None

    # ambiguous case: one bin per chromosome with different lengths
    bins = pd.DataFrame({
        'chrom': ['chr1', 'chr2', 'chr3'],
        'start': [0, 0, 0],
        'end': [100, 200, 300]
    })
    assert util.get_binsize(bins) is None


def test_get_chromsizes():
    chromsizes = util.read_chromsizes(op.join(datadir, 'toy.chrom.sizes'))
    bins = util.binnify(chromsizes, 10)
    assert np.allclose(util.get_chromsizes(bins), chromsizes)


def test_bedslice():
    chromsizes = util.read_chromsizes(op.join(datadir, 'toy.chrom.sizes'))
    bins = util.binnify(chromsizes, 10)
    grouped = bins.groupby('chrom')
    df = util.bedslice(grouped, chromsizes, 'chr1:0-12')
    assert df['chrom'].tolist() == ['chr1', 'chr1']
    assert df['start'].tolist() == [0, 10]


def test_cmd_exists():
    util.cmd_exists('ls')


def test_mad():
    from scipy.stats import median_abs_deviation
    x = np.arange(50)
    assert np.isclose(util.mad(x), median_abs_deviation(x, scale=1))


def test_hdf5_contextmanagers():
    path = op.join(datadir, 'toy.symm.upper.2.cool')

    # file path creates managed handle that gets closed on teardown
    with util.open_hdf5(path) as f:
        pass
    assert not f.id

    # allow appendable open file to pass through with mode='r'
    # might be good to raise a warning
    f = h5py.File(path, 'r+')
    with util.open_hdf5(f, 'r'):
        pass
    assert f.id
    f.close()

    # open file passes through without getting closed on teardown
    f = h5py.File(path, 'r')
    with util.open_hdf5(f):
        pass
    assert f.id

    # can't change mode on open file
    with pytest.raises(ValueError):
        with util.open_hdf5(f, 'r+'):
            pass

    # not allowed on open files
    for mode in ['w', 'w-', 'x']:
        with pytest.raises(ValueError):
            with util.open_hdf5(f, mode):
                pass

    # group's parent file gets closed on teardown
    with util.closing_hdf5(f['chroms']):
        pass
    assert not f.id

    # closing works as a standalone object, not only as a contextmanager
    f = h5py.File(path, 'r')
    grp = util.closing_hdf5(f['chroms'])
    grp.close()
    assert not f.id


def test_hdf5_attrs_to_jsonable_dict():
    b = BytesIO()
    f = h5py.File(b, 'a')
    f.attrs['a'] = np.array([1, 2, 3])
    f.attrs['b'] = 'hello'
    f.attrs['c'] = 3
    dct = util.attrs_to_jsonable(f.attrs)
    assert dct['a'] == [1, 2, 3]
    assert dct['b'] == 'hello'
    assert dct['c'] == 3


def test_check_bins():
    chromsizes = util.read_chromsizes(op.join(datadir, 'toy.chrom.sizes'))
    bins = util.binnify(chromsizes, 10)
    bins['chrom'] = bins['chrom'].astype(str)
    bins = util.check_bins(bins, chromsizes)
    assert pd.api.types.is_categorical_dtype(bins["chrom"])


def test_genome_segmentation():
    chromsizes = util.read_chromsizes(op.join(datadir, 'toy.chrom.sizes'))
    bins = util.binnify(chromsizes, 10)
    gs = util.GenomeSegmentation(chromsizes, bins)
    df = gs.fetch('chr1')
    assert len(df) == 4
    df = gs.fetch('chr1:2-30')
    assert len(df) == 3
    util.balanced_partition(gs, 2, ['chr1'])


def test_dataframe_meta():
    df = pd.DataFrame({
        'a': [1, 2, 3],
        'b': [4., 5., 6.]
    })
    util.infer_meta(df)
    # meta2 = util.get_meta(df.columns, df.dtypes)
    # assert meta1 == meta2
