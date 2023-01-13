import os.path as op

import h5py
import numpy as np
import pandas
import pandas as pd
import pytest
import scipy.sparse as sps

from cooler import api

testdir = op.realpath(op.dirname(__file__))
datadir = op.join(testdir, "data")


def test_info(mock_cooler):
    info = api.info(mock_cooler)
    assert isinstance(info, dict)


def test_chromtable(mock_cooler):
    table = api.chroms(mock_cooler)
    assert np.all(table["length"] == mock_cooler["chroms"]["length"])


def test_bintable(mock_cooler):
    chromID_lookup = pd.Series({"chr1": 0, "chr2": 1})
    lo, hi = 2, 10
    table = api.bins(mock_cooler, lo, hi)
    assert np.all(chromID_lookup[table["chrom"]] == mock_cooler["bins"]["chrom"][lo:hi])
    assert np.all(table["start"] == mock_cooler["bins"]["start"][lo:hi])
    assert np.all(table["end"] == mock_cooler["bins"]["end"][lo:hi])

    table = api.bins(mock_cooler, lo, hi, fields=["start", "end"])
    assert np.all(table["start"] == mock_cooler["bins"]["start"][lo:hi])
    assert np.all(table["end"] == mock_cooler["bins"]["end"][lo:hi])


def test_bintable_many_contigs():
    # In a file with many contigs, bins/chrom does not have an ENUM header,
    # so chromosome names are taken from the chroms/name
    clr = api.Cooler(op.join(datadir, "manycontigs.1.cool"))
    bins = clr.bins()[:10]
    assert pd.api.types.is_categorical_dtype(bins["chrom"].dtype)

    bins = clr.bins()[['chrom', 'start']][:10]
    assert pd.api.types.is_categorical_dtype(bins["chrom"].dtype)

    chroms = clr.bins()['chrom'][:10]
    clr.bins()['start'][:10]
    assert pd.api.types.is_categorical_dtype(chroms.dtype)


def test_pixeltable(mock_cooler):
    lo, hi = 2, 10
    table = api.pixels(mock_cooler, lo, hi, join=False)
    assert np.all(table["bin1_id"] == mock_cooler["pixels"]["bin1_id"][lo:hi])
    assert np.all(table["bin2_id"] == mock_cooler["pixels"]["bin2_id"][lo:hi])

    table = api.pixels(mock_cooler, lo, hi, join=True)
    assert table.shape == (hi - lo, len(mock_cooler["pixels"]) + 4)


def test_annotate(mock_cooler):
    clr = api.Cooler(mock_cooler)

    # works with full bin table / view or only required bins
    df = clr.matrix(as_pixels=True, balance=False).fetch("chr1")
    df1 = api.annotate(df, clr.bins()[:])
    df2 = api.annotate(df, clr.bins())
    df3 = api.annotate(df, clr.bins().fetch("chr1"))
    assert np.all(df1 == df2)
    assert np.all(df1 == df3)

    # works on empty dataframe
    df4 = api.annotate(df[0:0], clr.bins()[:])
    assert np.all(df4.columns == df3.columns)
    assert len(df4) == 0


def test_annotate_with_partial_bins():
    # Addresses a bug where partial bin-table dataframes were sliced incorrectly,
    # specifically when the pixel dataframe happened to be shorter than it. This
    # led to incorrect NaNs in the join output.
    #
    # This is different from the case where there are pixels from bins that do
    # not appear in the provided partial bin-table dataframe. This will lead to
    # NaNs in the join but the result will be correct because those bins were
    # missing from the input. However, we may want to raise an error in such
    # cases, or disallow partial bin-table inputs entirely.
    clr = api.Cooler(op.join(datadir, "hg19.GM12878-MboI.matrix.2000kb.cool"))
    pix = clr.matrix(as_pixels=True, balance=False).fetch("chr2").iloc[:50]

    bins_chr2 = clr.bins().fetch("chr2")
    assert len(bins_chr2) > len(pix)

    out = api.annotate(pix, bins_chr2)

    for col in ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']:
        assert out[col].notnull().all()


def test_matrix():
    clr = api.Cooler(op.join(datadir, "yeast.10kb.cool"))
    region = ("chrI:100345-220254", "chrII:200789-813183")
    # numpy array
    clr.matrix(balance=False).fetch(*region)
    clr.matrix(balance=True).fetch(*region)
    clr.matrix(balance="weight").fetch(*region)
    clr.matrix(balance="weight", divisive_weights=True).fetch(*region)
    # sparse coo_matrix
    clr.matrix(sparse=True, balance=False).fetch(*region)
    clr.matrix(sparse=True, balance=True).fetch(*region)
    clr.matrix(sparse=True, balance="weight").fetch(*region)
    clr.matrix(
        sparse=True, balance="weight", divisive_weights=True
    ).fetch(*region)
    # dataframe
    clr.matrix(as_pixels=True, join=False, balance=False).fetch(*region)
    clr.matrix(as_pixels=True, join=False, balance=True).fetch(*region)
    clr.matrix(as_pixels=True, join=True, balance=True).fetch(*region)
    clr.matrix(as_pixels=True, join=True, balance=True).fetch(*region)
    clr.matrix(as_pixels=True, join=True, balance="weight").fetch(*region)
    clr.matrix(
        as_pixels=True, join=True, balance="weight", divisive_weights=True
    ).fetch(*region)

    # Unbalanced and asymmetric cooler
    clr = api.Cooler(op.join(datadir, "toy.asymm.2.cool"))
    region = ("chr2", "chr1:2-24")
    # numpy array
    clr.matrix(balance=False).fetch(*region)
    with pytest.raises(ValueError):
        clr.matrix(balance=True).fetch(*region)
    # sparse coo_matrix
    clr.matrix(sparse=True, balance=False).fetch(*region)
    with pytest.raises(ValueError):
        clr.matrix(sparse=True, balance=True).fetch(*region)
    # dataframe
    clr.matrix(as_pixels=True, join=False, balance=False).fetch(*region)
    with pytest.raises(ValueError):
        clr.matrix(
            as_pixels=True, join=True, balance="weight", divisive_weights=True
        ).fetch(*region)


def test_cooler_class(mock_cooler):
    clr = api.Cooler(mock_cooler)
    assert clr.shape == (20, 20)

    # chrom table
    table = clr.chroms()[:]
    assert (
        table["name"].tolist()
        == mock_cooler["chroms"]["name"].astype('U').tolist()
    )
    assert np.all(table["length"] == mock_cooler["chroms"]["length"])

    # bin table
    table = clr.bins().fetch("chr1")
    assert np.all(table["start"] == mock_cooler["bins"]["start"][0:10])
    assert np.all(table["end"] == mock_cooler["bins"]["end"][0:10])

    # pixel table
    table = clr.pixels().fetch("chr1")

    # offsets
    assert clr.offset("chr1") == 0
    assert clr.extent("chr1") == (0, 10)

    # 2D range queries as rectangular or triangular
    A1 = np.triu(clr.matrix(balance=False).fetch("chr2"))
    df = clr.matrix(as_pixels=True, join=False, balance=False).fetch("chr2")
    i0 = clr.offset("chr2")
    i, j, v = df["bin1_id"], df["bin2_id"], df["count"]
    mat = sps.coo_matrix((v, (i - i0, j - i0)), (A1.shape))
    A2 = np.triu(mat.toarray())
    assert np.all(A1 == A2)


def test_cooler_class2():
    path = op.join(datadir, "toy.symm.upper.2.cool")
    with h5py.File(path, 'r') as f:
        clr = api.Cooler(f)
        repr(clr)
        assert clr.root == '/'
        assert clr.filename == path
        assert isinstance(clr.store, h5py.File)
        with clr.open('r') as f:
            pass

    with pytest.raises(KeyError):
        api.Cooler(path + '::/does/not/exist')

    clr = api.Cooler(path)
    clr._load_dset('indexes/chrom_offset')
    clr._load_dset('indexes/bin1_offset')
    clr._load_attrs('bins/chrom')

    with clr.open('r') as f:
        pass

    assert clr.storage_mode == 'symmetric-upper'
    assert clr.binsize == 2
    assert len(clr.chromsizes) == 2
    assert clr.info['nchroms'] == 2
    assert clr.chromnames == ['chr1', 'chr2']
    assert repr(clr) == '<Cooler "{}::{}">'.format('toy.symm.upper.2.cool', '/')
