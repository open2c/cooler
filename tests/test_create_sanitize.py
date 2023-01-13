import os.path as op
from io import StringIO

import numpy as np
import pandas as pd
import pytest

import cooler
from cooler.create import (
    BadInputError,
    aggregate_records,
    sanitize_pixels,
    sanitize_records,
    validate_pixels,
)

testdir = op.dirname(op.realpath(__file__))
datadir = op.join(testdir, "data")


columns = [
    "chrom1",
    "pos1",
    "strand1",
    "chrom2",
    "pos2",
    "strand2",
    "name",
    "pair_type",
    "triu",
]

valid_data = """chr1\t1\t+\tchr2\t100\t-\t.\tLL\t1
chr2\t99\t+\tchr1\t13\t-\t.\tLL\t0
chr2\t13\t+\tchr2\t60\t-\t.\tLL\t1
chr1\t200\t+\tchr2\t50\t-\t.\tLL\t1
chr3\t11\t+\tchr3\t40\t-\t.\tLL\t1
chr1\t234\t+\tchr3\t30\t-\t.\tLL\t1
chr3\t3\t+\tchr2\t20\t-\t.\tLL\t0
chr2\t23\t+\tchr3\t11\t-\t.\tLL\t1
chr1\t123\t+\tchr1\t200\t-\t.\tLL\t1
"""

nuisance_chroms = """chr1\t222\t+\tchr9\t200\t-\t.\tLL\t1
chr9\t222\t+\tchr9\t200\t-\t.\tLL\t1"""
oob_lower = """chr1\t-1\t+\tchr1\t10\t+\t.\tLL\t1"""
oob_upper = """chr1\t123\t+\tchr1\t301\t+\t.\tLL\t1"""

binsize = 10
chromsizes = pd.Series(index=["chr1", "chr2", "chr3"], data=[300, 300, 300])
bins = cooler.util.binnify(chromsizes, binsize)


def _insert_lines(d, new):
    lines = d.split("\n")
    for line in new.split("\n"):
        lines.insert(np.random.randint(len(lines)), line)
    return "\n".join(lines)


def test_sanitize_records():

    chunk = pd.read_csv(StringIO(valid_data), sep="\t", names=columns)
    with pytest.raises(ValueError):
        sanitize_records(
            bins,
            schema="doesnotexist",
            validate=True,
            tril_action="reflect",
            is_one_based=True,
            sort=True,
        )(chunk.copy())

    chunk = pd.read_csv(StringIO(valid_data), sep="\t", names=columns)
    sanitize_records(
        bins,
        schema="pairs",
        validate=True,
        tril_action="reflect",
        is_one_based=True,
        sort=True,
    )(chunk.copy())

    # variable-length bins
    chunk = pd.read_csv(StringIO(valid_data), sep="\t", names=columns)
    sanitize_records(
        pd.DataFrame({
            'chrom': ['chr1', 'chr1', 'chr2', 'chr2', 'chr3'],
            'start': [0, 150, 0, 100, 0],
            'end': [150, 300, 100, 300, 300],
        }),
        schema="pairs",
        validate=True,
        tril_action="reflect",
        is_one_based=True,
        sort=True,
    )(chunk.copy())

    # input with already enum-encoded chromosomes (decode_chroms=False)
    text = """0\t1\t+\t1\t100\t-\t.\tLL\t1
1\t99\t+\t0\t13\t-\t.\tLL\t0
1\t13\t+\t1\t60\t-\t.\tLL\t1
0\t200\t+\t1\t50\t-\t.\tLL\t1
2\t11\t+\t2\t40\t-\t.\tLL\t1
0\t234\t+\t2\t30\t-\t.\tLL\t1
2\t3\t+\t1\t20\t-\t.\tLL\t0
1\t23\t+\t2\t11\t-\t.\tLL\t1
0\t123\t+\t-1\t200\t-\t.\tLL\t1
"""
    chunk = pd.read_csv(StringIO(text), sep="\t", names=columns)
    sanitize_records(
        bins,
        schema="pairs",
        decode_chroms=False,
        validate=True,
        tril_action="reflect"
    )(chunk.copy())
    # fails on string chromosomes
    chunk = pd.read_csv(StringIO(valid_data), sep="\t", names=columns)
    with pytest.raises(BadInputError):
        sanitize_records(
            bins,
            schema="pairs",
            decode_chroms=False,
            validate=True,
            tril_action="reflect"
        )(chunk.copy())

    # empty chunk
    out = sanitize_records(
        bins,
        schema="pairs",
        validate=True,
        tril_action="reflect"
    )(chunk.iloc[0:0])
    assert len(out) == 0


def test_sanitize_records_triu_action():
    text = valid_data
    chunk = pd.read_csv(StringIO(text), sep="\t", names=columns)
    out = sanitize_records(bins, schema="pairs", validate=True, tril_action="reflect")(
        chunk.copy()
    )
    is_tril = ~np.array(out["triu"], dtype=bool)
    is_tril_ix = out.index[is_tril]
    assert np.all(out.loc[is_tril_ix, "chrom1"] == chunk.loc[is_tril_ix, "chrom2"])
    assert np.all(out.loc[is_tril_ix, "chrom2"] == chunk.loc[is_tril_ix, "chrom1"])
    assert np.all(out.loc[is_tril_ix, "strand1"] == "+")

    text = valid_data
    chunk = pd.read_csv(StringIO(text), sep="\t", names=columns)
    out = sanitize_records(bins, schema="pairs", validate=True, tril_action="drop")(
        chunk.copy()
    )
    is_tril = ~np.array(out["triu"], dtype=bool)
    is_tril_ix = out.index[is_tril]
    assert np.all(out.loc[is_tril_ix, "chrom1"] == chunk.loc[is_tril_ix, "chrom2"])
    assert np.all(out.loc[is_tril_ix, "chrom2"] == chunk.loc[is_tril_ix, "chrom1"])
    assert np.all(out.loc[is_tril_ix, "strand1"] == "+")
    assert len(out) == chunk["triu"].sum()

    func = sanitize_records(bins, schema="pairs", validate=True, tril_action="raise")
    text = valid_data
    chunk = pd.read_csv(StringIO(text), sep="\t", names=columns)
    with pytest.raises(BadInputError):
        func(chunk)


def test_sanitize_records_with_strand_column():
    text = valid_data
    chunk = pd.read_csv(StringIO(text), sep="\t", names=columns)
    out = sanitize_records(
        bins,
        schema="pairs",
        validate=True,
        tril_action="reflect",
        sided_fields=("chrom", "pos", "strand"),
    )(chunk.copy())
    is_tril = ~np.array(out["triu"], dtype=bool)
    assert np.all(out.loc[is_tril, "chrom1"] == chunk.loc[is_tril, "chrom2"])
    assert np.all(out.loc[is_tril, "chrom2"] == chunk.loc[is_tril, "chrom1"])
    assert np.all(out.loc[is_tril, "strand1"] == "-")


def test_sanitize_records_with_nuisance_records():
    text = _insert_lines(valid_data, nuisance_chroms)
    chunk = pd.read_csv(StringIO(text), sep="\t", names=columns)
    out = sanitize_records(bins, schema="pairs", validate=True, tril_action="reflect")(
        chunk.copy()
    )
    assert ("chr9" not in out["chrom1"]) and ("chr9" not in out["chrom2"])


def test_sanitize_records_with_bad_records():
    func = sanitize_records(bins, schema="pairs", validate=True, tril_action="reflect")

    text = _insert_lines(valid_data, oob_lower)
    chunk = pd.read_csv(StringIO(text), sep="\t", names=columns)
    with pytest.raises(BadInputError):
        func(chunk)

    text = _insert_lines(valid_data, oob_upper)
    chunk = pd.read_csv(StringIO(text), sep="\t", names=columns)
    with pytest.raises(BadInputError):
        func(chunk)


def test_sanitize_pixels():
    bins = cooler.binnify(
        cooler.util.read_chromsizes(op.join(datadir, "toy.chrom.sizes")), 1
    )
    chunk = pd.read_csv(
        op.join(datadir, "toy.symm.upper.1.zb.coo"),
        sep='\t',
        names=['bin1_id', 'bin2_id', 'count']
    )
    chunk['foo1'] = 4
    chunk['foo2'] = 2
    sanitize_pixels(
        bins,
    )(chunk.copy())

    # one-based bin IDs
    out = sanitize_pixels(
        bins,
        is_one_based=True,
    )(chunk.copy())
    assert (out['bin1_id'] == chunk['bin1_id'] - 1).all()

    # tril action: reflect (after swapping bin1, bin2)
    tril_chunk = chunk.copy()
    tril_chunk['bin2_id'] = chunk['bin1_id']
    tril_chunk['bin1_id'] = chunk['bin2_id']
    out = sanitize_pixels(
        bins,
        tril_action="reflect",
        sided_fields=['foo'],
    )(tril_chunk.copy())
    assert len(out) == len(chunk)
    assert (out['foo2'] == chunk['foo1']).all()
    assert (out['foo1'] == chunk['foo2']).all()
    assert (out['bin1_id'] == chunk['bin1_id']).all()
    assert (out['bin2_id'] == chunk['bin2_id']).all()

    # tril action: drop
    out = sanitize_pixels(
        bins,
        tril_action="drop",
    )(tril_chunk.copy())
    assert len(out) == 0

    # tril action: raise
    with pytest.raises(BadInputError):
        sanitize_pixels(
            bins,
            tril_action="raise",
        )(tril_chunk.copy())


def test_validate_pixels():
    bins = cooler.binnify(
        cooler.util.read_chromsizes(op.join(datadir, "toy.chrom.sizes")), 1
    )
    chunk = pd.read_csv(
        op.join(datadir, "toy.symm.upper.1.zb.coo"),
        sep='\t',
        names=['bin1_id', 'bin2_id', 'count']
    )
    validator = validate_pixels(
        len(bins),
        boundscheck=True,
        triucheck=True,
        dupcheck=True,
        ensure_sorted=True
    )
    validator(chunk.copy())
    validator(chunk.to_dict(orient='series'))

    # wrongly assume zero-based, producing -1 bins IDs
    chunk_ = sanitize_pixels(
        bins,
        is_one_based=True,
    )(chunk.copy())
    with pytest.raises(BadInputError):
        validator(chunk_)

    # out-of-bounds bin ID
    chunk_ = chunk.copy()
    chunk_.at[-1, 'bin1_id'] = len(bins) + 1
    with pytest.raises(BadInputError):
        validator(chunk_)

    # pass in non-triu data
    tril_chunk = chunk.copy()
    tril_chunk['bin2_id'] = chunk['bin1_id']
    tril_chunk['bin1_id'] = chunk['bin2_id']
    with pytest.raises(BadInputError):
        validator(tril_chunk)

    # pass in duplicates
    with pytest.raises(BadInputError):
        validator(pd.concat([chunk, chunk], ignore_index=True))


def test_aggregate_records():
    bins = cooler.binnify(
        cooler.util.read_chromsizes(op.join(datadir, "toy.chrom.sizes")), 1
    )
    records = pd.read_csv(
        op.join(datadir, "toy.pairs"),
        sep='\t',
        names=[
            "read_id",
            "chrom1", "pos1",
            "chrom2", "pos2",
            "strand1", "strand2",
            "value"
        ]
    )
    sanitizer = sanitize_records(
        bins,
        schema="pairs",
        validate=False,
        tril_action="reflect",
        is_one_based=False,
        sort=False,
    )
    chunk = sanitizer(records)

    aggregator = aggregate_records()
    aggregator(chunk)
