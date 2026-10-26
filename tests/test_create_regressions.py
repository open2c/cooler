"""
Regression tests for bugs in cooler creation and ingestion.
"""
import os.path as op
import warnings
import zipfile

import h5py
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

import cooler
from cooler.cli.cload import pairs as cload_pairs
from cooler.cli.load import load
from cooler.create import aggregate_records

CHROMSIZES = "chr1\t100\nchr2\t50\n"


@pytest.fixture
def chromsizes_file(tmp_path):
    path = tmp_path / "cs.txt"
    path.write_text(CHROMSIZES)
    return str(path)


def run_cload_pairs(chromsizes_file, pairs_file, out, *args):
    result = CliRunner().invoke(
        cload_pairs,
        [
            "-c1", "2", "-p1", "3", "-c2", "4", "-p2", "5",
            *args,
            chromsizes_file + ":10", pairs_file, out,
        ],
    )
    return result


def dense(clr):
    return clr.matrix(balance=False)[:]


@pytest.mark.parametrize(
    "zero_based, valid_pos, invalid_pos",
    [(False, 100, 101), (True, 99, 100)],
)
def test_cload_pairs_chromosome_end(
    tmp_path, chromsizes_file, zero_based, valid_pos, invalid_pos
):
    args = ["--zero-based"] if zero_based else []

    # The last base of chr1 goes in chr1's last bin, not chr2's first
    pairs = tmp_path / "ok.pairs"
    pairs.write_text(f"r1\tchr1\t{valid_pos}\tchr1\t{valid_pos}\n")
    out = str(tmp_path / "ok.cool")
    result = run_cload_pairs(chromsizes_file, str(pairs), out, *args)
    assert result.exit_code == 0, result.output
    px = cooler.Cooler(out).pixels(join=True)[:]
    assert px[["chrom1", "start1", "chrom2", "start2"]].values.tolist() == [
        ["chr1", 90, "chr1", 90]
    ]

    # One past the end is rejected, on any chromosome
    for chrom in ["chr1", "chr2"]:
        pos = invalid_pos if chrom == "chr1" else invalid_pos - 50
        pairs = tmp_path / "bad.pairs"
        pairs.write_text(f"r1\t{chrom}\t5\t{chrom}\t{pos}\n")
        result = run_cload_pairs(
            chromsizes_file, str(pairs), str(tmp_path / "bad.cool"), *args
        )
        assert result.exit_code != 0
        assert "exceeding chromosome length" in str(result.exception)


def write_value_pairs(path, n=500, seed=1):
    rng = np.random.default_rng(seed)
    p1 = rng.integers(1, 101, n)
    p2 = rng.integers(1, 101, n)
    p1, p2 = np.minimum(p1, p2), np.maximum(p1, p2)
    v = rng.integers(1, 10, n)
    df = pd.DataFrame(
        {"id": [f"r{i}" for i in range(n)], "c1": "chr1", "p1": p1, "c2": "chr1",
         "p2": p2, "s1": "+", "s2": "+", "v": v}
    )
    df.to_csv(path, sep="\t", header=False, index=False)
    return df


@pytest.mark.parametrize("func", ["max", "min", "sum"])
def test_cload_pairs_aggregation_across_chunks(tmp_path, chromsizes_file, func):
    pairs = str(tmp_path / "v.pairs")
    write_value_pairs(pairs)
    results = {}
    for chunksize in [50, 100_000]:
        out = str(tmp_path / f"v{chunksize}.cool")
        result = run_cload_pairs(
            chromsizes_file, pairs, out,
            "--field", f"v=8:agg={func}", "--chunksize", str(chunksize),
        )
        assert result.exit_code == 0, result.output
        results[chunksize] = cooler.Cooler(out).pixels()[:]
    # Chunking doesn't change the result
    pd.testing.assert_frame_equal(results[50], results[100_000])


def test_cload_pairs_mean_across_chunks(tmp_path, chromsizes_file):
    pairs = str(tmp_path / "v.pairs")
    records = write_value_pairs(pairs)
    records.loc[::7, "v"] = np.nan  # missing values are left out of the mean
    records.to_csv(pairs, sep="\t", header=False, index=False)

    results = {}
    for chunksize in [50, 100_000]:
        out = str(tmp_path / f"m{chunksize}.cool")
        result = run_cload_pairs(
            chromsizes_file, pairs, out,
            "--field", "v=8:agg=mean", "--chunksize", str(chunksize),
        )
        assert result.exit_code == 0, result.output
        results[chunksize] = cooler.Cooler(out).pixels()[:]
    pd.testing.assert_frame_equal(results[50], results[100_000])

    # Matches the mean over records of each pixel
    bins1 = (records["p1"] - 1) // 10
    bins2 = (records["p2"] - 1) // 10
    expected = records.groupby([bins1, bins2])["v"].mean().to_numpy()
    np.testing.assert_allclose(results[50]["v"].to_numpy(), expected)


def test_cload_pairs_non_combinable_aggregation(tmp_path, chromsizes_file):
    pairs = str(tmp_path / "v.pairs")
    write_value_pairs(pairs)

    # Fine within a single chunk
    result = run_cload_pairs(
        chromsizes_file, pairs, str(tmp_path / "one.cool"),
        "--field", "v=8:agg=median", "--chunksize", "100000",
    )
    assert result.exit_code == 0, result.output

    # Rejected across chunks
    result = run_cload_pairs(
        chromsizes_file, pairs, str(tmp_path / "many.cool"),
        "--field", "v=8:agg=median", "--chunksize", "50",
    )
    assert result.exit_code != 0
    assert "can't be combined across chunks" in str(result.exception)


@pytest.mark.parametrize("max_merge", [1, 2, 3, 200])
def test_create_unordered_unsorted_chunks(tmp_path, max_merge):
    rng = np.random.default_rng(0)
    bins = cooler.util.binnify(pd.Series({"chr1": 400}), 10)
    n = len(bins)
    mat = np.zeros((n, n), dtype=int)
    chunks = []
    for _ in range(5):
        i = rng.integers(0, n, 40)
        j = rng.integers(0, n, 40)
        df = pd.DataFrame({"bin1_id": np.minimum(i, j), "bin2_id": np.maximum(i, j)})
        df = df.drop_duplicates()
        df["count"] = rng.integers(1, 5, len(df))
        np.add.at(mat, (df["bin1_id"], df["bin2_id"]), df["count"])
        chunks.append(df.sample(frac=1, random_state=1))  # shuffled

    out = str(tmp_path / "u.cool")
    cooler.create_cooler(
        out, bins, iter(chunks), ordered=False, mergebuf=10, max_merge=max_merge
    )
    clr = cooler.Cooler(out)
    px = clr.pixels()[:]
    assert not px.duplicated(["bin1_id", "bin2_id"]).any()
    assert (px.sort_values(["bin1_id", "bin2_id"]).index == px.index).all()
    np.testing.assert_array_equal(np.triu(dense(clr)), mat)


def test_create_unordered_empty_input(tmp_path):
    bins = cooler.util.binnify(pd.Series({"chr1": 100}), 10)
    empty = pd.DataFrame(
        {"bin1_id": np.array([], int), "bin2_id": np.array([], int),
         "count": np.array([], int)}
    )
    for i, chunks in enumerate([[], [empty, empty]]):
        out = str(tmp_path / f"e{i}.cool")
        cooler.create_cooler(out, bins, iter(chunks), ordered=False)
        assert cooler.Cooler(out).info["nnz"] == 0


def test_cload_pairs_no_matching_chromosomes(tmp_path, chromsizes_file):
    pairs = tmp_path / "x.pairs"
    pairs.write_text("r1\tchrX\t5\tchrX\t10\n")
    out = str(tmp_path / "x.cool")
    with warnings.catch_warnings():
        # Unknown chromosomes must not rely on deprecated Categorical behavior
        warnings.filterwarnings("error", message=".*Categorical.*")
        result = run_cload_pairs(chromsizes_file, str(pairs), out)
    assert result.exit_code == 0, result.output
    assert cooler.Cooler(out).info["nnz"] == 0


def test_count_overflow_raises(tmp_path):
    bins = cooler.util.binnify(pd.Series({"chr1": 100}), 10)
    pixels = pd.DataFrame({"bin1_id": [0, 1], "bin2_id": [0, 1], "count": [1, 3e9]})
    pixels["count"] = pixels["count"].astype(np.int64)
    with pytest.raises(ValueError, match="int64"):
        cooler.create_cooler(str(tmp_path / "o.cool"), bins, pixels)

    out = str(tmp_path / "ok.cool")
    cooler.create_cooler(out, bins, pixels, dtypes={"count": "int64"})
    assert cooler.Cooler(out).pixels()["count"][:].tolist() == [1, 3_000_000_000]


def test_cload_pairs_comment_char(tmp_path, chromsizes_file):
    pairs = tmp_path / "c.pairs"
    pairs.write_text("@HD header\n@PG other\nr1\tchr1\t10\tchr1\t50\n")
    out = str(tmp_path / "c.cool")
    result = run_cload_pairs(
        chromsizes_file, str(pairs), out, "--comment-char", "@"
    )
    assert result.exit_code == 0, result.output
    assert cooler.Cooler(out).info["sum"] == 1


def test_cload_pairs_zip_input(tmp_path, chromsizes_file):
    pairs = tmp_path / "z.pairs"
    pairs.write_text("#comment\nr1\tchr1\t10\tchr1\t50\nr2\tchr2\t5\tchr2\t7\n")
    zipped = tmp_path / "z.pairs.zip"
    with zipfile.ZipFile(zipped, "w") as zf:
        zf.write(pairs, "z.pairs")
    out = str(tmp_path / "z.cool")
    result = run_cload_pairs(chromsizes_file, str(zipped), out)
    assert result.exit_code == 0, result.output
    assert cooler.Cooler(out).info["sum"] == 2


def test_load_extra_field_keeps_count(tmp_path, chromsizes_file):
    coo = tmp_path / "f.coo"
    coo.write_text("0\t1\t5\t0.5\n2\t3\t7\t1.5\n")
    out = str(tmp_path / "f.cool")
    result = CliRunner().invoke(
        load,
        ["-f", "coo", "--field", "score=4", chromsizes_file + ":10", str(coo), out],
    )
    assert result.exit_code == 0, result.output
    px = cooler.Cooler(out).pixels()[:]
    assert px["count"].tolist() == [5, 7]
    assert px["score"].tolist() == [0.5, 1.5]


def test_aggregate_records_does_not_modify_agg():
    agg = {"v": "max"}
    aggregate_records(agg=agg, count=True)
    assert agg == {"v": "max"}


def test_create_scool_sorts_and_names(tmp_path):
    bins = cooler.util.binnify(pd.Series({"chr1": 100}), 10)
    rng = np.random.default_rng(0)
    df = pd.DataFrame({"bin1_id": [0, 0, 1, 2, 4], "bin2_id": [0, 3, 1, 5, 9]})
    df["count"] = rng.integers(1, 5, len(df))
    shuffled = df.sample(frac=1, random_state=3)

    # DataFrame pixels are sorted
    path = str(tmp_path / "s.scool")
    cooler.create_scool(path, bins, {"cellA": shuffled})
    clr = cooler.Cooler(path + "::/cells/cellA")
    pd.testing.assert_frame_equal(clr.pixels()[:], df.reset_index(drop=True),
                                  check_dtype=False)

    # Names stored under the same basename are rejected before writing
    path = str(tmp_path / "d.scool")
    with pytest.raises(ValueError, match="same name"):
        cooler.create_scool(path, bins, {"x/a": df, "y/a": df})
    assert not op.exists(path)

    # A non-root scool group works and links its own tables
    path = str(tmp_path / "g.h5")
    cooler.create_scool(path + "::/grp", bins, {"cellA": df})
    clr = cooler.Cooler(path + "::/grp/cells/cellA")
    assert clr.pixels()[:]["count"].sum() == df["count"].sum()
    with h5py.File(path, "r") as f:
        root_addr = h5py.h5o.get_info(f["grp/bins/start"].id).addr
        cell_addr = h5py.h5o.get_info(f["grp/cells/cellA/bins/start"].id).addr
    assert root_addr == cell_addr


def test_field_order_does_not_swap_columns(tmp_path, chromsizes_file):
    # Fields given out of column order must keep their own values
    pairs = tmp_path / "ab.pairs"
    pairs.write_text("r1\tchr1\t10\tchr1\t50\t+\t+\t1\t100\n")
    out = str(tmp_path / "ab.cool")
    result = run_cload_pairs(
        chromsizes_file, str(pairs), out, "--field", "b=9", "--field", "a=8"
    )
    assert result.exit_code == 0, result.output
    px = cooler.Cooler(out).pixels()[:]
    assert px["a"].tolist() == [1.0]
    assert px["b"].tolist() == [100.0]

    coo = tmp_path / "f.coo"
    coo.write_text("0\t1\t5\t0.5\t9\n")
    out = str(tmp_path / "f.cool")
    result = CliRunner().invoke(
        load,
        ["-f", "coo", "--field", "y=5", "--field", "x=4",
         chromsizes_file + ":10", str(coo), out],
    )
    assert result.exit_code == 0, result.output
    px = cooler.Cooler(out).pixels()[:]
    assert (px["count"].tolist(), px["x"].tolist(), px["y"].tolist()) == (
        [5], [0.5], [9.0]
    )


def test_parse_bins_reports_parse_errors(tmp_path):
    from cooler.cli._util import parse_bins

    bed = tmp_path / "bad.bed"
    bed.write_text('chr1\t0\t10\n"chr1\t10\t20\n')
    with pytest.raises(ValueError, match="Failed to parse bins file"):
        parse_bins(str(bed))


def test_cload_tabix_rejects_out_of_range_positions(tmp_path, chromsizes_file):
    pysam = pytest.importorskip("pysam")
    from cooler.cli.cload import tabix as cload_tabix

    def run(lines):
        pairs = tmp_path / "t.pairs"
        pairs.write_text("".join(lines))
        gz = str(pairs) + ".gz"
        pysam.tabix_compress(str(pairs), gz, force=True)
        pysam.tabix_index(gz, seq_col=1, start_col=2, end_col=2, force=True)
        return CliRunner().invoke(
            cload_tabix,
            ["-c2", "4", "-p2", "5", chromsizes_file + ":10", gz,
             str(tmp_path / "t.cool")],
        )

    result = run(["r1\tchr1\t10\tchr1\t100\n"])  # last base of chr1
    assert result.exit_code == 0, result.output
    result = run(["r1\tchr1\t10\tchr1\t101\n"])  # one past the end
    assert result.exit_code != 0
    assert "outside of chr1" in str(result.exception)


def test_create_unordered_mean_of_pixel_occurrences(tmp_path):
    # Generic pixel chunks: each occurrence of a pixel counts as one value
    bins = cooler.util.binnify(pd.Series({"chr1": 100}), 10)
    chunks = [
        pd.DataFrame({"bin1_id": [0, 1], "bin2_id": [1, 2], "count": [1, 1],
                      "v": [2.0, 5.0]}),
        pd.DataFrame({"bin1_id": [0], "bin2_id": [1], "count": [1], "v": [4.0]}),
        pd.DataFrame({"bin1_id": [0], "bin2_id": [1], "count": [1], "v": [np.nan]}),
    ]
    out = str(tmp_path / "m.cool")
    cooler.create_cooler(
        out, bins, iter(chunks), columns=["count", "v"], ordered=False,
        agg={"v": "mean"},
    )
    px = cooler.Cooler(out).pixels()[:]
    assert px["v"].tolist() == [3.0, 5.0]
    assert px["count"].tolist() == [3, 1]


def test_storage_options(tmp_path):
    bins = cooler.util.binnify(pd.Series({"chr1": 100}), 10)
    px = pd.DataFrame({"bin1_id": [0, 1], "bin2_id": [3, 4], "count": [1, 2]})
    for opts in [{"fillvalue": 0}, {"scaleoffset": 0}]:
        path = str(tmp_path / "o.cool")
        cooler.create_cooler(path, bins, px, h5opts=opts)
        assert cooler.Cooler(path).pixels()[:]["count"].tolist() == [1, 2]
    with pytest.raises(ValueError, match="maxshape"):
        cooler.create_cooler(str(tmp_path / "m.cool"), bins, px, h5opts={"maxshape": 9})


@pytest.mark.parametrize("chunksize", ["1", "10"])
def test_load_rejects_duplicates_across_chunks(tmp_path, chromsizes_file, chunksize):
    coo = tmp_path / "d.coo"
    coo.write_text("0\t1\t5\n0\t1\t7\n2\t3\t1\n")
    out = tmp_path / "d.cool"
    result = CliRunner().invoke(
        load,
        ["-f", "coo", "--chunksize", chunksize, chromsizes_file + ":10", str(coo),
         str(out)],
    )
    assert result.exit_code != 0
    assert "duplicate pixel (0, 1)" in str(result.exception)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["cs.txt", "d.coo"]


def test_unordered_duplicates_option(tmp_path):
    bins = cooler.util.binnify(pd.Series({"chr1": 100}), 10)
    a = pd.DataFrame({"bin1_id": [0, 0], "bin2_id": [1, 1], "count": [2, 3]})
    b = pd.DataFrame({"bin1_id": [0], "bin2_id": [1], "count": [4]})

    # Combined by default, within and across chunks
    out = str(tmp_path / "c.cool")
    cooler.create_cooler(out, bins, iter([a, b]), ordered=False)
    assert cooler.Cooler(out).pixels()[:]["count"].tolist() == [9]

    for chunks in [[a], [a.iloc[:1], b]]:
        with pytest.raises(cooler.create.BadInputError, match="duplicate"):
            cooler.create_cooler(
                str(tmp_path / "e.cool"), bins, iter(chunks), ordered=False,
                duplicates="error",
            )


def test_failed_create_leaves_existing_file_intact(tmp_path):
    bins = cooler.util.binnify(pd.Series({"chr1": 100}), 10)
    path = str(tmp_path / "x.cool")
    good = pd.DataFrame({"bin1_id": [0], "bin2_id": [1], "count": [5]})
    cooler.create_cooler(path, bins, good)

    bad = pd.DataFrame({"bin1_id": [0], "bin2_id": [99], "count": [1]})
    with pytest.raises(cooler.create.BadInputError):
        cooler.create_cooler(path, bins, bad)
    assert cooler.Cooler(path).pixels()[:]["count"].tolist() == [5]
    assert sorted(p.name for p in tmp_path.iterdir()) == ["x.cool"]


def test_square_storage_does_not_warn(tmp_path):
    bins = cooler.util.binnify(pd.Series({"chr1": 100}), 10)
    px = pd.DataFrame({"bin1_id": [3], "bin2_id": [1], "count": [1]})
    with warnings.catch_warnings():
        warnings.filterwarnings("error", message=".*triucheck.*")
        cooler.create_cooler(
            str(tmp_path / "s.cool"), bins, px, symmetric_upper=False
        )


def test_skipped_chromosomes_are_reported_once(tmp_path, chromsizes_file):
    import logging

    records = []
    handler = logging.Handler()
    handler.emit = records.append
    logger = logging.getLogger("cooler.create")
    logger.addHandler(handler)
    try:
        pairs = tmp_path / "s.pairs"
        pairs.write_text(
            "r1\tchr1\t5\tchrZ\t10\nr2\tchrZ\t5\tchr1\t9\nr3\tchr1\t5\tchr1\t9\n"
        )
        result = run_cload_pairs(
            chromsizes_file, str(pairs), str(tmp_path / "s.cool"), "--chunksize", "1"
        )
        assert result.exit_code == 0, result.output
    finally:
        logger.removeHandler(handler)
    skipped = [r.getMessage() for r in records if "Skipping" in r.getMessage()]
    assert len(skipped) == 1 and "chrZ" in skipped[0]


def test_validate_pixels_duplicates():
    from cooler.create import BadInputError, validate_pixels

    check = validate_pixels(10, True, True, True, False)
    check({"bin1_id": np.array([0, 1]), "bin2_id": np.array([1, 1])})
    with pytest.raises(BadInputError, match="duplicate"):
        check(pd.DataFrame({"bin1_id": [1, 0, 1], "bin2_id": [1, 1, 1]}))
