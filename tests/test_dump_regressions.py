"""
Regression tests for bugs in cooler dump, cooler ls and multi-column queries.
"""

import gzip
import io

import h5py
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

import cooler
from cooler.cli.dump import dump
from cooler.cli.fileops import ls
from cooler.core import CSRReader, DirectRangeQuery2D, FillLowerRangeQuery2D

CHROMSIZES = pd.Series({"chr1": 100, "chr2": 60})


def make_cooler(path):
    rng = np.random.default_rng(0)
    bins = cooler.binnify(CHROMSIZES, 10)
    n = len(bins)
    i, j = rng.integers(0, n, 100), rng.integers(0, n, 100)
    px = pd.DataFrame({"bin1_id": np.minimum(i, j), "bin2_id": np.maximum(i, j)})
    px = px.drop_duplicates().sort_values(["bin1_id", "bin2_id"])
    px["count"] = rng.integers(1, 9, len(px))
    px["foo"] = rng.random(len(px))
    cooler.create_cooler(path, bins, px, columns=["count", "foo"])
    cooler.balance_cooler(
        cooler.Cooler(path), store=True, mad_max=0, min_nnz=0, ignore_diags=0
    )
    return path


def run_dump(args):
    result = CliRunner().invoke(dump, args)
    assert result.exit_code == 0, result.output
    return pd.read_csv(io.StringIO(result.output), sep="\t")


def test_dump_all_value_columns(tmp_path):
    path = make_cooler(str(tmp_path / "c.cool"))
    expected = cooler.Cooler(path).pixels()[:]
    out = run_dump([path, "-H"])
    assert list(out.columns) == ["bin1_id", "bin2_id", "count", "foo"]
    np.testing.assert_allclose(out.to_numpy(), expected.to_numpy(), rtol=1e-5)


def test_dump_columns(tmp_path):
    path = make_cooler(str(tmp_path / "c.cool"))
    expected = cooler.Cooler(path).pixels()[:]

    out = run_dump([path, "-H", "-c", "bin1_id,foo"])
    assert list(out.columns) == ["bin1_id", "foo"]
    np.testing.assert_allclose(out["foo"], expected["foo"], rtol=1e-5)

    # Balancing reads the counts without writing them
    out = run_dump([path, "-H", "-c", "bin1_id,bin2_id", "--balanced"])
    assert list(out.columns) == ["bin1_id", "bin2_id", "balanced"]

    out = run_dump([path, "-H", "-c", "foo", "--join"])
    assert list(out.columns) == [
        "chrom1",
        "start1",
        "end1",
        "chrom2",
        "start2",
        "end2",
        "foo",
    ]

    result = CliRunner().invoke(dump, [path, "-c", "nope"])
    assert result.exit_code == 2
    assert "nope" in result.output


def test_dump_one_based_ids(tmp_path):
    path = make_cooler(str(tmp_path / "c.cool"))
    expected = cooler.Cooler(path).pixels()[:]
    out = run_dump([path, "-H", "--one-based-ids"])
    np.testing.assert_array_equal(out["bin1_id"], expected["bin1_id"] + 1)
    np.testing.assert_array_equal(out["bin2_id"], expected["bin2_id"] + 1)


def test_dump_closes_files(tmp_path):
    path = make_cooler(str(tmp_path / "c.cool"))
    out = str(tmp_path / "out.tsv.gz")
    before = h5py.h5f.get_obj_count(h5py.h5f.OBJ_ALL, h5py.h5f.OBJ_FILE)
    result = CliRunner().invoke(dump, [path, "-o", out])
    assert result.exit_code == 0, result.output
    assert h5py.h5f.get_obj_count(h5py.h5f.OBJ_ALL, h5py.h5f.OBJ_FILE) == before
    with gzip.open(out, "rt") as f:
        assert sum(1 for _ in f) == len(cooler.Cooler(path).pixels())


@pytest.mark.parametrize("engine", [DirectRangeQuery2D, FillLowerRangeQuery2D])
def test_query_several_fields(tmp_path, engine):
    path = make_cooler(str(tmp_path / "c.cool"))
    clr = cooler.Cooler(path)
    with clr.open("r") as h5:
        reader = CSRReader(h5["pixels"], h5["indexes/bin1_offset"][:])
        bbox = (3, 12, 5, 16)
        query = engine(reader, ["foo", "count"], bbox, 7)
        frame = query.to_frame()
        assert list(frame.columns) == ["bin1_id", "bin2_id", "foo", "count"]
        for field in ["foo", "count"]:
            single = engine(reader, field, bbox, 7).to_frame()
            np.testing.assert_array_equal(frame[field], single[field])
            np.testing.assert_array_equal(frame["bin1_id"], single["bin1_id"])
        with pytest.raises(ValueError, match="single value column"):
            query.to_array()
        empty = engine(reader, ["foo", "count"], (3, 3, 5, 16), 7).to_frame()
        assert list(empty.columns) == ["bin1_id", "bin2_id", "foo", "count"]


def test_ls_uri(tmp_path):
    src = make_cooler(str(tmp_path / "c.cool"))
    mcool = str(tmp_path / "z.mcool")
    cooler.zoomify_cooler(src, mcool, [10, 20, 40], chunksize=50)

    def listing(uri):
        result = CliRunner().invoke(ls, [uri])
        assert result.exit_code == 0, result.output
        return [line.split("::")[1] for line in result.output.splitlines()]

    assert listing(mcool) == ["/resolutions/10", "/resolutions/20", "/resolutions/40"]
    assert listing(mcool + "::resolutions") == listing(mcool)
    assert listing(mcool + "::/resolutions/20") == ["/resolutions/20"]
    assert listing(mcool + "::/resolutions/2") == []  # not a prefix match
