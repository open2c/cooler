"""
Regression tests for bugs in coarsening, zoomifying and merging.
"""

import h5py
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

import cooler
from cooler._reduce import (
    CoolerMerger,
    coarsen_cooler,
    merge_breakpoints,
    zoomify_cooler,
)
from cooler.cli.coarsen import coarsen as coarsen_cli
from cooler.cli.zoomify import zoomify as zoomify_cli

CHROMSIZES = pd.Series({"chr1": 230, "chr2": 170})


def make_cooler(path, binsize, seed=0, **kwargs):
    """A random symmetric-upper cooler with a float value column `foo`."""
    rng = np.random.default_rng(seed)
    bins = cooler.util.binnify(CHROMSIZES, binsize)
    n = len(bins)
    i = rng.integers(0, n, 300)
    j = rng.integers(0, n, 300)
    px = pd.DataFrame({"bin1_id": np.minimum(i, j), "bin2_id": np.maximum(i, j)})
    px = px.drop_duplicates().sort_values(["bin1_id", "bin2_id"])
    px["count"] = rng.integers(1, 9, len(px))
    px["foo"] = rng.random(len(px))
    cooler.create_cooler(
        path, bins, px, columns=["count", "foo"], dtypes={"foo": float}, **kwargs
    )
    return path


def dense(clr, field):
    return clr.matrix(balance=False, field=field)[:]


def coarsen_dense(mat, bins, factor):
    """
    Coarsen a symmetric matrix stored as its upper triangle: sum blocks of
    factor x factor bins (within chromosomes) of the stored triangle, then
    fill in the lower triangle.
    """
    mat = np.triu(mat)
    chroms = bins["chrom"].astype(str).to_numpy()
    groups = []
    for chrom in pd.unique(chroms):
        idx = np.flatnonzero(chroms == chrom)
        groups += [idx[k : k + factor] for k in range(0, len(idx), factor)]
    out = np.zeros((len(groups), len(groups)))
    for a, ga in enumerate(groups):
        for b, gb in enumerate(groups):
            out[a, b] = mat[np.ix_(ga, gb)].sum()
    out = np.triu(out)
    return out + np.triu(out, 1).T


def test_coarsen_keeps_value_columns(tmp_path):
    src = make_cooler(str(tmp_path / "c.cool"), 10)
    clr = cooler.Cooler(src)
    out = str(tmp_path / "k.cool")
    dtypes = {"foo": np.float64}
    coarsen_cooler(src, out, 3, chunksize=50, columns=["count", "foo"], dtypes=dtypes)
    assert dtypes == {"foo": np.float64}  # not modified
    res = cooler.Cooler(out)
    assert list(res.pixels().columns) == ["bin1_id", "bin2_id", "count", "foo"]
    for field in ["count", "foo"]:
        expected = coarsen_dense(dense(clr, field), clr.bins()[:], 3)
        np.testing.assert_allclose(dense(res, field), expected)

    # Only a non-count column
    out = str(tmp_path / "f.cool")
    coarsen_cooler(src, out, 2, chunksize=50, columns=["foo"])
    assert list(cooler.Cooler(out).pixels().columns) == ["bin1_id", "bin2_id", "foo"]


def test_coarsen_keeps_assembly_and_metadata(tmp_path):
    src = make_cooler(
        str(tmp_path / "a.cool"), 10, assembly="hg38", metadata={"lab": "x"}
    )
    out = str(tmp_path / "a2.cool")
    coarsen_cooler(src, out, 2, chunksize=50)
    info = cooler.Cooler(out).info
    assert info["genome-assembly"] == "hg38"
    assert info["metadata"] == {"lab": "x"}


def test_coarsen_factor(tmp_path):
    src = make_cooler(str(tmp_path / "c.cool"), 10)
    with pytest.raises(ValueError, match="at least 2"):
        coarsen_cooler(src, str(tmp_path / "o.cool"), 1, chunksize=50)
    result = CliRunner().invoke(
        coarsen_cli, [src, "-k", "1", "-o", str(tmp_path / "o.cool")]
    )
    assert result.exit_code == 2  # usage error, not a crash


def test_zoomify_keeps_value_columns(tmp_path):
    src = make_cooler(str(tmp_path / "c.cool"), 10)
    out = str(tmp_path / "z.mcool")
    zoomify_cooler(src, out, [10, 20, 40], chunksize=50, columns=["count", "foo"])
    clr = cooler.Cooler(src)
    for res, factor in [(20, 2), (40, 4)]:
        level = cooler.Cooler(f"{out}::resolutions/{res}")
        expected = coarsen_dense(dense(clr, "foo"), clr.bins()[:], factor)
        np.testing.assert_allclose(dense(level, "foo"), expected)


def test_zoomify_multiple_base_coolers(tmp_path):
    b10 = make_cooler(str(tmp_path / "b10.cool"), 10, seed=1)
    b30 = make_cooler(str(tmp_path / "b30.cool"), 30, seed=2)
    out = str(tmp_path / "z.mcool")
    zoomify_cooler([b10, b30], out, [10, 20, 30, 60], chunksize=50)
    assert cooler.fileops.list_coolers(out) == [
        "/resolutions/10",
        "/resolutions/20",
        "/resolutions/30",
        "/resolutions/60",
    ]
    # Both bases are kept, and 60 is coarsened from the 30 base
    c30 = cooler.Cooler(b30)
    np.testing.assert_array_equal(
        dense(cooler.Cooler(out + "::resolutions/30"), "count"), dense(c30, "count")
    )
    np.testing.assert_array_equal(
        dense(cooler.Cooler(out + "::resolutions/60"), "count"),
        coarsen_dense(dense(c30, "count"), c30.bins()[:], 2),
    )

    with pytest.raises(ValueError, match="different resolutions"):
        zoomify_cooler([b10, b10], str(tmp_path / "d.mcool"), [10, 20], chunksize=50)


def test_zoomify_cli_binary_resolutions_and_balance(tmp_path):
    src = make_cooler(str(tmp_path / "c.cool"), 10)

    # Binary progression suffix (used to crash)
    result = CliRunner().invoke(
        zoomify_cli, [src, "-r", "10b", "-o", str(tmp_path / "b.mcool")]
    )
    assert result.exit_code == 0, result.output

    out = str(tmp_path / "z.mcool")
    result = CliRunner().invoke(
        zoomify_cli,
        [
            src,
            "-r",
            "20,40",
            "--balance",
            "--balance-args",
            "--name ice --ignore-diags 0 --mad-max 0 --min-nnz 0",
            "-o",
            out,
        ],
    )
    assert result.exit_code == 0, result.output
    levels = cooler.fileops.list_coolers(out)
    assert levels == ["/resolutions/10", "/resolutions/20", "/resolutions/40"]
    # Every level is balanced, including the base, under the requested name
    for level in levels:
        assert "ice" in cooler.Cooler(f"{out}::{level}").bins().columns


def test_multiplier_sequence_keeps_bases():
    from cooler._reduce import get_multiplier_sequence

    resn, pred, mult = get_multiplier_sequence([10, 20, 30, 60], {10, 30})
    assert list(resn) == [10, 20, 30, 60]
    assert list(pred) == [-1, 0, -1, 2]  # 30 is a base, not coarsened from 10
    assert list(mult) == [-1, 2, -1, 2]


def test_merge_breakpoints_respect_bufsize():
    partition, cum = merge_breakpoints([np.arange(101)], 10)
    assert np.diff(cum).tolist() == [10] * 10
    assert partition.tolist() == list(range(0, 101, 10))


def test_merger_closes_files(tmp_path):
    a = make_cooler(str(tmp_path / "a.cool"), 10, seed=1)
    b = make_cooler(str(tmp_path / "b.cool"), 10, seed=2)
    before = h5py.h5f.get_obj_count(h5py.h5f.OBJ_ALL, h5py.h5f.OBJ_FILE)
    merger = CoolerMerger([cooler.Cooler(a), cooler.Cooler(b)], mergebuf=50)
    total = sum(chunk["count"].sum() for chunk in merger)
    after = h5py.h5f.get_obj_count(h5py.h5f.OBJ_ALL, h5py.h5f.OBJ_FILE)
    assert after == before
    assert total == (cooler.Cooler(a).info["sum"] + cooler.Cooler(b).info["sum"])


def test_merge_keeps_callers_dtypes(tmp_path):
    path = str(tmp_path / "a.cool")
    bins = cooler.binnify(pd.Series({"chr1": 100}), 10)
    px = pd.DataFrame({"bin1_id": [0], "bin2_id": [1], "count": [1]})
    cooler.create_cooler(path, bins, px)
    dtypes = {}
    cooler.merge_coolers(str(tmp_path / "m.cool"), [path, path], 100, dtypes=dtypes)
    assert dtypes == {}
