"""
Tests for the commands that read coolers: info, dump and show.
"""
import gzip
import json
import shutil
from io import StringIO

import numpy as np
import pandas as pd
import pytest
from _common import DATA, cooler_cmp

import cooler
from cooler.create._writer import _add_columns

TOY = DATA / "toy.symm.upper.2.cool"
SQUARE = DATA / "toy.asymm.2.cool"


def read_tsv(result):
    return pd.read_csv(StringIO(result.stdout), sep="\t")


def with_custom_weights(tmp_path):
    path = tmp_path / "test.cool"
    shutil.copy(TOY, path)
    n_bins = cooler.Cooler(str(path)).info["nbins"]
    _add_columns(str(path), "bins", {"custom": np.full(n_bins, 2.0)}, replace=True)
    return path


def test_info(run_cli, tmp_path):
    info = json.loads(run_cli("info", TOY).stdout)
    assert info["nbins"] == cooler.Cooler(str(TOY)).info["nbins"]
    assert run_cli("info", TOY, "--field", "bin-type").stdout.strip() == "fixed"
    out = tmp_path / "info.json"
    result = run_cli("info", TOY, "--field", "doesnotexist", "-o", out, exit_code=2)
    assert result.stdout == "" and not out.exists()
    json.loads(run_cli("info", TOY, "--metadata").stdout)

    run_cli("info", TOY, "-o", out)
    assert json.loads(out.read_text()) == info


def test_dump(run_cli, tmp_path):
    clr = cooler.Cooler(str(TOY))
    n_bins = clr.info["nbins"]

    assert len(read_tsv(run_cli("dump", TOY, "-H"))) == clr.info["nnz"]
    chroms = read_tsv(run_cli("dump", TOY, "-H", "-t", "chroms", "--columns", "length"))
    assert list(chroms.columns) == ["length"]
    bins = read_tsv(run_cli("dump", TOY, "-H", "-t", "bins", "-c", "chrom,start"))
    assert list(bins.columns) == ["chrom", "start"] and len(bins) == n_bins

    # Ranges
    px = read_tsv(run_cli("dump", TOY, "-H", "-r", "chr1"))
    lo, hi = clr.extent("chr1")
    assert px["bin1_id"].between(lo, hi - 1).all()
    assert px["bin2_id"].between(lo, hi - 1).all()
    run_cli("dump", TOY, "-r", "chr1:0-16", "-r2", "chr1:10-25")
    run_cli("dump", TOY, "-r", "chr1:10-25", "-r2", "chr1:0-5")

    # Joining and annotating
    joined = read_tsv(run_cli("dump", TOY, "-H", "--join"))
    assert list(joined.columns[:6]) == [
        "chrom1", "start1", "end1", "chrom2", "start2", "end2"
    ]
    starts = read_tsv(run_cli("dump", TOY, "-H", "--join", "--one-based-starts"))
    assert (starts["start1"] == joined["start1"] + 1).all()
    run_cli("dump", TOY, "--annotate", "chrom", "--one-based-starts")
    result = run_cli("dump", TOY, "--annotate", "nosuchcolumn", exit_code=2)
    assert result.stdout == ""
    run_cli("dump", TOY, "-t", "chroms", "-c", "nosuchcolumn", exit_code=2)
    run_cli("dump", TOY, "-t", "bins", "-c", "nosuchcolumn", exit_code=2)
    run_cli("dump", TOY, "-r2", "chr2", exit_code=2)  # needs --range

    # Unbalanced file
    run_cli("dump", TOY, "-b", exit_code=1)

    # Output to a file, compressed by extension
    out = tmp_path / "pixels.tsv.gz"
    run_cli("dump", TOY, "-o", out)
    with gzip.open(out, "rt") as f:
        assert sum(1 for _ in f) == clr.info["nnz"]


def test_dump_roundtrip(run_cli, tmp_path):
    # Symmetric-upper storage
    bins = read_tsv(run_cli("dump", TOY, "-H", "-t", "bins"))
    pixels = read_tsv(run_cli("dump", TOY, "-H"))
    cooler.create_cooler(str(tmp_path / "out.cool"), bins, pixels)
    cooler_cmp(str(TOY), str(tmp_path / "out.cool"))

    # Filling in the lower triangle
    pixels2 = read_tsv(run_cli("dump", TOY, "--fill-lower", "-H"))
    assert len(pixels2) > len(pixels)
    upper = pixels2[pixels2["bin1_id"] <= pixels2["bin2_id"]].reset_index(drop=True)
    assert np.allclose(pixels, upper)
    trans_lower = read_tsv(run_cli("dump", TOY, "-H", "-r", "chr2", "-r2", "chr1"))
    assert len(trans_lower) == 0
    trans_lower = read_tsv(run_cli("dump", TOY, "-fH", "-r", "chr2", "-r2", "chr1"))
    assert len(trans_lower) > 0

    # Square storage: --fill-lower is a no-op
    bins = read_tsv(run_cli("dump", SQUARE, "-H", "-t", "bins"))
    pixels = read_tsv(run_cli("dump", SQUARE, "-H"))
    cooler.create_cooler(
        str(tmp_path / "square.cool"), bins, pixels, symmetric_upper=False
    )
    cooler_cmp(str(SQUARE), str(tmp_path / "square.cool"))
    assert np.allclose(pixels, read_tsv(run_cli("dump", SQUARE, "--fill-lower", "-H")))
    lower1 = read_tsv(run_cli("dump", SQUARE, "-H", "-r", "chr2", "-r2", "chr1"))
    lower2 = read_tsv(run_cli("dump", SQUARE, "-f", "-H", "-r", "chr2", "-r2", "chr1"))
    assert np.allclose(lower1, lower2)


def test_dump_weight_name(run_cli, tmp_path):
    path = with_custom_weights(tmp_path)
    # --weight-name implies --balanced
    pixels = read_tsv(run_cli("dump", path, "-H", "--weight-name", "custom"))
    assert np.allclose(pixels["balanced"], 4 * pixels["count"])
    run_cli("dump", path, "--weight-name", "doesnotexist", exit_code=1)


def test_show(run_cli, tmp_path):
    pytest.importorskip("matplotlib")
    out = tmp_path / "plot.png"
    run_cli("show", TOY, "chr1", "-o", out)
    assert out.stat().st_size > 0
    for scale in ["linear", "log2", "log10"]:
        run_cli("show", TOY, "chr1", "-r2", "chr2", "--scale", scale, "-o", out)

    path = with_custom_weights(tmp_path)
    run_cli("show", path, "chr1", "-o", out, "--weight-name", "custom")
    run_cli(
        "show", path, "chr1", "-o", out, "--weight-name", "doesnotexist",
        exit_code=1,
    )
