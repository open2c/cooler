"""
Tests for the commands that create coolers from contacts or pixels: cload
pairs and load.
"""
import json

import numpy as np
import pandas as pd
import pytest
from _common import DATA
from pandas.api import types

import cooler

CHROMSIZES = DATA / "toy.chrom.sizes"


def cload_pairs(run_cli, out, *args, binsize=2, exit_code=0):
    return run_cli(
        "cload", "pairs", f"{CHROMSIZES}:{binsize}", DATA / "toy.pairs", out,
        "-c1", 2, "-p1", 3, "-c2", 4, "-p2", 5,
        "--assembly", "toy", "--chunksize", 10, *args,
        exit_code=exit_code,
    )


def load(run_cli, out, matrix_file, format, *args, binsize=2):
    return run_cli(
        "load", "-f", format, f"{CHROMSIZES}:{binsize}", matrix_file, out,
        "--assembly", "toy", "--chunksize", 10, *args,
    )


def assert_pixels_match_bg2(path, ref, one_based_ref=True):
    out = cooler.Cooler(str(path)).pixels(join=True)[:]
    if one_based_ref:
        out["start1"] += 1
        out["start2"] += 1
    expected = pd.read_csv(
        ref,
        sep="\t",
        names=["chrom1", "start1", "end1", "chrom2", "start2", "end2", "count"],
    )
    assert np.all(out == expected)


@pytest.mark.parametrize(
    "ref,args",
    [
        ("symm.upper", []),  # reflect pairs in the lower triangle
        ("symm.upper", ["--input-copy-status", "unique"]),
        ("asymm", ["--no-symmetric-upper"]),
    ],
)
def test_cload_symm_asymm(run_cli, tmp_path, ref, args):
    out = tmp_path / "toy.2.cool"
    cload_pairs(run_cli, out, *args)
    assert_pixels_match_bg2(out, DATA / f"toy.{ref}.2.bg2")
    assert cooler.Cooler(str(out)).info["genome-assembly"] == "toy"


def test_cload_mergepass(run_cli, tmp_path):
    out = tmp_path / "toy.2.cool"
    cload_pairs(run_cli, out, "--temp-dir", tmp_path, "--no-delete-temp")
    assert_pixels_match_bg2(out, DATA / "toy.symm.upper.2.bg2")
    temp_files = [p for p in tmp_path.glob("*.cool") if p != out]
    assert temp_files and cooler.fileops.list_coolers(temp_files[0])


def test_cload_field(run_cli, tmp_path):
    out = tmp_path / "toy.2.cool"

    cload_pairs(run_cli, out, "--field", "score=8")
    pixels = cooler.Cooler(str(out)).pixels()[:]
    assert types.is_integer_dtype(pixels.dtypes["count"])
    assert types.is_float_dtype(pixels.dtypes["score"])

    # Fields named count replace the record counts
    cload_pairs(run_cli, out, "--field", "count=8")
    pixels = cooler.Cooler(str(out)).pixels()[:]
    assert types.is_integer_dtype(pixels.dtypes["count"])
    assert np.allclose(pixels["count"], 0)

    cload_pairs(run_cli, out, "--field", "count=8:dtype=float")
    pixels = cooler.Cooler(str(out)).pixels()[:]
    assert types.is_float_dtype(pixels.dtypes["count"])
    assert np.allclose(pixels["count"], 0.2)

    cload_pairs(run_cli, out, "--field", "count=8:agg=min,dtype=float")
    assert np.allclose(cooler.Cooler(str(out)).pixels()["count"][:], 0.1)

    # Standard columns' dtypes can be set without a field number
    cload_pairs(run_cli, out, "--field", "bin1_id:dtype=int32")
    assert cooler.Cooler(str(out)).pixels()["bin1_id"][:].dtype == np.int32
    cload_pairs(run_cli, out, "--field", "score", exit_code=2)
    cload_pairs(run_cli, out, "--field", "score=8:agg=mean", exit_code=2)


def test_cload_other_options(run_cli, tmp_path):
    meta = {"foo": "bar", "number": 42}
    (tmp_path / "meta.json").write_text(json.dumps(meta))
    out = tmp_path / "toy.2.cool"
    cload_pairs(
        run_cli, out,
        "--metadata", tmp_path / "meta.json",
        "--storage-options", "shuffle=True,fletcher32=True,compression=lzf",
    )
    clr = cooler.Cooler(str(out))
    assert clr.info["metadata"] == meta
    with clr.open("r") as h5:
        dset = h5["bins/start"]
        assert dset.shuffle
        assert dset.fletcher32
        assert dset.compression == "lzf"


@pytest.mark.parametrize(
    "ref,args",
    [
        ("symm.upper", []),  # reflect pixels in the lower triangle
        ("symm.upper", ["--one-based", "--input-copy-status", "unique"]),
        ("asymm", ["--one-based", "--no-symmetric-upper"]),
    ],
)
def test_load_symm_asymm(run_cli, tmp_path, ref, args):
    ref = DATA / f"toy.{ref}.2.bg2"
    out = tmp_path / "toy.2.cool"
    load(run_cli, out, ref, "bg2", *args)
    assert_pixels_match_bg2(out, ref)


@pytest.mark.parametrize(
    "args", [["--field", "count=7:dtype=float"], ["--count-as-float"]]
)
def test_load_float_counts(run_cli, tmp_path, args):
    out = tmp_path / "toy.2.cool"
    load(run_cli, out, DATA / "toy.symm.upper.2.bg2", "bg2", *args)
    expected = cooler.Cooler(str(DATA / "toy.symm.upper.2.cool")).pixels()[:]
    pixels = cooler.Cooler(str(out)).pixels()[:]
    assert types.is_float_dtype(pixels.dtypes["count"])
    assert np.allclose(expected["count"], pixels["count"])


def test_load_other_options(run_cli, tmp_path):
    meta = {"foo": "bar"}
    (tmp_path / "meta.json").write_text(json.dumps(meta))
    out = tmp_path / "toy.2.cool"
    load(
        run_cli, out, DATA / "toy.symm.upper.2.bg2", "bg2",
        "--metadata", tmp_path / "meta.json",
        "--field", "bin1_id:dtype=int32",
        "--field", "count:dtype=float32",
        "--storage-options", "compression=lzf",
    )
    clr = cooler.Cooler(str(out))
    assert clr.info["metadata"] == meta
    with clr.open("r") as h5:
        assert h5["pixels/bin1_id"].dtype == np.int32
        assert h5["pixels/count"].dtype == np.float32
        assert h5["pixels/count"].compression == "lzf"

    # The input copy status of a matrix with both triangles
    load(
        run_cli, out, DATA / "toy.symm.upper.2.bg2", "bg2",
        "--input-copy-status", "duplex",
    )


def test_load_field_from_count_column(run_cli, tmp_path):
    # A field can read the same column as another, e.g. a float copy of count
    out = tmp_path / "toy.1.cool"
    run_cli(
        "load", "-f", "coo", "--field", "v=3:dtype=float64",
        f"{CHROMSIZES}:1", DATA / "toy.symm.upper.1.zb.coo", out,
    )
    pixels = cooler.Cooler(str(out)).pixels()[:]
    assert types.is_float_dtype(pixels.dtypes["v"])
    assert (pixels["v"] == pixels["count"]).all()
