# from io import StringIO
import os.path as op
import tempfile
from glob import glob

import numpy as np
import pandas as pd
import pytest
import simplejson as json

# from _common import cooler_cmp
from click.testing import CliRunner
from pandas.api import types

import cooler

### INGEST AND AGGREGATION ###
from cooler.cli.cload import pairs as cload_pairs
from cooler.cli.load import load

tmp = tempfile.gettempdir()
testdir = op.realpath(op.dirname(__file__))
datadir = op.join(testdir, "data")


def _run_cload_pairs(runner, binsize, extra_args):
    args = [
        op.join(datadir, "toy.chrom.sizes") + ":" + str(binsize),
        op.join(datadir, "toy.pairs"),
        f"toy.{binsize}.cool",
        "-c1", "2",
        "-p1", "3",
        "-c2", "4",
        "-p2", "5",
        "--assembly", "toy",
        "--chunksize", "10",
    ] + extra_args
    return runner.invoke(cload_pairs, args)


def _cmp_pixels_2_bg(f_out, f_ref, one_based_ref=True):
    # output, 1-based starts
    out_df = cooler.Cooler(f_out).pixels(join=True)[:]
    if one_based_ref:
        out_df["start1"] += 1
        out_df["start2"] += 1

    # reference
    ref_df = pd.read_csv(
        f_ref,
        sep="\t",
        names=["chrom1", "start1", "end1", "chrom2", "start2", "end2", "count"],
    )

    assert np.all(out_df == ref_df)


# '--no-symmetric-upper'
# '--input-copy-status', 'unique|duplex',
@pytest.mark.parametrize(
    "ref,extra_args",
    [
        ("symm.upper", []),  # reflect triu pairs
        ("symm.upper", ["--input-copy-status", "unique"]),  # reflect triu pairs
        ("asymm", ["--no-symmetric-upper"]),
    ],
)
def test_cload_symm_asymm(ref, extra_args):
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = _run_cload_pairs(runner, 2, extra_args)
        assert result.exit_code == 0
        _cmp_pixels_2_bg("toy.2.cool", op.join(datadir, f"toy.{ref}.2.bg2"))


# '--temp-dir', '',
# '--no-delete-temp',
# '--max-merge', '',
@pytest.mark.parametrize(
    "ref,extra_args", [("symm.upper", ["--temp-dir", ".", "--no-delete-temp"])]
)
def test_cload_mergepass(ref, extra_args):
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = _run_cload_pairs(runner, 2, extra_args)
        assert result.exit_code == 0
        _cmp_pixels_2_bg("toy.2.cool", op.join(datadir, f"toy.{ref}.2.bg2"))
        assert len(cooler.fileops.list_coolers(glob("*.cool")[0])) > 0


# '--field', '',
# '--no-count', '',
def test_cload_field():
    runner = CliRunner()
    with runner.isolated_filesystem():
        extra_args = ["--field", "score=8"]
        result = _run_cload_pairs(runner, 2, extra_args)
        assert result.exit_code == 0
        pixels = cooler.Cooler("toy.2.cool").pixels()[:]
        assert "count" in pixels.columns and types.is_integer_dtype(
            pixels.dtypes["count"]
        )
        assert "score" in pixels.columns and types.is_float_dtype(
            pixels.dtypes["score"]
        )

        extra_args = ["--field", "count=8"]
        result = _run_cload_pairs(runner, 2, extra_args)
        assert result.exit_code == 0
        pixels = cooler.Cooler("toy.2.cool").pixels()[:]
        assert "count" in pixels.columns and types.is_integer_dtype(
            pixels.dtypes["count"]
        )
        assert np.allclose(pixels["count"][:], 0)

        extra_args = ["--field", "count=8:dtype=float"]
        result = _run_cload_pairs(runner, 2, extra_args)
        assert result.exit_code == 0
        pixels = cooler.Cooler("toy.2.cool").pixels()[:]
        assert "count" in pixels.columns and types.is_float_dtype(
            pixels.dtypes["count"]
        )
        assert np.allclose(pixels["count"][:], 0.2)

        extra_args = ["--field", "count=8:agg=min,dtype=float"]
        result = _run_cload_pairs(runner, 2, extra_args)
        assert result.exit_code == 0
        pixels = cooler.Cooler("toy.2.cool").pixels()[:]
        assert "count" in pixels.columns and types.is_float_dtype(
            pixels.dtypes["count"]
        )
        assert np.allclose(pixels["count"][:], 0.1)

        ## don't implement the --no-count for now
        # extra_args =  ['--field', 'score=7:dtype=float', '--no-count']
        # result = _run_cload_pairs(runner, 2, extra_args)
        # assert result.exit_code == 0
        # pixels = cooler.Cooler('toy.2.cool').pixels()[:]
        # assert 'count' not in pixels.columns
        # assert 'score' in pixels.columns and types.is_float_dtype(pixels.dtypes['score'])


# '--metadata', '',
# '--zero-based',
# '--comment-char', '',
# '--storage-options', '',
def test_cload_other_options():
    runner = CliRunner()
    with runner.isolated_filesystem():
        meta = {"foo": "bar", "number": 42}
        with open("meta.json", "w") as f:
            json.dump(meta, f)
        extra_args = [
            "--metadata",
            "meta.json",
            "--zero-based",
            "--storage-options",
            "shuffle=True,fletcher32=True,compression=lzf",
        ]
        result = _run_cload_pairs(runner, 2, extra_args)
        assert result.exit_code == 0
        c = cooler.Cooler("toy.2.cool")
        assert c.info["metadata"] == meta
        with c.open("r") as h5:
            dset = h5["bins/start"]
            assert dset.shuffle
            assert dset.fletcher32
            assert dset.compression == "lzf"


def _run_load(runner, matrix_file, format, binsize, extra_args):
    args = [
        "-f",
        format,
        op.join(datadir, "toy.chrom.sizes") + ":" + str(binsize),
        op.join(datadir, matrix_file),
        f"toy.{binsize}.cool",
        "--assembly",
        "toy",
        "--chunksize",
        "10",
    ] + extra_args
    return runner.invoke(load, args)


# '--no-symmetric-upper'
# '--input-copy-status', 'unique|duplex',
@pytest.mark.parametrize(
    "ref,extra_args",
    [
        ("symm.upper", []),  # reflect tril pairs
        ("symm.upper", ["--one-based", "--input-copy-status", "unique"]),  # reflect tril pairs
        ("asymm", ["--one-based", "--no-symmetric-upper"]),
    ],
)
def test_load_symm_asymm(ref, extra_args):
    runner = CliRunner()
    with runner.isolated_filesystem():
        ref = op.join(datadir, f"toy.{ref}.2.bg2")
        result = _run_load(runner, ref, "bg2", 2, extra_args)
        assert result.exit_code == 0
        _cmp_pixels_2_bg("toy.2.cool", ref)


# '--field', '',
def test_load_field():
    runner = CliRunner()
    with runner.isolated_filesystem():
        extra_args = ["--field", "count=7:dtype=float"]
        result = _run_load(runner, "toy.symm.upper.2.bg2", "bg2", 2, extra_args)
        assert result.exit_code == 0
        pixels1 = cooler.Cooler(op.join(datadir, "toy.symm.upper.2.cool")).pixels()[:]
        pixels2 = cooler.Cooler("toy.2.cool").pixels()[:]
        assert "count" in pixels2.columns and types.is_float_dtype(
            pixels2.dtypes["count"]
        )
        assert np.allclose(pixels1["count"][:], pixels2["count"][:])


def test_load_field2():
    runner = CliRunner()
    with runner.isolated_filesystem():
        extra_args = ["--count-as-float"]
        result = _run_load(runner, "toy.symm.upper.2.bg2", "bg2", 2, extra_args)
        assert result.exit_code == 0
        pixels1 = cooler.Cooler(op.join(datadir, "toy.symm.upper.2.cool")).pixels()[:]
        pixels2 = cooler.Cooler("toy.2.cool").pixels()[:]
        assert "count" in pixels2.columns and types.is_float_dtype(
            pixels2.dtypes["count"]
        )
        assert np.allclose(pixels1["count"][:], pixels2["count"][:])
