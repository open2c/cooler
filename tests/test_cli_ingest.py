"""
Tests for the commands that create coolers from contacts or pixels: cload
pairs and load.
"""
import json

import h5py
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


def test_cload_pairs_options(run_cli, tmp_path):
    expected = tmp_path / "expected.cool"
    cload_pairs(run_cli, expected)

    # Pairs from stdin
    out = tmp_path / "stdin.cool"
    run_cli(
        "cload", "pairs", f"{CHROMSIZES}:2", "-", out,
        "-c1", 2, "-p1", 3, "-c2", 4, "-p2", 5,
        input=(DATA / "toy.pairs").read_text(),
    )
    pd.testing.assert_frame_equal(
        cooler.Cooler(str(out)).pixels()[:], cooler.Cooler(str(expected)).pixels()[:]
    )

    # A storage option with a list value
    cload_pairs(run_cli, out, "--storage-options", "chunks=[1]")
    with h5py.File(out, "r") as f:
        assert f["pixels/count"].chunks == (1,)

    # Field numbers start at 1
    cload_pairs(run_cli, out, "-c1", 0, exit_code=2)


def write_hiclib_pairs(path):
    # A hiclib "frag" file: chromosome IDs and zero-based positions of both
    # sides, upper triangular and sorted by the first side
    records = pd.read_csv(DATA / "toy.pairs", sep="\t", header=None)
    chrom_ids = {"chr1": 0, "chr2": 1}
    c1, p1 = records[1].map(chrom_ids), records[2] - 1
    c2, p2 = records[3].map(chrom_ids), records[4] - 1
    flip = (c1 > c2) | ((c1 == c2) & (p1 > p2))
    df = pd.DataFrame({
        "chrms1": np.where(flip, c2, c1), "cuts1": np.where(flip, p2, p1),
        "chrms2": np.where(flip, c1, c2), "cuts2": np.where(flip, p1, p2),
    }).sort_values(["chrms1", "cuts1"])
    with h5py.File(path, "w") as f:
        for name in df.columns:
            f[name] = df[name].to_numpy()


def test_cload_hiclib(run_cli, tmp_path):
    frag = tmp_path / "frag.h5"
    write_hiclib_pairs(frag)
    (tmp_path / "meta.json").write_text(json.dumps({"foo": "bar"}))
    out = tmp_path / "hiclib.cool"
    run_cli(
        "cload", "hiclib", f"{CHROMSIZES}:2", frag, out,
        "--metadata", tmp_path / "meta.json", "--chunksize", 7,
    )
    expected = tmp_path / "pairs.cool"
    cload_pairs(run_cli, expected)
    pd.testing.assert_frame_equal(
        cooler.Cooler(str(out)).pixels()[:], cooler.Cooler(str(expected)).pixels()[:]
    )
    assert cooler.Cooler(str(out)).info["metadata"] == {"foo": "bar"}


def test_cload_tabix_metadata(run_cli, tmp_path):
    pytest.importorskip("pysam")
    (tmp_path / "meta.json").write_text(json.dumps({"foo": "bar"}))
    out = tmp_path / "tabix.cool"
    run_cli(
        "cload", "tabix", DATA / "hg19.bins.2000kb.bed.gz",
        DATA / "hg19.GM12878-MboI.pairs.subsample.sorted.txt.gz", out,
        "--metadata", tmp_path / "meta.json",
    )
    assert cooler.Cooler(str(out)).info["metadata"] == {"foo": "bar"}


def test_get_header():
    import io

    from cooler.cli.cload import get_header

    # A binary stream with peek() but no text wrapper
    stream = io.BufferedReader(io.BytesIO(b"#a\n#b\n1\t2\n"))
    header, rest = get_header(stream)
    assert header == ["#a", "#b"]
    assert rest.read() == b"1\t2\n"

    with pytest.raises(ValueError, match="comment char"):
        get_header(stream, comment_char="")
    with pytest.raises(ValueError, match="peek"):
        get_header(io.StringIO("#a\n"))


def test_load_more_options(run_cli, tmp_path):
    bg2 = DATA / "toy.symm.upper.2.bg2"
    out = tmp_path / "toy.2.cool"
    run_cli(
        "load", "-f", "bg2", f"{CHROMSIZES}:2", "-", out,
        "--storage-options", "chunks=[1]", input=bg2.read_text(),
    )
    assert_pixels_match_bg2(out, bg2)
    with h5py.File(out, "r") as f:
        assert f["pixels/count"].chunks == (1,)
    # A value field needs a field number
    run_cli(
        "load", "-f", "bg2", f"{CHROMSIZES}:2", bg2, out, "--field", "score",
        exit_code=2,
    )
