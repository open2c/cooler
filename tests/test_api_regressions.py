"""
Regression tests for bugs in the API, selectors and table operations.
"""
import copy
import pathlib
import pickle

import h5py
import numpy as np
import pandas as pd
import pytest

import cooler
from cooler.core import get, put
from cooler.parallel import split
from cooler.util import parse_humanized, parse_region_string

CHROMSIZES = pd.Series({"chr1": 100, "chr2": 60})


def make_cooler(path, bins=None, **kwargs):
    if bins is None:
        bins = cooler.binnify(CHROMSIZES, 10)
    rng = np.random.default_rng(0)
    n = len(bins)
    i, j = rng.integers(0, n, 100), rng.integers(0, n, 100)
    px = pd.DataFrame({"bin1_id": np.minimum(i, j), "bin2_id": np.maximum(i, j)})
    px = px.drop_duplicates().sort_values(["bin1_id", "bin2_id"])
    px["count"] = rng.integers(1, 9, len(px))
    cooler.create_cooler(path, bins, px, **kwargs)
    return path


def test_column_names_containing_chrom(tmp_path):
    bins = cooler.binnify(CHROMSIZES, 10)
    bins["chrom_arm"] = np.arange(len(bins)) % 2
    bins["name"] = [f"bin{i}" for i in range(len(bins))]  # used to crash
    clr = cooler.Cooler(make_cooler(str(tmp_path / "c.cool"), bins))
    np.testing.assert_array_equal(clr.bins()["chrom_arm"][:], bins["chrom_arm"])
    assert clr.bins()["name"][:3].tolist() == ["bin0", "bin1", "bin2"]
    assert clr.bins()["chrom"][:1].tolist() == ["chr1"]


def test_parse_humanized_decimals():
    values = range(1, 100_000)  # 0.01 kb to 999.99 kb
    wrong = [x for x in values if parse_humanized(f"{x / 100:.2f}kb") != x * 10]
    assert wrong == []
    assert parse_humanized("4.1Mb") == 4_100_000
    with pytest.raises(ValueError, match="whole number"):
        parse_humanized("1.0001kb")


@pytest.mark.parametrize("region", ["chr1:100-200-300", "chr1:100-200 x", "a:b:1-2"])
def test_parse_region_string_rejects_junk(region):
    with pytest.raises(ValueError):
        parse_region_string(region)


@pytest.mark.parametrize(
    "bins, match",
    [
        ([("chr1", 5, 10), ("chr1", 10, 20)], "does not start at 0"),
        ([("chr1", 0, 10), ("chr1", 15, 20)], "gaps or overlapping"),
        ([("chr1", 0, 10), ("chr1", 5, 20)], "gaps or overlapping"),
        ([("chr1", 0, 10), ("chr1", 10, 10)], "end is not after"),
        ([("chr1", 0, 10), ("chr2", 0, 10), ("chr1", 10, 20)], "not consecutive"),
    ],
)
def test_bin_table_must_tile(tmp_path, bins, match):
    # #254
    bins = pd.DataFrame(bins, columns=["chrom", "start", "end"])
    px = pd.DataFrame({"bin1_id": [0], "bin2_id": [1], "count": [1]})
    with pytest.raises(ValueError, match=match):
        cooler.create_cooler(str(tmp_path / "c.cool"), bins, px)


def test_variable_bins_region_before_first_bin(tmp_path):
    # Files written before bin tables were checked can have a first bin that
    # doesn't start at 0. A query must not reach the previous chromosome.
    bins = pd.DataFrame({
        "chrom": ["chr1"] * 3 + ["chr2"] * 2,
        "start": [0, 30, 70, 0, 25],
        "end": [30, 70, 100, 25, 60],
    })
    path = make_cooler(str(tmp_path / "v.cool"), bins)
    with h5py.File(path, "r+") as f:
        f["bins/start"][3] = 10
    clr = cooler.Cooler(path)
    assert clr.extent("chr2:0-20") == (3, 4)
    assert clr.bins().fetch("chr2:0-20")["chrom"].tolist() == ["chr2"]


def test_out_of_range_indices(tmp_path):
    clr = cooler.Cooler(make_cooler(str(tmp_path / "c.cool")))
    cooler.balance_cooler(clr, store=True, mad_max=0, min_nnz=0, ignore_diags=0)
    n = len(clr.bins())
    with pytest.raises(IndexError):
        clr.bins()[-n - 1]
    with pytest.raises(IndexError):
        clr.matrix()[-n - 1]
    with pytest.raises(IndexError):
        clr.pixels()[len(clr.pixels())]
    assert clr.matrix()[n - 2 : n + 5, :3].shape == (2, 3)
    assert clr.matrix(balance=True)[n - 2 : n + 5].shape == (2, n)
    assert clr.matrix()[10:5].shape == (0, n)
    assert clr.bins()[10:5].shape == (0, 3 + 1)
    with pytest.raises(TypeError):
        clr.bins()[1.5]
    with pytest.raises(TypeError):
        clr.bins()[[0, 1]]


def test_single_column_selector(tmp_path):
    clr = cooler.Cooler(make_cooler(str(tmp_path / "c.cool")))
    sel = clr.bins()["start"]
    assert sel.columns.tolist() == ["start"]
    assert sel.keys() == ["start"]
    assert "start" in sel
    assert list(sel.dtypes.index) == ["start"]


def test_put_appends(tmp_path):
    with h5py.File(tmp_path / "t.h5", "w") as f:
        grp = f.create_group("t")
        put(grp, pd.DataFrame({
            "cat": pd.Categorical(["a", "b", "a"]), "s": ["x", "yy", "z"]
        }))
        more = pd.DataFrame({"cat": pd.Categorical(["b", "b"]), "s": ["p", "q"]})
        put(grp, more, lo=3)
        out = get(grp)
        assert out["cat"].tolist() == ["a", "b", "a", "b", "b"]
        assert out["s"].tolist() == ["x", "yy", "z", "p", "q"]
        with pytest.raises(ValueError, match="longer"):
            put(grp, pd.DataFrame({"s": ["toolong"]}), lo=5)
        with pytest.raises(ValueError, match="Categories"):
            put(grp, pd.DataFrame({"cat": pd.Categorical(["c"])}), lo=5)


def test_strings_stored_as_bytes(tmp_path):
    # Files from other writers: bytes attributes, variable-length strings
    path = make_cooler(str(tmp_path / "c.cool"))
    expected = cooler.Cooler(path).matrix(balance=False)[:]
    with h5py.File(path, "r+") as f:
        for key in ["format", "storage-mode", "bin-type"]:
            f.attrs[key] = np.bytes_(f.attrs[key])
        names = f["chroms/name"][:].astype(str).astype(object)
        del f["chroms/name"]
        f.create_dataset("chroms/name", data=names, dtype=h5py.string_dtype())
    assert cooler.fileops.is_cooler(path)
    clr = cooler.Cooler(path)
    assert clr.storage_mode == "symmetric-upper"
    assert clr.chromnames == ["chr1", "chr2"]
    np.testing.assert_array_equal(clr.matrix(balance=False)[:], expected)
    assert clr.matrix(balance=False).fetch("chr2").shape == (6, 6)


def test_path_arguments(tmp_path):
    path = tmp_path / "c.cool"
    make_cooler(path)
    clr = cooler.Cooler(path)
    assert clr.info["nbins"] == 16
    with cooler.util.open_hdf5(path) as f:
        assert "pixels" in f


def test_annotate_uint64_bin_ids(tmp_path):
    clr = cooler.Cooler(make_cooler(str(tmp_path / "c.cool")))
    px = clr.pixels()[:5]
    px[["bin1_id", "bin2_id"]] = px[["bin1_id", "bin2_id"]].astype(np.uint64)
    out = cooler.annotate(px, clr.bins()[:])
    assert out["chrom1"].tolist() == clr.pixels(join=True)[:5]["chrom1"].tolist()


def test_pipe_can_be_pickled_and_copied(tmp_path):
    clr = cooler.Cooler(make_cooler(str(tmp_path / "c.cool")))
    pipe = split(clr, chunksize=10).pipe(len)
    for other in [copy.deepcopy(pipe), pickle.loads(pickle.dumps(pipe))]:
        assert other.map is map
        assert list(other) == list(pipe)


def test_coarsest_zoom_level_has_a_binsize(tmp_path):
    src = make_cooler(str(tmp_path / "c.cool"))
    out = str(tmp_path / "z.mcool")
    # At 200 bp, every chromosome is a single bin
    cooler.zoomify_cooler(src, out, [10, 20, 200], chunksize=50)
    clr = cooler.Cooler(out + "::resolutions/200")
    assert clr.info["bin-type"] == "fixed"
    assert clr.binsize == 200

    bins = cooler.binnify(CHROMSIZES, 10)
    px = pd.DataFrame({"bin1_id": [0], "bin2_id": [1], "count": [1]})
    with pytest.raises(ValueError, match="uniform size"):
        cooler.create.create(str(tmp_path / "x.cool"), bins, px, binsize=20)


def test_path_object_uri(tmp_path):
    # PathLike URIs work wherever a URI is parsed
    assert cooler.util.parse_cooler_uri(pathlib.Path("a.mcool")) == ("a.mcool", "/")
