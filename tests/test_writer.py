"""Tests of the chunked cooler writer underneath the create functions."""

import threading

import h5py
import numpy as np
import pandas as pd
import pytest

import cooler
from cooler.create import BadInputError
from cooler.create._writer import _CoolerWriter

BINS = cooler.util.binnify(pd.Series({"chr1": 60, "chr2": 40}), 10)  # 10 bins


def pixels(pairs, count=1):
    bin1, bin2 = zip(*pairs) if pairs else ((), ())
    return pd.DataFrame(
        {
            "bin1_id": np.array(bin1, dtype=int),
            "bin2_id": np.array(bin2, dtype=int),
            "count": np.full(len(pairs), count, dtype=int),
        }
    )


def tmp_groups(path):
    with h5py.File(path, "r") as f:
        found = []
        f.visit(lambda name: found.append(name) if "__writing__" in name else None)
        return found


def test_chunks_match_single_write(tmp_path):
    chunks = [pixels([(0, 0), (0, 5)]), pixels([]), pixels([(1, 1), (3, 9)], 2)]
    path = str(tmp_path / "w.cool")
    with _CoolerWriter(path, BINS) as w:
        for chunk in chunks:
            w.write(chunk)

    ref = str(tmp_path / "ref.cool")
    cooler.create_cooler(ref, BINS, pd.concat(chunks))
    a, b = cooler.Cooler(path), cooler.Cooler(ref)
    pd.testing.assert_frame_equal(a.pixels()[:], b.pixels()[:])
    assert a.info["nnz"] == 4
    assert a.info["sum"] == 6
    assert tmp_groups(path) == []


def test_pixel_datasets_have_nnz_rows(tmp_path):
    path = str(tmp_path / "e.cool")
    with _CoolerWriter(path, BINS):
        pass
    with h5py.File(path, "r") as f:
        assert all(len(f["pixels"][col]) == 0 for col in f["pixels"])
    assert cooler.Cooler(path).info["nnz"] == 0


@pytest.mark.parametrize(
    "first, second, problem",
    [
        ([(2, 3), (2, 5)], [(1, 9)], "out-of-order"),
        ([(2, 3), (2, 5)], [(2, 5), (4, 4)], "duplicate"),
        ([(2, 5), (2, 3)], [], "out-of-order"),
    ],
)
def test_order_is_checked_across_chunks(tmp_path, first, second, problem):
    path = str(tmp_path / "o.cool")
    with pytest.raises(BadInputError, match=problem), _CoolerWriter(path, BINS) as w:
        w.write(pixels(first))
        w.write(pixels(second))


def test_duplicates_allowed_without_dupcheck(tmp_path):
    path = str(tmp_path / "d.cool")
    with _CoolerWriter(path, BINS, dupcheck=False) as w:
        w.write(pixels([(2, 3), (2, 5)]))
        w.write(pixels([(2, 5)]))
    assert cooler.Cooler(path).info["nnz"] == 3


def test_ensure_sorted_sorts_within_chunks(tmp_path):
    path = str(tmp_path / "s.cool")
    with _CoolerWriter(path, BINS, ensure_sorted=True) as w:
        w.write(pixels([(0, 5), (0, 0)]))
        w.write(pixels([(3, 9), (1, 1)]))
    px = cooler.Cooler(path).pixels()[:]
    assert px[["bin1_id", "bin2_id"]].values.tolist() == [
        [0, 0],
        [0, 5],
        [1, 1],
        [3, 9],
    ]


@pytest.mark.parametrize("group", ["/", "/resolutions/10"])
def test_failure_leaves_destination_untouched(tmp_path, group):
    path = str(tmp_path / "f.h5")
    uri = path + "::" + group
    cooler.create_cooler(uri, BINS, pixels([(0, 0)]), mode="w")
    with h5py.File(path, "r+") as f:
        f.create_group("other")

    with pytest.raises(BadInputError), _CoolerWriter(uri, BINS, mode="a") as w:
        w.write(pixels([(1, 1), (2, 2)]))
        w.write(pixels([(0, 0)]))  # out of order

    # The previous cooler and unrelated groups are intact
    assert cooler.Cooler(uri).pixels()[:][["bin1_id", "bin2_id"]].values.tolist() == [
        [0, 0]
    ]
    with h5py.File(path, "r") as f:
        assert "other" in f
    assert tmp_groups(path) == []

    # A successful write replaces only the destination group
    with _CoolerWriter(uri, BINS, mode="a") as w:
        w.write(pixels([(1, 1)]))
    assert cooler.Cooler(uri).info["nnz"] == 1
    with h5py.File(path, "r") as f:
        assert "other" in f


def test_leftover_temporary_group_is_replaced(tmp_path):
    path = str(tmp_path / "l.h5")
    with h5py.File(path, "w") as f:
        f.create_group("x" + ".__writing__").attrs["stale"] = True
    with _CoolerWriter(path + "::/x", BINS, mode="a") as w:
        w.write(pixels([(0, 1)]))
    assert cooler.Cooler(path + "::/x").info["nnz"] == 1
    assert tmp_groups(path) == []


def test_overflow_leaves_nothing_behind(tmp_path):
    path = str(tmp_path / "v.h5")
    with h5py.File(path, "w"):
        pass
    big = pixels([(0, 0)])
    big["count"] = np.array([3_000_000_000], dtype=np.int64)
    writer = _CoolerWriter(path + "::/c", BINS, mode="a")
    with pytest.raises(ValueError, match="int64"), writer as w:
        w.write(big)
    with h5py.File(path, "r") as f:
        assert list(f) == []


def test_lock_mode(tmp_path):
    path = str(tmp_path / "k.cool")
    lock = threading.Lock()
    with _CoolerWriter(path, BINS, lock=lock) as w:
        assert w._file is None  # the file is only opened under the lock
        w.write(pixels([(0, 0), (1, 2)]))
        assert not lock.locked()
    assert cooler.Cooler(path).info["nnz"] == 2


def test_create_rejects_unknown_arguments(tmp_path):
    path = str(tmp_path / "u.cool")
    with pytest.raises(TypeError, match="ordred"):
        cooler.create.create(path, BINS, pixels([(0, 0)]), ordred=True)
    with pytest.warns(FutureWarning):
        cooler.create.create(path, BINS, pixels([(0, 0)]), dtype={"count": float})
    assert cooler.Cooler(path).pixels()[:]["count"].dtype.kind == "f"


def test_unordered_rejects_unknown_arguments_up_front(tmp_path):
    from cooler.create._writer import _UnorderedCoolerWriter

    with pytest.raises(TypeError, match="ordred"):
        _UnorderedCoolerWriter(str(tmp_path / "x.cool"), BINS, ordred=True)


def test_unordered_cleans_up_temporary_files(tmp_path):
    from cooler.create._writer import _UnorderedCoolerWriter

    out = tmp_path / "out"
    out.mkdir()
    path = str(out / "u.cool")
    with _UnorderedCoolerWriter(path, BINS) as w:
        w.write(pixels([(3, 9), (0, 0)]))
        w.write(pixels([(0, 0), (1, 1)]))
    assert sorted(p.name for p in out.iterdir()) == ["u.cool"]
    px = cooler.Cooler(path).pixels()[:]
    assert px.values.tolist() == [[0, 0, 2], [1, 1, 1], [3, 9, 1]]

    # On error, nothing is written and temporary files are removed
    path = str(out / "bad.cool")
    with pytest.raises(BadInputError), _UnorderedCoolerWriter(path, BINS) as w:
        w.write(pixels([(0, 0)]))
        w.write(pixels([(0, 99)]))  # out of bounds
    assert sorted(p.name for p in out.iterdir()) == ["u.cool"]


def _addr(f, path):
    return h5py.h5o.get_info(f[path].id).addr


def test_scool_writer_adds_cells_incrementally(tmp_path):
    from cooler.create._writer import _ScoolWriter

    path = str(tmp_path / "s.scool")
    with _ScoolWriter(path, BINS) as w:
        # Sorted chunks
        w.add_cell("a", iter([pixels([(0, 0), (0, 3)]), pixels([(2, 2)])]))
        # Unsorted chunks, with an external sort
        w.add_cell(
            "b", iter([pixels([(4, 5), (0, 1)]), pixels([(0, 1)])]), ordered=False
        )
        # Per-cell bin columns
        w.add_cell("c", pixels([(1, 1)]), bins=BINS.assign(weight=np.arange(10.0)))

    assert cooler.fileops.is_scool_file(path)
    assert cooler.fileops.list_scool_cells(path) == ["/cells/a", "/cells/b", "/cells/c"]
    b = cooler.Cooler(path + "::/cells/b").pixels()[:]
    assert b.values.tolist() == [[0, 1, 2], [4, 5, 1]]
    c = cooler.Cooler(path + "::/cells/c")
    assert c.bins()["weight"][:].tolist() == list(np.arange(10.0))
    with h5py.File(path, "r") as f:
        assert f.attrs["ncells"] == 3
        for cell in "abc":
            assert _addr(f, f"cells/{cell}/bins/start") == _addr(f, "bins/start")
            assert _addr(f, f"cells/{cell}/chroms/name") == _addr(f, "chroms/name")
        assert "weight" not in f["bins"]


def test_scool_writer_rejects_bad_cells(tmp_path):
    from cooler.create._writer import _ScoolWriter

    path = str(tmp_path / "s.scool")
    with _ScoolWriter(path, BINS) as w:
        w.add_cell("a", pixels([(0, 0)]))
        with pytest.raises(ValueError, match="Invalid cell name"):
            w.add_cell("x/y", pixels([(0, 0)]))
        with pytest.raises(ValueError, match="already been added"):
            w.add_cell("a", pixels([(1, 1)]))
        other = cooler.util.binnify(pd.Series({"chr1": 60, "chr2": 40}), 20)
        with pytest.raises(ValueError, match="doesn't match"):
            w.add_cell("b", pixels([(0, 0)]), bins=other)
        # A failed cell leaves no trace, and earlier cells are intact
        with pytest.raises(BadInputError):
            w.add_cell("c", iter([pixels([(2, 2)]), pixels([(1, 1)])]))
        w.add_cell("d", pixels([(3, 3)]))

    assert cooler.fileops.list_scool_cells(path) == ["/cells/a", "/cells/d"]
    assert tmp_groups(path) == []
    with h5py.File(path, "r") as f:
        assert f.attrs["ncells"] == 2


def test_scool_writer_append_mode(tmp_path):
    from cooler.create._writer import _ScoolWriter

    path = str(tmp_path / "s.scool")
    with _ScoolWriter(path, BINS) as w:
        w.add_cell("a", pixels([(0, 0)]))
    with _ScoolWriter(path, BINS, mode="a") as w:
        w.add_cell("b", pixels([(1, 1)]))
        w.add_cell("a", pixels([(2, 2)]))  # replaces the existing cell
    assert cooler.fileops.list_scool_cells(path) == ["/cells/a", "/cells/b"]
    assert cooler.Cooler(path + "::/cells/a").pixels()[:]["bin1_id"].tolist() == [2]
    with h5py.File(path, "r") as f:
        assert f.attrs["ncells"] == 2

    other = cooler.util.binnify(pd.Series({"chr1": 60, "chr2": 40}), 20)
    with pytest.raises(ValueError, match="different bin table"):
        _ScoolWriter(path, other, mode="a")


def test_create_scool_unordered_chunks_and_keys(tmp_path):
    path = str(tmp_path / "s.scool")
    chunks = [pixels([(4, 5), (0, 1)]), pixels([(0, 1)])]
    cooler.create_scool(path, BINS, {"a": iter(chunks)}, ordered=False)
    px = cooler.Cooler(path + "::/cells/a").pixels()[:]
    assert px.values.tolist() == [[0, 1, 2], [4, 5, 1]]

    with pytest.raises(ValueError, match="matching keys"):
        cooler.create_scool(
            str(tmp_path / "k.scool"),
            {"a": BINS, "b": BINS},
            {"a": pixels([(0, 0)])},
        )


@pytest.fixture
def small_cool(tmp_path):
    path = str(tmp_path / "c.cool")
    cooler.create_cooler(path, BINS, pixels([(0, 0), (0, 5), (1, 1), (3, 9)]))
    return path


def test_column_writer_batches(small_cool):
    from cooler.create._writer import _ColumnWriter

    cols = {"a": np.float64, "b": np.int32}
    with _ColumnWriter(small_cool, "bins", cols, attrs={"a": {"note": "hi"}}) as w:
        w.write({"a": np.arange(4.0), "b": np.arange(4)})
        w.write(pd.DataFrame({"a": np.arange(4.0, 10.0), "b": np.arange(4, 10)}))
    bins = cooler.Cooler(small_cool).bins()[:]
    assert bins["a"].tolist() == list(np.arange(10.0))
    assert bins["b"].tolist() == list(range(10))
    with h5py.File(small_cool, "r") as f:
        assert f["bins/a"].attrs["note"] == "hi"
        assert f["bins/b"].dtype == np.int32
        assert ".__writing__columns" not in f


def test_column_writer_offsets(small_cool):
    from cooler.create._writer import _ColumnWriter

    with _ColumnWriter(small_cool, "pixels", {"x": None}) as w:
        w.write({"x": [2.0, 3.0]}, lo=2)
        w.write({"x": [0.0, 1.0]}, lo=0)
        with pytest.raises(ValueError, match="overlap"):
            w.write({"x": [9.0]}, lo=1)
        with pytest.raises(ValueError, match="out of range"):
            w.write({"x": [9.0]}, lo=4)
    assert cooler.Cooler(small_cool).pixels()[:]["x"].tolist() == [0.0, 1.0, 2.0, 3.0]


def test_column_writer_requires_every_row(small_cool):
    from cooler.create._writer import _ColumnWriter

    writer = _ColumnWriter(small_cool, "pixels", {"x": None})
    with pytest.raises(ValueError, match=r"\[2, 4\)"), writer as w:
        w.write({"x": [0.0, 1.0]})
    with h5py.File(small_cool, "r") as f:
        assert "x" not in f["pixels"]
        assert ".__writing__columns" not in f


def test_column_writer_replace(small_cool):
    from cooler.create._writer import _add_columns, _ColumnWriter

    _add_columns(small_cool, "bins", {"w": np.ones(10)})
    with pytest.raises(ValueError, match="already exist"):
        _add_columns(small_cool, "bins", {"w": np.zeros(10)})

    # A failed replacement keeps the existing column
    writer = _ColumnWriter(small_cool, "bins", {"w": None}, replace=True)
    with pytest.raises(RuntimeError), writer as w:
        w.write({"w": np.zeros(5)})
        raise RuntimeError("interrupted")
    assert cooler.Cooler(small_cool).bins()[:]["w"].tolist() == [1.0] * 10

    _add_columns(small_cool, "bins", {"w": np.zeros(10)}, replace=True)
    assert cooler.Cooler(small_cool).bins()[:]["w"].tolist() == [0.0] * 10


def test_column_writer_rejects_core_columns_and_tables(small_cool):
    from cooler.create._writer import _ColumnWriter

    for table, col in [("bins", "start"), ("pixels", "count"), ("chroms", "name")]:
        with pytest.raises(ValueError, match="core columns"):
            _ColumnWriter(small_cool, table, {col: None}, replace=True)
    with pytest.raises(ValueError, match="Unknown table"):
        _ColumnWriter(small_cool, "nope", {"x": None})


def test_column_writer_lock_and_empty_table(tmp_path):
    from cooler.create._writer import _ColumnWriter

    path = str(tmp_path / "e.cool")
    cooler.create_cooler(path, BINS, pixels([]))
    lock = threading.Lock()
    with _ColumnWriter(path, "pixels", {"x": None}, lock=lock) as w:
        assert w._file is None
        w.write({"x": []})
    assert "x" in cooler.Cooler(path).pixels().columns


def test_column_writer_polars_batches(small_cool):
    pl = pytest.importorskip("polars")
    from cooler.create._writer import _ColumnWriter

    lf = pl.LazyFrame({"count": [1, 1, 1, 1]})
    with _ColumnWriter(small_cool, "pixels", {"double": np.int64}) as w:
        lf.select(double=pl.col("count") * 2).sink_batches(w.write, chunk_size=3)
    assert cooler.Cooler(small_cool).pixels()[:]["double"].tolist() == [2, 2, 2, 2]


def test_balance_force_keeps_weights_on_failure(small_cool):
    from click.testing import CliRunner

    from cooler.cli.balance import balance

    runner = CliRunner()
    args = ["--ignore-diags", "0", "--mad-max", "0", "--min-nnz", "0", "-p", "1"]
    result = runner.invoke(balance, [small_cool, *args])
    assert result.exit_code == 0, result.output
    before = cooler.Cooler(small_cool).bins()[:]["weight"].to_numpy()

    # A forced rebalance that is discarded leaves the previous weights
    result = runner.invoke(
        balance,
        [
            small_cool,
            *args,
            "--force",
            "--max-iters",
            "1",
            "--convergence-policy",
            "discard",
        ],
    )
    assert result.exit_code == 0, result.output
    after = cooler.Cooler(small_cool).bins()[:]["weight"].to_numpy()
    np.testing.assert_array_equal(before, after)


def test_writer_guards(tmp_path):
    from cooler.create._writer import (
        _ColumnWriter,
        _ScoolWriter,
        _UnorderedCoolerWriter,
    )

    path = str(tmp_path / "g.cool")
    with pytest.raises(ValueError, match="DataFrame"):
        _CoolerWriter(path, BINS.to_dict())

    w = _CoolerWriter(path, BINS)
    with pytest.raises(ValueError, match="not found in input"):
        w.write({"bin1_id": [0], "bin2_id": [1]})
    w.finish()
    w.finish()  # finishing again does nothing
    w.abort()  # neither does aborting a finished writer
    with pytest.raises(ValueError, match="already finished"):
        w.write(pixels([(0, 1)]))

    with pytest.raises(ValueError, match="duplicates"):
        _UnorderedCoolerWriter(path, BINS, duplicates="sum")
    with pytest.raises(ValueError, match="unknown columns"):
        _UnorderedCoolerWriter(path, BINS, agg={"nope": "sum"})
    u = _UnorderedCoolerWriter(path, BINS)
    u.write(pixels([(0, 1)]))
    u.finish()
    u.finish()
    with pytest.raises(ValueError, match="already finished"):
        u.write(pixels([(0, 1)]))

    scool = str(tmp_path / "s.scool")
    with pytest.raises(ValueError, match="DataFrame"):
        _ScoolWriter(scool, BINS.to_dict())
    s = _ScoolWriter(scool, BINS, assembly="toy", metadata={"a": 1})
    s.add_cell("c", pixels([(0, 1)]))
    s.finish()
    s.finish()
    with pytest.raises(ValueError, match="already finished"):
        s.add_cell("d", pixels([(0, 1)]))
    cell = cooler.Cooler(scool + "::/cells/c")
    assert cell.info["genome-assembly"] == "toy"
    assert cell.info["metadata"] == {"a": 1}

    with pytest.raises(ValueError, match="No columns"):
        _ColumnWriter(path, "bins", {})
    with pytest.raises(ValueError, match="Unknown table"):
        _ColumnWriter(path, "nope", {"x": float})
    c = _ColumnWriter(path, "bins", {"x": float, "y": float})
    with pytest.raises(ValueError, match="missing from batch"):
        c.write({"x": np.zeros(10)})
    with pytest.raises(ValueError, match="same length"):
        c.write({"x": np.zeros(10), "y": np.zeros(9)})
    c.write({"x": np.zeros(10), "y": np.ones(10)})
    c.finish()
    c.finish()
    with pytest.raises(ValueError, match="already finished"):
        c.write({"x": np.zeros(10), "y": np.ones(10)})
    assert cooler.Cooler(path).bins()["y"][:].sum() == 10
