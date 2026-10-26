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
        [0, 0], [0, 5], [1, 1], [3, 9]
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
