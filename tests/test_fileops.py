import os.path as op
import shutil

import h5py
import numpy as np
import pytest
from _common import cooler_cmp, isolated_filesystem

from cooler import fileops

testdir = op.realpath(op.dirname(__file__))


def test_is_cooler():
    assert not fileops.is_cooler(op.join(testdir, "data", "toy.chrom.sizes"))
    assert fileops.is_cooler(op.join(testdir, "data", "toy.symm.upper.2.cool"))
    assert not fileops.is_cooler(op.join(testdir, "data", "toy.symm.upper.2.mcool"))
    assert fileops.is_cooler(
        op.join(testdir, "data", "toy.symm.upper.2.mcool") + "::resolutions/2"
    )


def test_is_multires_file():
    assert not fileops.is_multires_file(op.join(testdir, "data", "toy.chrom.sizes"))
    assert not fileops.is_multires_file(
        op.join(testdir, "data", "toy.symm.upper.2.cool")
    )
    assert fileops.is_multires_file(op.join(testdir, "data", "toy.symm.upper.2.mcool"))


def test_list_coolers():
    listing = fileops.list_coolers(op.join(testdir, "data", "toy.symm.upper.2.mcool"))
    paths = set(listing)
    for path in (
        "/resolutions/2",
        "/resolutions/4",
        "/resolutions/8",
        "/resolutions/16",
        "/resolutions/32",
    ):
        assert path in paths


def test_ls_attr_tree():
    src_file = op.join(testdir, "data", "toy.symm.upper.2.mcool")
    with h5py.File(src_file, "r") as f:
        attr_dict = fileops.read_attr_tree(f, level=None)
    assert attr_dict["@attrs"]["format"] == "HDF5::MCOOL"


def test_ls_data_tree():
    with isolated_filesystem():
        src_file = op.join(testdir, "data", "toy.symm.upper.2.mcool")
        listing = fileops.ls(src_file + "::" + "resolutions/2")
        for path in [
            "/resolutions/2",
            "/resolutions/2/chroms",
            "/resolutions/2/chroms/name",
            "/resolutions/2/chroms/length",
            "/resolutions/2/bins",
            "/resolutions/2/bins/chrom",
            "/resolutions/2/bins/start",
            "/resolutions/2/bins/end",
            "/resolutions/2/pixels",
            "/resolutions/2/pixels/bin1_id",
            "/resolutions/2/pixels/bin2_id",
            "/resolutions/2/pixels/count",
        ]:
            assert path in listing


def test_cp():
    with isolated_filesystem():
        src_file = op.join(testdir, "data", "toy.symm.upper.2.mcool")

        # file-to-file
        src_uri = src_file + "::resolutions/2"
        fileops.cp(src_uri, "test.2.cool")
        cooler_cmp(src_uri, "test.2.cool")
        fileops.cp(src_uri, "test.2.cool::/nested/")
        cooler_cmp(src_uri, "test.2.cool::nested")

        # within-file
        test_file = "test.src.mcool"
        shutil.copyfile(src_file, test_file)
        fileops.cp(test_file + "::resolutions/2", test_file + "::abc/d")
        cooler_cmp(test_file + "::resolutions/2", test_file + "::abc/d")
        with h5py.File(test_file, mode="r") as f:
            assert "resolutions/2" in f
            assert "abc/d" in f
            assert f["resolutions/2"].id != f["abc/d"].id


def test_mv():
    with isolated_filesystem():
        ref_file = "test.ref.mcool"
        src_file = "test.src.mcool"
        shutil.copyfile(op.join(testdir, "data", "toy.symm.upper.2.mcool"), ref_file)
        shutil.copyfile(op.join(testdir, "data", "toy.symm.upper.2.mcool"), src_file)
        fileops.mv(src_file + "::resolutions/2", src_file + "::abc/d")
        with h5py.File(src_file, mode="r") as f:
            assert "resolutions/2" not in f
            assert "abc/d" in f
        cooler_cmp(ref_file + "::resolutions/2", src_file + "::abc/d")


def test_ln():
    with isolated_filesystem():
        src_file = op.join(testdir, "data", "toy.symm.upper.2.mcool")

        # within-file hard link
        test_file = "test.hardlink.mcool"
        shutil.copyfile(src_file, test_file)
        fileops.ln(test_file + "::resolutions/2", test_file + "::abc/d")
        with h5py.File(test_file, mode="r") as f:
            assert "resolutions/2" in f
            assert "abc/d" in f
            assert f["resolutions/2"].id == f["abc/d"].id
        cooler_cmp(test_file + "::resolutions/2", test_file + "::abc/d")

        # within-file soft link
        test_file = "test.softlink.mcool"
        shutil.copyfile(src_file, test_file)
        fileops.ln(test_file + "::resolutions/2", test_file + "::abc/d", soft=True)
        with h5py.File(test_file, mode="r") as f:
            assert "resolutions/2" in f
            assert "abc/d" in f
            assert f["resolutions/2"].id == f["abc/d"].id
        cooler_cmp(test_file + "::resolutions/2", test_file + "::abc/d")

        # between-file external link
        test_file = "test.extlink.mcool"
        dst_file = "test.dst.cool"
        shutil.copyfile(src_file, test_file)
        fileops.ln(test_file + "::resolutions/2", dst_file + "::abc/d", soft=True)
        cooler_cmp(test_file + "::resolutions/2", dst_file + "::abc/d")


def test_print_trees():
    src_file = op.join(testdir, "data", "toy.symm.upper.2.mcool")
    fileops.pprint_attr_tree(src_file, level=3)
    fileops.pprint_data_tree(src_file, level=3)

    with h5py.File(src_file) as f:
        t = fileops.TreeViewer(f)
        t._repr_mimebundle_()


def test_is_scool_file():
    src_file = op.join(testdir, "data", "scool_test_file.scool")
    assert fileops.is_scool_file(src_file)


def test_list_scool_cells():
    src_file = op.join(testdir, "data", "scool_test_file.scool")
    paths = [
        "/cells/GSM2687248_41669_ACAGTG-R1-DpnII.100000.cool",
        "/cells/GSM2687249_41670_GGCTAC-R1-DpnII.100000.cool",
        "/cells/GSM2687250_41671_TTAGGC-R1-DpnII.100000.cool",
        "/cells/GSM2687251_41672_AGTTCC-R1-DpnII.100000.cool",
        "/cells/GSM2687252_41673_CCGTCC-R1-DpnII.100000.cool",
    ]
    cell_paths = fileops.list_scool_cells(src_file)
    assert len(cell_paths) == 5
    assert all(cell in paths for cell in cell_paths)


MCOOL = op.join(testdir, "data", "toy.symm.upper.2.mcool")
COOL = op.join(testdir, "data", "toy.symm.upper.2.cool")


def test_failed_copy_leaves_destination_intact(tmp_path):
    dst = str(tmp_path / "dst.mcool")
    shutil.copy(MCOOL, dst)
    before = fileops.list_coolers(dst)

    # Missing source group
    with pytest.raises(KeyError):
        fileops.cp(COOL + "::/typo", dst + "::/x", overwrite=True)
    # Hard link across files
    with pytest.raises(OSError):
        fileops.ln(COOL, dst + "::/x", overwrite=True)
    # Move across files
    with pytest.raises(OSError):
        fileops.mv(COOL, dst + "::/x", overwrite=True)

    assert fileops.list_coolers(dst) == before

    # A failed operation doesn't create the destination file either
    new = str(tmp_path / "new.cool")
    with pytest.raises(OSError):
        fileops.ln(COOL, new + "::/x")
    assert not op.exists(new)


def test_overwrite_replaces_only_destination_group(tmp_path):
    dst = str(tmp_path / "dst.mcool")
    shutil.copy(MCOOL, dst)
    before = fileops.list_coolers(dst)

    # Across files: only the target group is replaced
    with pytest.raises(FileExistsError):
        fileops.cp(COOL, dst + "::resolutions/2")
    fileops.cp(COOL, dst + "::resolutions/2", overwrite=True)
    assert fileops.list_coolers(dst) == before

    # Within a file
    with pytest.raises(FileExistsError):
        fileops.cp(dst + "::resolutions/4", dst + "::resolutions/2")
    fileops.cp(dst + "::resolutions/4", dst + "::resolutions/2", overwrite=True)
    assert fileops.list_coolers(dst) == before
    with h5py.File(dst, "r") as f:
        assert f["resolutions/2"].attrs["bin-size"] == 4

    # Replacing the root of another file replaces the whole file
    with pytest.raises(FileExistsError):
        fileops.cp(COOL, dst)
    fileops.cp(COOL, dst, overwrite=True)
    assert fileops.list_coolers(dst) == ["/"]


def test_copy_within_file_by_different_path(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    shutil.copy(MCOOL, "e.mcool")
    before = fileops.list_coolers("e.mcool")
    fileops.cp("./e.mcool::resolutions/2", "e.mcool::/x")
    fileops.ln("./e.mcool::/x", "e.mcool::/y")
    assert sorted(fileops.list_coolers("e.mcool")) == sorted([*before, "/x", "/y"])


def test_copy_into_itself_is_rejected(tmp_path):
    path = str(tmp_path / "f.mcool")
    shutil.copy(MCOOL, path)
    with pytest.raises(ValueError):
        fileops.cp(path + "::resolutions", path + "::resolutions/nested")
    with pytest.raises(ValueError):
        fileops.ln(path + "::/", path + "::/loop")
    with pytest.raises(ValueError):
        fileops.mv(path + "::resolutions/2", path + "::resolutions", overwrite=True)
    assert len(fileops.list_coolers(path)) == 5


def test_external_link_relative_to_linking_file(tmp_path, monkeypatch):
    import cooler

    sub = tmp_path / "sub"
    sub.mkdir()
    shutil.copy(COOL, sub / "src.cool")

    # Create the link from the parent directory, then read it from elsewhere
    monkeypatch.chdir(tmp_path)
    fileops.ln("sub/src.cool", "sub/dst.cool::/lnk", soft=True)
    monkeypatch.chdir(sub)
    assert fileops.list_coolers("dst.cool") == ["/lnk"]
    assert cooler.Cooler("dst.cool::/lnk").info["nbins"] > 0


def test_traversal_with_links(tmp_path):
    path = str(tmp_path / "f.mcool")
    shutil.copy(MCOOL, path)
    before = fileops.list_coolers(path)

    # Dangling links are skipped and links back to an ancestor aren't followed
    with h5py.File(path, "r+") as f:
        f["dangling"] = h5py.SoftLink("/does/not/exist")
        f["resolutions/2/loop"] = h5py.SoftLink("/resolutions")
    assert fileops.list_coolers(path) == before
    fileops.pprint_data_tree(path, None)
    fileops.pprint_attr_tree(path, None)

    # A hard-linked group is listed under every path it can be reached by
    with h5py.File(path, "r+") as f:
        del f["resolutions/2/loop"]
        f["hard"] = f["resolutions/2"]
    assert sorted(fileops.list_coolers(path)) == sorted([*before, "/hard"])


def test_attr_tree_with_array_attributes(tmp_path):
    path = str(tmp_path / "a.cool")
    shutil.copy(COOL, path)
    with h5py.File(path, "r+") as f:
        f.attrs["arr"] = [1, 2, 3]
        f.attrs["one"] = np.array([7])
    with h5py.File(path, "r") as f:
        attrs = fileops.read_attr_tree(f, None)["@attrs"]
    assert attrs["arr"] == [1, 2, 3]
    assert attrs["one"] == [7]


def test_detection_returns_false(tmp_path):
    assert not fileops.is_cooler(MCOOL + "::/nope")
    not_hdf5 = tmp_path / "x.txt"
    not_hdf5.write_text("hello")
    assert not fileops.is_scool_file(str(not_hdf5))
    assert not fileops.is_scool_file(str(tmp_path / "missing.scool"))
