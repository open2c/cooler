import os.path as op
import shutil

import h5py
from _common import cooler_cmp, isolated_filesystem

from cooler import fileops

testdir = op.realpath(op.dirname(__file__))


def test_is_cooler():
    assert not fileops.is_cooler(op.join(testdir, "data", "toy.chrom.sizes"))
    assert fileops.is_cooler(op.join(testdir, "data", "toy.symm.upper.2.cool"))
    assert not fileops.is_cooler(op.join(testdir, "data", "toy.symm.upper.2.mcool"))
    assert fileops.is_cooler(
        op.join(testdir, "data", "toy.symm.upper.2.mcool") + '::resolutions/2'
    )


def test_is_multires_file():
    assert not fileops.is_multires_file(
        op.join(testdir, "data", "toy.chrom.sizes")
    )
    assert not fileops.is_multires_file(
        op.join(testdir, "data", "toy.symm.upper.2.cool")
    )
    assert fileops.is_multires_file(
        op.join(testdir, "data", "toy.symm.upper.2.mcool")
    )


def test_list_coolers():
    listing = fileops.list_coolers(
        op.join(testdir, "data", "toy.symm.upper.2.mcool")
    )
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
        listing = fileops.ls(src_file + '::' + 'resolutions/2')
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
        with h5py.File(test_file, mode='r') as f:
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
        with h5py.File(src_file, mode='r') as f:
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
        with h5py.File(test_file, mode='r') as f:
            assert "resolutions/2" in f
            assert "abc/d" in f
            assert f["resolutions/2"].id == f["abc/d"].id
        cooler_cmp(test_file + "::resolutions/2", test_file + "::abc/d")

        # within-file soft link
        test_file = "test.softlink.mcool"
        shutil.copyfile(src_file, test_file)
        fileops.ln(test_file + "::resolutions/2", test_file + "::abc/d", soft=True)
        with h5py.File(test_file, mode='r') as f:
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
    src_file = op.join(testdir, "data", 'scool_test_file.scool')
    assert fileops.is_scool_file(src_file)


def test_list_scool_cells():
    src_file = op.join(testdir, "data", 'scool_test_file.scool')
    paths = ['/cells/GSM2687248_41669_ACAGTG-R1-DpnII.100000.cool', '/cells/GSM2687249_41670_GGCTAC-R1-DpnII.100000.cool',
             '/cells/GSM2687250_41671_TTAGGC-R1-DpnII.100000.cool', '/cells/GSM2687251_41672_AGTTCC-R1-DpnII.100000.cool',
             '/cells/GSM2687252_41673_CCGTCC-R1-DpnII.100000.cool']
    cell_paths = fileops.list_scool_cells(src_file)
    assert len(cell_paths) == 5
    for cell in paths:
        if cell not in cell_paths:
            assert False
