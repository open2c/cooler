"""
Tests for the commands that work on files and groups: cp, mv, ln, ls, tree and
attrs.
"""
import yaml
from _common import DATA

import cooler

TOY = DATA / "toy.symm.upper.2.cool"


def test_cp_mv_ln(run_cli, tmp_path):
    test = tmp_path / "test.cool"
    test2 = tmp_path / "test2.cool"
    run_cli("cp", TOY, test)
    run_cli("cp", test, f"{test2}::some/old_path")
    run_cli("mv", f"{test2}::some/old_path", f"{test2}::some/path")
    assert cooler.fileops.list_coolers(test2) == ["/some/path"]

    # mv only renames within a file
    run_cli("mv", test, f"{tmp_path / 'test4.cool'}::x", exit_code=1)
    assert not (tmp_path / "test4.cool").exists()

    run_cli("ln", f"{test2}::some/path", f"{test2}::hard/link")
    run_cli("ln", "-s", f"{test2}::some/path", f"{test2}::soft/link")
    run_cli("ln", "-s", f"{test2}::some/path", f"{tmp_path / 'test3.cool'}::ext/link")
    assert cooler.fileops.list_coolers(test2) == [
        "/hard/link", "/soft/link", "/some/path"
    ]
    assert cooler.fileops.list_coolers(tmp_path / "test3.cool") == ["/ext/link"]


def test_ls(run_cli):
    assert run_cli("ls", TOY).stdout == f"{TOY}::/\n"
    assert run_cli("ls", "-l", TOY).stdout == f"{TOY}::/\t2\n"
    var = DATA / "toy.symm.upper.var.cool"
    assert run_cli("ls", "-l", var).stdout == f"{var}::/\t<variable>\n"


def test_tree(run_cli):
    output = run_cli("tree", TOY).stdout
    for name in ["chroms", "bins", "pixels", "indexes"]:
        assert name in output


def test_attrs(run_cli):
    attrs = yaml.safe_load(run_cli("attrs", TOY).stdout)
    assert attrs["@attrs"]["format"] == "HDF5::Cooler"
