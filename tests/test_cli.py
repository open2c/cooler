"""
Tests for the top-level CLI, the small commands (makebins, digest, csort) and
the parsing helpers in cooler.cli._util.

To debug a failing command, run it with ``run_cli``: an exception raised by
the command is re-raised with its traceback, so ``pytest --pdb`` works.
"""
import os

import click
import numpy as np
import pytest
from _common import DATA

from cooler.cli import _util as util
from cooler.util import cmd_exists

TOY = DATA / "toy.symm.upper.2.cool"
CHROMSIZES = DATA / "toy.chrom.sizes"


def test_cli(run_cli):
    assert "Usage" in run_cli("-h").stdout
    run_cli("-V")
    run_cli("info", TOY)
    run_cli("-v", "info", TOY)
    run_cli("no-such-command", exit_code=2)


def test_cli_process_info(run_cli):
    pytest.importorskip("psutil")
    result = run_cli("-vv", "info", TOY)
    assert "PSINFO:'pid'" in result.stderr


def test_makebins(run_cli):
    result = run_cli("makebins", CHROMSIZES, 8)
    assert result.stdout == (
        "chr1\t0\t8\nchr1\t8\t16\nchr1\t16\t24\nchr1\t24\t32\n"
        "chr2\t0\t8\nchr2\t8\t16\nchr2\t16\t24\nchr2\t24\t32\n"
    )

    result = run_cli("makebins", "--header", "--rel-ids", 1, CHROMSIZES, 8)
    assert result.stdout == (
        "chrom\tstart\tend\tid\n"
        "chr1\t0\t8\t1\nchr1\t8\t16\t2\nchr1\t16\t24\t3\nchr1\t24\t32\t4\n"
        "chr2\t0\t8\t1\nchr2\t8\t16\t2\nchr2\t16\t24\t3\nchr2\t24\t32\t4\n"
    )


@pytest.mark.skipif(
    os.name == "nt", reason="Non-reproducible off-by-one bug on Win32 GitHub runner"
)
def test_digest(run_cli):
    fasta = DATA / "toy.fasta"
    result = run_cli("digest", CHROMSIZES, fasta, "HindIII")
    assert result.stdout == "chr1\t0\t32\nchr2\t0\t32\n"

    result = run_cli(
        "digest", "--header", "--rel-ids", 1, CHROMSIZES, fasta, "HindIII"
    )
    assert result.stdout == "chrom\tstart\tend\tid\nchr1\t0\t32\t1\nchr2\t0\t32\t1\n"

    result = run_cli("digest", CHROMSIZES, fasta, "NoSuchEnzyme", exit_code=1)


@pytest.mark.skipif(
    not cmd_exists("bgzip") or not cmd_exists("pairix"),
    reason="missing command line requirements for csort",
)
def test_csort(run_cli, tmp_path):
    out = tmp_path / "test.pairs.blksrt"
    run_cli(
        "csort", "-c1", 2, "-p1", 3, "-c2", 4, "-p2", 5, "-p", 1,
        DATA / "toy.pairs", CHROMSIZES, "-o", out,
    )
    assert out.exists()


def test_delimited_tuple():
    param = util.DelimitedTuple(sep=",")
    assert param.convert("spam,eggs", None, None) == ("spam", "eggs")
    assert param.convert(("spam", "eggs"), None, None) == ("spam", "eggs")
    assert param.convert(None, None, None) is None
    assert param.name == "separated[,]"

    assert util.DelimitedTuple(sep=";").convert("a;b", None, None) == ("a", "b")

    param = util.DelimitedTuple(sep=",", type=int)
    assert param.convert("10,20", None, None) == (10, 20)
    with pytest.raises(click.BadParameter):
        param.convert("spam,20", None, None)


def test_parse_kv_list_param():
    result = util.parse_kv_list_param("spam=1,eggs=2,apples=oranges,x=[0,0,0]")
    assert result == {"spam": 1, "eggs": 2, "apples": "oranges", "x": [0, 0, 0]}
    assert util.parse_kv_list_param("a=1;b=2", item_sep=";") == {"a": 1, "b": 2}
    with pytest.raises(click.BadParameter):
        util.parse_kv_list_param("a=[1,")


def test_parse_field_param():
    assert util.parse_field_param("count=7:dtype=int,agg=sum") == (
        "count", 6, np.dtype(int), "sum"
    )
    assert util.parse_field_param("count:agg=max", includes_colnum=False) == (
        "count", None, None, "max"
    )
    assert util.parse_field_param("count") == ("count", None, None, None)

    for arg, kwargs in [
        ("a:b:c", {}),  # more than one ':'
        ("a=x", {}),  # not a number
        ("a=0", {}),  # field numbers start at 1
        ("a=1=2", {}),
        ("a=1:dtype", {}),  # property without a value
        ("a=1:color=red", {}),  # unknown property
        ("a=1:dtype=flaot", {}),  # unknown dtype
        ("a=1:agg=sum", {"includes_agg": False}),
        ("a=1:agg=mean", {"aggs": ("sum", "max")}),
    ]:
        with pytest.raises(click.BadParameter):
            util.parse_field_param(arg, **kwargs)


def test_parse_bins(tmp_path):
    cs1, bins1 = util.parse_bins(f"{CHROMSIZES}:8")
    cs2, bins2 = util.parse_bins(str(DATA / "toy.bins.8.bed"))
    assert cs1.equals(cs2)
    assert (bins1 == bins2).all().all()

    with pytest.raises(ValueError, match="not found"):
        util.parse_bins(f"{tmp_path / 'missing.sizes'}:8")
    with pytest.raises(ValueError, match="integer binsize"):
        util.parse_bins(f"{CHROMSIZES}:8kb")
    with pytest.raises(ValueError, match="Expected BINS"):
        util.parse_bins(str(tmp_path / "missing.bed"))
