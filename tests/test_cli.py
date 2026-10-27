"""
Tests for the top-level CLI, the small commands (makebins, digest) and
the parsing helpers in cooler.cli._util.

To debug a failing command, run it with ``run_cli``: an exception raised by
the command is re-raised with its traceback, so ``pytest --pdb`` works.
"""
import os
import pdb  # noqa: T100 (patched below, not used)
import sys

import click
import numpy as np
import pytest
from _common import DATA

from cooler.cli import _util as util

TOY = DATA / "toy.symm.upper.2.cool"
CHROMSIZES = DATA / "toy.chrom.sizes"


def test_cli(run_cli):
    assert "Usage" in run_cli("-h").stdout
    run_cli("-V")
    run_cli("info", TOY)
    run_cli("-v", "info", TOY)
    run_cli("no-such-command", exit_code=2)


def test_cli_process_info(run_cli, monkeypatch):
    psutil = pytest.importorskip("psutil")
    result = run_cli("-vv", "info", TOY)
    assert "PSINFO:'pid'" in result.stderr

    def gone():
        raise psutil.NoSuchProcess(pid=0)

    monkeypatch.setattr(psutil, "Process", gone)
    result = run_cli("-vv", "info", TOY)
    assert "Process no longer exists" in result.stderr


def test_cli_process_info_without_psutil(run_cli, monkeypatch):
    monkeypatch.setitem(sys.modules, "psutil", None)  # import fails
    result = run_cli("-vv", "info", TOY)
    assert "Install psutil" in result.stderr


def test_cli_debug(run_cli, monkeypatch, capsys):
    # -d installs a hook that prints the traceback and starts pdb post-mortem
    monkeypatch.setattr(sys, "excepthook", sys.excepthook)
    run_cli("-d", "info", TOY)
    assert sys.excepthook is not sys.__excepthook__

    calls = []
    monkeypatch.setattr(pdb, "pm", lambda: calls.append("pm"))
    try:
        raise ValueError("boom")
    except ValueError:
        sys.excepthook(*sys.exc_info())
    assert calls == ["pm"]
    assert "ValueError: boom" in capsys.readouterr().err


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
    pytest.importorskip("pyfaidx")
    pytest.importorskip("Bio")
    fasta = DATA / "toy.fasta"
    result = run_cli("digest", CHROMSIZES, fasta, "HindIII")
    assert result.stdout == "chr1\t0\t32\nchr2\t0\t32\n"

    result = run_cli(
        "digest", "--header", "--rel-ids", 1, CHROMSIZES, fasta, "HindIII"
    )
    assert result.stdout == "chrom\tstart\tend\tid\nchr1\t0\t32\t1\nchr2\t0\t32\t1\n"

    result = run_cli("digest", CHROMSIZES, fasta, "NoSuchEnzyme", exit_code=1)



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
