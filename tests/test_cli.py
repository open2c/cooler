"""

Useful things for debugging CLI failures:
pytest --pdb
print(result.output)
import traceback; traceback.print_tb(result.exception.__traceback__)

"""
import os.path as op

import click
import pytest
from click.testing import CliRunner

from cooler.cli import _util as util
from cooler.cli import cli
from cooler.util import cmd_exists

testdir = op.realpath(op.dirname(__file__))
datadir = op.join(testdir, "data")


def test_cli():
    runner = CliRunner()

    result = runner.invoke(cli, ["-h"])
    assert result.exit_code == 0

    result = runner.invoke(cli, ["-V"])
    assert result.exit_code == 0


def test_verbose():
    runner = CliRunner()

    result = runner.invoke(cli, ["info", op.join(datadir, "toy.symm.upper.2.cool")])
    assert result.exit_code == 0

    result = runner.invoke(cli, ["-v", "info", op.join(datadir, "toy.symm.upper.2.cool")])
    assert result.exit_code == 0

    pytest.importorskip("psutil")
    result = runner.invoke(cli, ["-vv", "info", op.join(datadir, "toy.symm.upper.2.cool")])
    assert result.exit_code == 0


### HELPERS ###
def test_makebins():
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "makebins",
            op.join(datadir, "toy.chrom.sizes"),
            '8'
        ]
    )
    assert result.exit_code == 0
    assert result.output == "chr1\t0\t8\nchr1\t8\t16\nchr1\t16\t24\nchr1\t24\t32\nchr2\t0\t8\nchr2\t8\t16\nchr2\t16\t24\nchr2\t24\t32\n"

    result = runner.invoke(
        cli,
        [
            "makebins",
            "--header",
            "--rel-ids", "1",
            op.join(datadir, "toy.chrom.sizes"),
            '8'
        ]
    )
    assert result.exit_code == 0
    assert result.output == "chrom\tstart\tend\tid\nchr1\t0\t8\t1\nchr1\t8\t16\t2\nchr1\t16\t24\t3\nchr1\t24\t32\t4\nchr2\t0\t8\t1\nchr2\t8\t16\t2\nchr2\t16\t24\t3\nchr2\t24\t32\t4\n"


def test_digest():
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "digest",
            op.join(datadir, "toy.chrom.sizes"),
            op.join(datadir, "toy.fasta"),
            'HindIII'
        ]
    )
    assert result.exit_code == 0
    assert result.output == "chr1\t0\t32\nchr2\t0\t32\n"

    result = runner.invoke(
        cli,
        [
            "digest",
            "--header",
            "--rel-ids", "1",
            op.join(datadir, "toy.chrom.sizes"),
            op.join(datadir, "toy.fasta"),
            'HindIII'
        ]
    )
    assert result.exit_code == 0
    assert result.output == "chrom\tstart\tend\tid\nchr1\t0\t32\t1\nchr2\t0\t32\t1\n"


@pytest.mark.skipif(
    not cmd_exists('bgzip') or not cmd_exists('pairix'),
    reason="missing command line requirements for csort"
)
def test_csort():
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(
            cli,
            [
                "csort",
                '-c1', '2',
                '-p1', '3',
                '-c2', '4',
                '-p2', '5',
                '-p', 1,
                op.join(datadir, "toy.pairs"),
                op.join(datadir, "toy.chrom.sizes"),
                '-o', 'test.pairs.blksrt'
            ]
        )
        assert result.exit_code == 0


def test_util_delimited_tuple():
    param = util.DelimitedTuple(sep=',')
    result = param.convert('spam,eggs', None, None)
    assert result == ('spam', 'eggs')

    param = util.DelimitedTuple(sep=',', type=int)
    result = param.convert('10,20', None, None)
    assert result == (10, 20)
    with pytest.raises(click.BadParameter):
        param.convert('spam,20', None, None)


def test_util_parse_kv_list_param():
    result = util.parse_kv_list_param(
        'spam=1,eggs=2,apples=oranges,x=[0,0,0]'
    )
    assert result == {
        'spam': 1,
        'eggs': 2,
        'apples': 'oranges',
        'x': [0, 0, 0],
    }


def test_util_parse_field_param():
    util.parse_field_param('count=7:dtype=int,agg=sum')
    util.parse_field_param('count:dtype=int,agg=sum', includes_colnum=False)
    util.parse_field_param('count:dtype=int', includes_agg=False)


def test_parse_bins():
    cs1, bins1 = util.parse_bins(op.join(datadir, 'toy.chrom.sizes') + ':8')
    cs2, bins2 = util.parse_bins(op.join(datadir, 'toy.bins.8.bed'))
    assert cs1.equals(cs2)
    assert (bins1 == bins2).all().all()


def test_check_ncpus():
    assert util.check_ncpus(10) <= 10
