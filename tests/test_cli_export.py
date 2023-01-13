import os.path as op
from io import StringIO

import numpy as np
import pandas as pd
from _common import cooler_cmp
from click.testing import CliRunner

import cooler
from cooler.cli.dump import dump

### EXPORT ###
from cooler.cli.info import info
from cooler.cli.show import show

# import pytest


testdir = op.realpath(op.dirname(__file__))
datadir = op.join(testdir, "data")


def test_info():
    runner = CliRunner()
    with runner.isolated_filesystem():
        f_in = op.join(datadir, "toy.symm.upper.2.cool")
        result = runner.invoke(info, [f_in])
        assert result.exit_code == 0

        result = runner.invoke(info, [f_in, "--field", "bin-type"])
        assert result.exit_code == 0

        result = runner.invoke(info, [f_in, "--field", "doesnotexist"])
        assert result.exit_code > 0

        result = runner.invoke(info, [f_in, "--metadata"])
        assert result.exit_code == 0


def test_dump():
    runner = CliRunner()
    with runner.isolated_filesystem():
        f_in = op.join(datadir, "toy.symm.upper.2.cool")
        result = runner.invoke(dump, [f_in])
        assert result.exit_code == 0
        result = runner.invoke(dump, [f_in, "-t", "chroms", "--columns", "length"])
        assert result.exit_code == 0
        result = runner.invoke(dump, [f_in, "-t", "bins", "--columns", "chrom,start"])
        assert result.exit_code == 0
        result = runner.invoke(dump, [f_in, "-r", "chr1"])
        assert result.exit_code == 0
        result = runner.invoke(dump, [f_in, "-r", "chr1:0-16", "-r2", "chr1:10-25"])
        assert result.exit_code == 0
        result = runner.invoke(dump, [f_in, "-r", "chr1:10-25", "-r2", "chr1:0-5"])
        assert result.exit_code == 0
        result = runner.invoke(dump, [f_in, "--join"])
        assert result.exit_code == 0
        result = runner.invoke(dump, [f_in, "--join", "--one-based-ids"])
        assert result.exit_code == 0
        result = runner.invoke(dump, [f_in, "--join", "--one-based-starts"])
        assert result.exit_code == 0
        result = runner.invoke(dump, [f_in, "--annotate", "chrom", "--one-based-starts"])
        assert result.exit_code == 0

        # unbalanced file
        result = runner.invoke(dump, [f_in, "-b"])
        assert result.exit_code == 1

        # roundtrip symm-upper data
        result = runner.invoke(dump, [f_in, "-H", "-t", "bins"])
        bins = pd.read_csv(StringIO(result.output), sep="\t")
        result = runner.invoke(dump, [f_in, "-H"])
        pixels = pd.read_csv(StringIO(result.output), sep="\t")
        cooler.create_cooler("out.cool", bins, pixels, symmetric_upper=True)
        cooler_cmp(f_in, "out.cool")

        # duplexed output
        result = runner.invoke(dump, [f_in, "--fill-lower", "-H"])
        pixels2 = pd.read_csv(StringIO(result.output), sep="\t")
        assert len(pixels2) > len(pixels)
        upper = pixels2[pixels2["bin1_id"] <= pixels2["bin2_id"]].reset_index(drop=True)
        assert np.allclose(pixels, upper)

        # lower triangle
        result = runner.invoke(dump, [f_in, "-H", "-r", "chr2", "-r2", "chr1"])
        trans_lower = pd.read_csv(StringIO(result.output), sep="\t")
        assert len(trans_lower) == 0
        result = runner.invoke(dump, [f_in, "-m", "-H", "-r", "chr2", "-r2", "chr1"])
        trans_lower = pd.read_csv(StringIO(result.output), sep="\t")
        assert len(trans_lower) > 0

        # roundtrip square data
        f_in = op.join(datadir, "toy.asymm.2.cool")
        result = runner.invoke(dump, [f_in, "-H", "-t", "bins"])
        bins = pd.read_csv(StringIO(result.output), sep="\t")
        result = runner.invoke(dump, [f_in, "-H"])
        pixels = pd.read_csv(StringIO(result.output), sep="\t")
        cooler.create_cooler("out.cool", bins, pixels, symmetric_upper=False)
        cooler_cmp(f_in, "out.cool")
        result = runner.invoke(dump, [f_in, "--fill-lower", "-H"])
        pixels2 = pd.read_csv(StringIO(result.output), sep="\t")
        assert np.allclose(pixels, pixels2)

        # for square data, -f (--fill-lower) is a no-op
        result = runner.invoke(dump, [f_in, "-H", "-r", "chr2", "-r2", "chr1"])
        lower1 = pd.read_csv(StringIO(result.output), sep="\t")
        result = runner.invoke(dump, [f_in, "-f", "-H", "-r", "chr2", "-r2", "chr1"])
        lower2 = pd.read_csv(StringIO(result.output), sep="\t")
        assert np.allclose(lower1, lower2)


def test_show():
    runner = CliRunner()
    with runner.isolated_filesystem():
        f_in = op.join(datadir, "toy.symm.upper.2.cool")
        result = runner.invoke(show, [f_in, 'chr1', '-o', 'bla'])
        assert result.exit_code == 0
