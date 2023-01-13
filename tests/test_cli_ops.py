import os.path as op

import numpy as np
import pandas as pd
from click.testing import CliRunner

import cooler

### COMPUTE ###
from cooler.cli.balance import balance
from cooler.cli.coarsen import coarsen

## REDUCTION ###
from cooler.cli.merge import merge
from cooler.cli.zoomify import zoomify

# import pytest


testdir = op.realpath(op.dirname(__file__))
datadir = op.join(testdir, "data")


def test_merge():
    runner = CliRunner()
    with runner.isolated_filesystem():
        f_in = op.join(datadir, "toy.symm.upper.2.cool")
        result = runner.invoke(
            merge, ["toy.2.double.cool", f_in, f_in, "--field", "count:dtype=int"]
        )
        assert result.exit_code == 0
        total1 = cooler.Cooler(f_in).pixels()["count"][:].sum()
        total2 = cooler.Cooler("toy.2.double.cool").pixels()["count"][:].sum()
        assert total2 == 2 * total1


def test_coarsen():
    runner = CliRunner()
    with runner.isolated_filesystem():
        f_in = op.join(datadir, "toy.symm.upper.2.cool")
        f_ref = op.join(datadir, "toy.symm.upper.4.cool")
        result = runner.invoke(
            coarsen, [
                f_in,
                "--factor", "2",
                "--nproc", "2",
                "-o", "toy.2.coarsen_2.cool"
            ],
        )
        assert result.exit_code == 0
        pix1 = cooler.Cooler(f_ref).pixels()["count"][:]
        pix2 = cooler.Cooler("toy.2.coarsen_2.cool").pixels()["count"][:]
        assert np.allclose(pix1, pix2)

        result = runner.invoke(
            coarsen, [
                f_in,
                "--factor", "2",
                "--field", "count:dtype=float,agg=mean",
                "-o", "toy.2.coarsen_2_mean.cool"
            ],
        )
        assert result.exit_code == 0
        pix2 = cooler.Cooler("toy.2.coarsen_2_mean.cool").pixels()["count"][:]
        assert pix2.dtype.kind == 'f'


def test_zoomify():
    runner = CliRunner()

    with runner.isolated_filesystem():
        f_in = op.join(datadir, "toy.symm.upper.2.cool")
        result = runner.invoke(
            zoomify, [f_in, "--balance", "--legacy", "-o", "toy.2.mcool"]
        )
        assert result.exit_code == 0

        f_in = op.join(datadir, "toy.symm.upper.2.cool")
        result = runner.invoke(
            zoomify, [
                f_in,
                "--balance",
                "--nproc", "2",
                "-o", "toy.2.mcool"
            ]
        )
        assert result.exit_code == 0

        f_in = op.join(datadir, "toy.symm.upper.2.cool")
        result = runner.invoke(
            zoomify, [
                f_in,
                "--balance",
                "--resolutions", "2,4,8",
                "-o", "toy.2.mcool"
            ]
        )
        assert result.exit_code == 0

        f_in = op.join(datadir, "toy.symm.upper.2.cool")
        result = runner.invoke(
            zoomify, [
                f_in,
                "--balance",
                "--resolutions", "2,4,8",
                "--field", "count:dtype=float,agg=mean",
                "-o", "toy.2.mcool"
            ]
        )
        assert result.exit_code == 0
        # pix1 = cooler.Cooler(f_ref).pixels()['count'][:]
        # pix2 = cooler.Cooler('toy.4.cool').pixels()['count'][:]
        # assert np.allclose(pix1, pix2)


def test_balance():
    runner = CliRunner()
    with runner.isolated_filesystem():
        f_in = op.join(datadir, "toy.symm.upper.2.cool")
        result = runner.invoke(
            balance, [
                f_in,
                "--ignore-diags", "2",
                "--mad-max", "0",
                "--min-nnz", "0",
                "--tol", "0.05",
                "--nproc", "2",
                "--stdout",
            ],
        )
        assert result.exit_code == 0
        assert len(result.output.split("\n")) == 32

        # convergence
        result = runner.invoke(
            balance, [
                f_in,
                "--ignore-diags", "2",
                "--mad-max", "0",
                "--min-nnz", "0",
                "--tol", "0.05",
                "--max-iters", "1",
                "--convergence-policy", "store_final",
                "--stdout",
            ],
        )
        assert result.exit_code == 0
        assert len(result.output)
        result = runner.invoke(
            balance, [
                f_in,
                "--ignore-diags", "2",
                "--mad-max", "0",
                "--min-nnz", "0",
                "--tol", "0.05",
                "--max-iters", "1",
                "--convergence-policy", "discard",
                "--stdout",
            ],
        )
        assert result.exit_code == 0
        assert not result.output
        result = runner.invoke(
            balance, [
                f_in,
                "--ignore-diags", "2",
                "--mad-max", "0",
                "--min-nnz", "0",
                "--tol", "0.05",
                "--max-iters", "1",
                "--convergence-policy", "error",
                "--stdout",
            ],
        )
        assert result.exit_code == 1

        # file is unbalanced
        result = runner.invoke(
            balance, [
                f_in,
                "--check"
            ],
        )
        assert result.exit_code == 1

        # blacklisting regions
        blacklist = pd.DataFrame({
            'chrom': ['chr1', 'chr2'],
            'start': [5, 10],
            'end': [10, 20],
        }, columns=['chrom', 'start', 'end'])
        blacklist.to_csv(
            'blacklist.bed', sep='\t', index=False, header=False
        )
        result = runner.invoke(
            balance, [
                f_in,
                "--ignore-diags", "2",
                "--mad-max", "0",
                "--min-nnz", "0",
                "--tol", "0.05",
                "--blacklist", "blacklist.bed",
                "--stdout",
            ],
        )
        assert result.exit_code == 0
