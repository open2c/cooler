import os
import shutil
import tempfile

import h5py
import numpy as np
import pytest
from click.testing import CliRunner

import cooler
from cooler.cli.balance import balance

testdir = os.path.dirname(os.path.realpath(__file__))


@pytest.fixture
def temp_cool(tmpdir):
    def _copy_cool(fp):
        temp_fp = os.path.join(tmpdir, os.path.basename(fp))
        shutil.copy(fp, temp_fp)
        return temp_fp
    return _copy_cool


@pytest.mark.parametrize(
    "fp,tol",
    [(os.path.join(testdir, "data", "hg19.GM12878-MboI.matrix.2000kb.cool"), 1e-2)],
)
def test_balancing_genomewide(temp_cool, fp, tol):
    runner = CliRunner()
    temp_fp = temp_cool(fp)

    args = [
        temp_fp, "--ignore-diags", "1", "--min-nnz", "10",
        "--tol", str(tol), "--force", "--nproc", "1"
    ]
    result = runner.invoke(balance, args)
    assert result.exit_code == 0, f"Command failed: {result.output}"

    with h5py.File(temp_fp, "r") as h5:
        weights = h5["/bins/weight"][:]
        assert not np.all(np.isnan(weights)), "Weights are all NaN"

        clr = cooler.Cooler(temp_fp)
        mat = clr.matrix(balance=False, sparse=True)[:, :]
        mat.data = weights[mat.row] * weights[mat.col] * mat.data
        arr = mat.toarray()

        mask = np.isnan(weights)
        arr[mask, :] = 0
        arr[:, mask] = 0
        np.fill_diagonal(arr, 0)

        marg = np.sum(arr, axis=0)
        var = np.var(marg[marg != 0])
        assert var < tol, f"Variance {var} exceeds tolerance {tol}"

        conv_marg = marg[~np.isnan(marg)].mean()
        err_marg = marg[~np.isnan(marg)].std()
        assert np.isclose(conv_marg, 1, atol=err_marg), "Marginal not unity"


@pytest.mark.filterwarnings("ignore")
@pytest.mark.parametrize(
    "fp,tol",
    [(os.path.join(testdir, "data", "hg19.GM12878-MboI.matrix.2000kb.cool"), 1e-2)],
)
def test_balancing_cisonly(temp_cool, fp, tol):
    runner = CliRunner()
    temp_fp = temp_cool(fp)

    args = [
        temp_fp, "--cis-only", "--ignore-diags", "1",
        "--min-nnz", "10", "--tol", str(tol), "--force", "--nproc", "1"
    ]
    result = runner.invoke(balance, args)
    assert result.exit_code == 0, f"Command failed: {result.output}"

    with h5py.File(temp_fp, "r") as h5:
        weights = h5["/bins/weight"][:]
        chrom_offsets = h5["/indexes/chrom_offset"][:]

        clr = cooler.Cooler(temp_fp)
        mat = clr.matrix(balance=False, sparse=True)[:, :]
        mat.data = weights[mat.row] * weights[mat.col] * mat.data
        arr = mat.toarray()

        mask = np.isnan(weights)
        arr[mask, :] = 0
        arr[:, mask] = 0
        np.fill_diagonal(arr, 0)

        from scipy.linalg import block_diag
        spans = list(zip(chrom_offsets[:-1], chrom_offsets[1:]))
        blocks = [np.ones((hi - lo,) * 2) for lo, hi in spans]
        mask = block_diag(*blocks).astype(bool)
        arr[~mask] = 0

        marg = np.sum(arr, axis=0)
        for lo, hi in spans:
            m = marg[lo:hi]
            m = m[m != 0]
            if len(m):
                var = np.var(m)
                assert var < tol, f"Cis variance {var} exceeds tolerance {tol}"
                conv_marg = m[~np.isnan(m)].mean()
                err_marg = m[~np.isnan(m)].std()
                assert np.isclose(conv_marg, 1,
                                  atol=err_marg), "Cis marginal not unity"


@pytest.mark.filterwarnings("ignore")
@pytest.mark.parametrize(
    "fp,tol",
    [(os.path.join(testdir, "data", "hg19.GM12878-MboI.matrix.2000kb.cool"), 1e-2)],
)
def test_balancing_transonly(temp_cool, fp, tol):
    runner = CliRunner()
    temp_fp = temp_cool(fp)

    args = [
        temp_fp, "--trans-only", "--ignore-diags", "1",
        "--min-nnz", "10", "--tol", str(tol), "--force", "--nproc", "1"
    ]
    result = runner.invoke(balance, args)
    assert result.exit_code == 0, f"Command failed: {result.output}"

    with h5py.File(temp_fp, "r") as h5:
        weights = h5["/bins/weight"][:]
        chrom_offsets = h5["/indexes/chrom_offset"][:]

        clr = cooler.Cooler(temp_fp)
        mat = clr.matrix(balance=False, sparse=True)[:, :]
        mat.data = weights[mat.row] * weights[mat.col] * mat.data
        arr = mat.toarray()

        mask = np.isnan(weights)
        arr[mask, :] = 0
        arr[:, mask] = 0

        for lo, hi in zip(chrom_offsets[:-1], chrom_offsets[1:]):
            arr[lo:hi, lo:hi] = np.nan

        marg = np.nansum(arr, axis=0)
        var = np.nanvar(marg[marg != 0])
        assert var < tol, f"Trans variance {var} exceeds tolerance {tol}"
        conv_marg = marg[~np.isnan(marg)].mean()
        err_marg = marg[~np.isnan(marg)].std()
        assert np.isclose(conv_marg, 1, atol=err_marg), "Trans marginal not unity"


@pytest.mark.parametrize(
    "fp,tol",
    [(os.path.join(testdir, "data", "toy.symm.upper.2.cool"), 1e-1)],
)
def test_balancing_with_blacklist(temp_cool, fp, tol):
    runner = CliRunner()
    temp_fp = temp_cool(fp)

    blacklist_cases = [
        "chr1\t0\t2\n",  # Single row
        "chr1\t0\t2\nchr2\t0\t2\n",  # Multi-row
        "track=test\nchr1\t0\t2\n",  # With header
        "",  # Empty
    ]

    for bed_content in blacklist_cases:
        with tempfile.NamedTemporaryFile(mode="w", delete=False,
                                         suffix=".bed") as f:
            f.write(bed_content)
            blacklist_path = f.name

        try:
            args = [
                temp_fp, "--blacklist", blacklist_path, "--ignore-diags", "1",
                "--min-nnz", "1", "--tol", str(tol), "--force",
                "--nproc", "1", "--max-iters", "1000"
            ]
            result = runner.invoke(balance, args)
            err_msg = f"Command failed for {bed_content!r}: {result.output}"
            assert result.exit_code == 0, err_msg

            with h5py.File(temp_fp, "r") as h5:
                weights = h5["/bins/weight"][:]
                err_msg = f"Weights all NaN for {bed_content!r}"
                assert not np.all(np.isnan(weights)), err_msg

                if bed_content.strip():
                    clr = cooler.Cooler(temp_fp)
                    bins = clr.bins()[:]
                    mask = (bins["chrom"] == "chr1") & (bins["start"] >= 0)
                    masked = bins[mask & (bins["end"] <= 2)].index
                    err_msg = f"Blacklist not applied for {bed_content!r}"
                    assert np.all(np.isnan(weights[masked])), err_msg

        finally:
            os.unlink(blacklist_path)
