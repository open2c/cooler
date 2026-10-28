import os

import h5py
import numpy as np
import pytest

import cooler
from cooler import _balance

testdir = os.path.dirname(os.path.realpath(__file__))


@pytest.mark.parametrize(
    "fp,tol",
    [(os.path.join(testdir, "data", "hg19.GM12878-MboI.matrix.2000kb.cool"), 1e-2)],
)
def test_balancing_genomewide(fp, tol):
    clr = cooler.Cooler(fp)
    weights, _stats = _balance.balance_cooler(
        clr, ignore_diags=1, min_nnz=10, tol=tol
    )

    # Extract matrix and apply weights
    mat = cooler.Cooler(fp).matrix(balance=False, sparse=True)[:, :]
    mat.data = weights[mat.row] * weights[mat.col] * mat.data
    arr = mat.toarray()

    # Re-apply bin level filters
    mask = np.isnan(weights)
    arr[mask, :] = 0
    arr[:, mask] = 0

    # Apply diagonal filter
    np.fill_diagonal(arr, 0)

    # Check that the balanced marginal is flat
    marg = np.sum(arr, axis=0)
    var = np.var(marg[marg != 0])
    assert var < tol

    # Check that the balanced marginal is unity
    conv_marg = marg[~np.isnan(marg)].mean()
    err_marg = marg[~np.isnan(marg)].std()
    assert np.isclose(conv_marg, 1, atol=err_marg)


def test_marginalize_counts_diagonal_once():
    # Upper triangle of a symmetric 3x3 matrix, with diagonal entries
    pixels = {
        "bin1_id": np.array([0, 0, 1, 1, 2]),
        "bin2_id": np.array([0, 1, 1, 2, 2]),
    }
    data = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    chunk = {"bins": {"chrom": np.zeros(3, dtype=int)}, "pixels": pixels}
    dense = np.array([[1, 2, 0], [2, 3, 4], [0, 4, 5]])
    assert np.allclose(_balance._marginalize(chunk, data), dense.sum(axis=0))


@pytest.mark.parametrize(
    "fp,tol",
    [(os.path.join(testdir, "data", "hg19.GM12878-MboI.matrix.2000kb.cool"), 1e-5)],
)
def test_balancing_no_ignore_diags(fp, tol):
    # The main diagonal must be counted once per bin in the marginals
    clr = cooler.Cooler(fp)
    weights, _stats = _balance.balance_cooler(
        clr, ignore_diags=0, min_nnz=10, tol=tol, max_iters=1000
    )

    mat = cooler.Cooler(fp).matrix(balance=False, sparse=True)[:, :]
    mat.data = weights[mat.row] * weights[mat.col] * mat.data
    arr = mat.toarray()

    mask = np.isnan(weights)
    arr[mask, :] = 0
    arr[:, mask] = 0

    # Check that the balanced marginal, including the diagonal, is flat
    marg = np.sum(arr, axis=0)
    var = np.var(marg[marg != 0])
    assert var < tol


@pytest.mark.filterwarnings("ignore")
@pytest.mark.parametrize(
    "fp,tol",
    [(os.path.join(testdir, "data", "hg19.GM12878-MboI.matrix.2000kb.cool"), 1e-2)],
)
def test_balancing_cisonly(fp, tol):
    with h5py.File(fp, "r") as h5:
        clr = cooler.Cooler(h5)
        chrom_offsets = h5["indexes/chrom_offset"][:]
        weights, _stats = _balance.balance_cooler(
            clr, ignore_diags=1, min_nnz=10, tol=tol, cis_only=True
        )

    # Extract matrix and apply weights
    mat = cooler.Cooler(fp).matrix(balance=False, sparse=True)[:, :]
    mat.data = weights[mat.row] * weights[mat.col] * mat.data
    arr = mat.toarray()

    # Re-apply bin level filters
    mask = np.isnan(weights)
    arr[mask, :] = 0
    arr[:, mask] = 0

    # Apply diagonal filter
    np.fill_diagonal(arr, 0)

    # Filter out trans data
    spans = list(zip(chrom_offsets[:-1], chrom_offsets[1:]))
    from scipy.linalg import block_diag

    blocks = [np.ones((hi - lo,) * 2) for lo, hi in spans]
    mask = block_diag(*blocks).astype(bool)
    arr[~mask] = 0

    # Check that the balanced marginal is flat
    marg = np.sum(arr, axis=0)
    for lo, hi in spans:
        m = marg[lo:hi]
        m = m[m != 0]
        if len(m):
            print(lo, hi)
            var = np.var(m[m != 0])
            assert var < tol

            # Check that the balanced marginal is unity
            conv_marg = m[~np.isnan(m)].mean()
            err_marg = m[~np.isnan(m)].std()
            assert np.isclose(conv_marg, 1, atol=err_marg)


@pytest.mark.filterwarnings("ignore")
@pytest.mark.parametrize(
    "fp,tol",
    [(os.path.join(testdir, "data", "hg19.GM12878-MboI.matrix.2000kb.cool"), 1e-2)],
)
def test_balancing_transonly(fp, tol):
    with h5py.File(fp, "r") as h5:
        clr = cooler.Cooler(h5)
        chrom_offsets = h5["indexes/chrom_offset"][:]
        weights, _stats = _balance.balance_cooler(
            clr, ignore_diags=1, min_nnz=10, tol=tol, trans_only=True
        )

    # Extract matrix and apply weights
    mat = cooler.Cooler(fp).matrix(balance=False, sparse=True)[:, :]
    mat.data = weights[mat.row] * weights[mat.col] * mat.data
    arr = mat.toarray()

    # Re-apply bin level filters
    mask = np.isnan(weights)
    arr[mask, :] = 0
    arr[:, mask] = 0

    # Filter out cis data
    for lo, hi in zip(chrom_offsets[:-1], chrom_offsets[1:]):
        arr[lo:hi, lo:hi] = np.nan

    # Check that the balanced marginal is flat
    marg = np.nansum(arr, axis=0)
    var = np.nanvar(marg[marg != 0])
    assert var < tol

    # Check that the balanced marginal is unity
    conv_marg = marg[~np.isnan(marg)].mean()
    err_marg = marg[~np.isnan(marg)].std()
    assert np.isclose(conv_marg, 1, atol=err_marg)


@pytest.mark.parametrize(
    "fp,tol",
    [(os.path.join(testdir, "data", "hg19.GM12878-MboI.matrix.2000kb.cool"), 1e-5)],
)
def test_balancing_other_options(fp, tol):
    clr = cooler.Cooler(fp)
    x0 = np.random.default_rng(0).random(len(clr.bins())) + 0.5
    x0_copy = x0.copy()
    _weights, stats = _balance.balance_cooler(
        clr, ignore_diags=1, min_nnz=10, tol=tol, x0=x0
    )
    assert stats['converged']
    np.testing.assert_array_equal(x0, x0_copy)  # not modified

    _weights, stats = _balance.balance_cooler(
        clr,
        chunksize=3000,
        ignore_diags=1,
        min_nnz=10,
        tol=tol,
    )
    assert stats['converged']

    _weights, stats = _balance.balance_cooler(
        clr,
        ignore_diags=1,
        min_nnz=10,
        min_count=1,
        mad_max=12,
        tol=tol,
    )
    assert stats['converged']

    _weights, stats = _balance.balance_cooler(
        clr, ignore_diags=1, min_nnz=10, tol=tol, blacklist=[0, 4]
    )
    assert stats['converged']


def test_balancing_without_chunks():
    clr = cooler.Cooler(os.path.join(testdir, "data", "toy.symm.upper.2.cool"))
    w1, _ = cooler.balance_cooler(clr, chunksize=None, mad_max=0, min_nnz=0)
    w2, _ = cooler.balance_cooler(clr, chunksize=7, mad_max=0, min_nnz=0)
    np.testing.assert_allclose(w1, w2)
