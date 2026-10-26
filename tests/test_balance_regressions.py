"""
Regression tests for bugs in matrix balancing.
"""
import h5py
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

import cooler
from cooler._balance import balance_cooler
from cooler.cli.balance import balance as balance_cli

CHROMSIZES = pd.Series({"chr1": 300, "chr2": 200, "chr3": 100})


def make_cooler(path, bins, seed=0, density=0.5, **kwargs):
    """A random symmetric cooler, stored as its upper triangle by default."""
    rng = np.random.default_rng(seed)
    n = len(bins)
    mat = np.triu(rng.integers(1, 20, (n, n)) * (rng.random((n, n)) < density))
    i, j = np.nonzero(mat)
    if kwargs.get("symmetric_upper", True):
        px = pd.DataFrame({"bin1_id": i, "bin2_id": j, "count": mat[i, j]})
    else:
        full = mat + np.triu(mat, 1).T
        i, j = np.nonzero(full)
        px = pd.DataFrame({"bin1_id": i, "bin2_id": j, "count": full[i, j]})
    cooler.create_cooler(path, bins, px, **kwargs)
    return path


def dense(clr):
    return clr.matrix(balance=False)[:].astype(float)


def dense_balance(mat, tol=1e-14, max_iters=10_000):
    """Reference: balance a dense symmetric matrix, dropping empty bins."""
    bias = (mat.sum(axis=0) > 0).astype(float)
    for _ in range(max_iters):
        marg = (bias[:, None] * mat * bias[None, :]).sum(axis=0)
        nz = marg != 0
        if marg[nz].var() / marg[nz].mean() ** 2 < tol:
            break
        bias[nz] /= marg[nz] / marg[nz].mean()
    marg = (bias[:, None] * mat * bias[None, :]).sum(axis=0)
    bias[marg == 0] = np.nan
    return bias / np.sqrt(np.nanmean(marg[marg != 0]))


def cis_mask(bins):
    chroms = bins["chrom"].to_numpy()
    return chroms[:, None] == chroms[None, :]


NOFILTER = {"mad_max": 0, "min_nnz": 0, "tol": 1e-12}


def test_square_storage(tmp_path):
    bins = cooler.binnify(CHROMSIZES, 20)
    square = make_cooler(
        str(tmp_path / "s.cool"), bins, symmetric_upper=False, triucheck=False
    )
    assert cooler.Cooler(square).storage_mode == "square"
    # Not defined for a possibly asymmetric matrix; it would need separate row
    # and column weights
    with pytest.raises(ValueError, match="symmetric-upper"):
        balance_cooler(cooler.Cooler(square), **NOFILTER)


def test_ignore_diags_keeps_trans_pixels(tmp_path):
    # Bins next to each other across a chromosome boundary are trans (#211)
    bins = cooler.binnify(CHROMSIZES, 20)
    clr = cooler.Cooler(make_cooler(str(tmp_path / "c.cool"), bins))
    mat = dense(clr)
    i, j = np.indices(mat.shape)
    mat[(np.abs(i - j) < 3) & cis_mask(bins)] = 0

    weights, _ = balance_cooler(clr, ignore_diags=3, **NOFILTER)
    np.testing.assert_allclose(weights, dense_balance(mat), rtol=1e-5)


def test_ignore_dist(tmp_path):
    # Uniform bins: the same as the equivalent number of diagonals
    bins = cooler.binnify(CHROMSIZES, 20)
    clr = cooler.Cooler(make_cooler(str(tmp_path / "u.cool"), bins))
    w_diags, _ = balance_cooler(clr, ignore_diags=3, **NOFILTER)
    w_dist, stats = balance_cooler(clr, ignore_diags=0, ignore_dist=41, **NOFILTER)
    np.testing.assert_allclose(w_dist, w_diags)
    assert stats["ignore_dist"] == 41

    # Variable-size bins (#464)
    rng = np.random.default_rng(1)
    parts = []
    for chrom, size in CHROMSIZES.items():
        edges = np.unique(np.r_[0, rng.choice(np.arange(1, size), 9), size])
        parts.append(pd.DataFrame(
            {"chrom": chrom, "start": edges[:-1], "end": edges[1:]}
        ))
    vbins = pd.concat(parts, ignore_index=True)
    clr = cooler.Cooler(make_cooler(str(tmp_path / "v.cool"), vbins))
    mat = dense(clr)
    mids = (vbins["start"] + vbins["end"]).to_numpy() / 2
    mat[(np.abs(mids[:, None] - mids[None, :]) < 50) & cis_mask(vbins)] = 0
    weights, _ = balance_cooler(clr, ignore_diags=0, ignore_dist=50, **NOFILTER)
    np.testing.assert_allclose(weights, dense_balance(mat), rtol=1e-5)

    result = CliRunner().invoke(
        balance_cli,
        [clr.uri, "--ignore-dist", "50", "--mad-max", "0", "--min-nnz", "0"],
    )
    assert result.exit_code == 0, result.output


def test_bins_without_data_get_nan(tmp_path):
    bins = cooler.binnify(CHROMSIZES, 20)
    path = str(tmp_path / "c.cool")
    rng = np.random.default_rng(0)
    # No contacts on chr3 at all, and none for bin 3 of chr1
    px = pd.DataFrame({
        "bin1_id": rng.integers(0, 25, 200), "bin2_id": rng.integers(0, 25, 200)
    })
    px = px[(px.bin1_id != 3) & (px.bin2_id != 3)]
    px[["bin1_id", "bin2_id"]] = np.sort(px[["bin1_id", "bin2_id"]].to_numpy(), 1)
    px = px.groupby(["bin1_id", "bin2_id"]).size().rename("count").reset_index()
    cooler.create_cooler(path, bins, px)
    clr = cooler.Cooler(path)
    empty = np.r_[3, np.arange(25, 30)]

    for kwargs in [
        {"mad_max": 0, "min_nnz": 0},  # no filtering at all
        {"mad_max": 5, "min_nnz": 0},  # MAD filter with an empty chromosome
    ]:
        weights, stats = balance_cooler(clr, ignore_diags=0, **kwargs)
        assert np.isnan(weights[empty]).all()
        assert np.isfinite(np.delete(weights, empty)).all()
        assert stats["converged"]

    # An empty chromosome on its own counts as converged in cis-only mode
    weights, stats = balance_cooler(clr, ignore_diags=0, cis_only=True, **NOFILTER)
    assert np.isnan(weights[empty]).all()
    assert np.all(stats["converged"])


def test_every_bin_filtered(tmp_path):
    bins = cooler.binnify(CHROMSIZES, 20)
    clr = cooler.Cooler(make_cooler(str(tmp_path / "c.cool"), bins))
    for kwargs in [{}, {"cis_only": True}, {"trans_only": True}]:
        weights, stats = balance_cooler(clr, min_nnz=1000, **kwargs)
        assert np.isnan(weights).all()
        assert not np.any(stats["converged"])


def test_max_iters_and_x0(tmp_path):
    bins = cooler.binnify(CHROMSIZES, 20)
    clr = cooler.Cooler(make_cooler(str(tmp_path / "c.cool"), bins))

    with pytest.warns(UserWarning, match="Iteration limit"):
        weights, stats = balance_cooler(clr, max_iters=0, **NOFILTER)
    assert not stats["converged"]
    np.testing.assert_array_equal(weights, 1.0)  # filtered, but not balanced

    x0 = np.ones(len(bins), dtype=int)
    _weights, stats = balance_cooler(clr, x0=x0, **NOFILTER)
    assert stats["converged"]
    assert x0.dtype == int
    assert (x0 == 1).all()
    with pytest.raises(ValueError, match="x0"):
        balance_cooler(clr, x0=np.ones(3))


def test_large_stats_attributes(tmp_path):
    # Per-chromosome stats of cis-only balancing with many contigs are larger
    # than the 64 KiB limit of compact HDF5 attributes (#226)
    from cooler.create._writer import _add_columns

    bins = cooler.binnify(CHROMSIZES, 20)
    path = make_cooler(str(tmp_path / "c.cool"), bins)
    stats = {"scale": np.random.rand(10_000), "converged": np.ones(10_000, bool)}
    _add_columns(path, "bins", {"weight": np.ones(len(bins))}, attrs={"weight": stats})
    with h5py.File(path, "r") as f:
        attrs = f["bins/weight"].attrs
        np.testing.assert_array_equal(attrs["scale"], stats["scale"])
        np.testing.assert_array_equal(attrs["converged"], stats["converged"])


def test_cli_blacklist_unknown_chromosome(tmp_path):
    bins = cooler.binnify(CHROMSIZES, 20)
    path = make_cooler(str(tmp_path / "c.cool"), bins)
    bed = tmp_path / "bl.bed"
    bed.write_text("chr1\t0\t40\nchrY\t0\t100\nchrM\t0\t10\n")
    result = CliRunner().invoke(
        balance_cli,
        [path, "--blacklist", str(bed), "--mad-max", "0", "--min-nnz", "0"],
    )
    assert result.exit_code == 0, result.output
    weights = cooler.Cooler(path).bins()["weight"][:]
    assert weights[:2].isna().all()
    assert weights[2:].notna().all()


def test_cli_rescale_marginals(tmp_path):
    # --no-rescale-marginals (#448)
    bins = cooler.binnify(CHROMSIZES, 20)
    path = make_cooler(str(tmp_path / "c.cool"), bins)
    args = [path, "--mad-max", "0", "--min-nnz", "0", "--force"]
    assert CliRunner().invoke(balance_cli, args).exit_code == 0
    rescaled = cooler.Cooler(path).bins()["weight"][:].to_numpy()
    result = CliRunner().invoke(balance_cli, [*args, "--no-rescale-marginals"])
    assert result.exit_code == 0, result.output
    with h5py.File(path, "r") as f:
        raw = f["bins/weight"][:]
        scale = f["bins/weight"].attrs["scale"]
    np.testing.assert_allclose(raw, rescaled * np.sqrt(scale))
