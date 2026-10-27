"""
Tests for the commands that compute new data: merge, coarsen, zoomify and
balance.
"""
import shutil

import multiprocess as mp
import numpy as np
import pandas as pd
import pytest
from _common import DATA

import cooler

TOY = DATA / "toy.symm.upper.2.cool"
NOFILTER = ["--ignore-diags", 2, "--mad-max", 0, "--min-nnz", 0, "--tol", 0.05]


def test_merge(run_cli, tmp_path):
    out = tmp_path / "double.cool"
    run_cli("merge", out, TOY, TOY, "--field", "count:dtype=int")
    total = cooler.Cooler(str(TOY)).pixels()["count"][:].sum()
    assert cooler.Cooler(str(out)).pixels()["count"][:].sum() == 2 * total

    run_cli("merge", out, TOY, TOY, "--field", "count:agg=max")
    assert cooler.Cooler(str(out)).pixels()["count"][:].sum() == total
    run_cli("merge", out, TOY, TOY, "--field", "count:agg=mean", exit_code=2)
    run_cli("merge", out, exit_code=2)  # no inputs


def test_coarsen(run_cli, tmp_path):
    out = tmp_path / "coarse.cool"
    run_cli("coarsen", TOY, "--factor", 2, "--nproc", 2, "-o", out)
    expected = cooler.Cooler(str(DATA / "toy.symm.upper.4.cool")).pixels()[:]
    pd.testing.assert_frame_equal(cooler.Cooler(str(out)).pixels()[:], expected)

    run_cli(
        "coarsen", TOY, "--factor", 2, "--field", "count:dtype=float,agg=max",
        "-o", out,
    )
    assert cooler.Cooler(str(out)).pixels()["count"][:].dtype.kind == "f"

    result = run_cli(
        "coarsen", TOY, "--factor", 2, "--field", "count:agg=mean", "-o", out,
        exit_code=2,
    )
    assert "Unsupported aggregation 'mean'" in result.stderr


def resolutions(path):
    return sorted(int(uri.split("/")[-1]) for uri in cooler.fileops.list_coolers(path))


def test_zoomify(run_cli, tmp_path):
    out = tmp_path / "toy.mcool"
    run_cli("zoomify", TOY, "--balance", "--nproc", 2, "-o", out)
    assert resolutions(out)[0] == 2
    for res in resolutions(out):
        assert "weight" in cooler.Cooler(f"{out}::resolutions/{res}").bins()

    run_cli(
        "zoomify", TOY, "--resolutions", "2,4,8", "--field",
        "count:dtype=float,agg=max", "-o", out,
    )
    assert resolutions(out) == [2, 4, 8]
    clr = cooler.Cooler(f"{out}::resolutions/8")
    assert clr.pixels()["count"][:].dtype.kind == "f"

    # Without -o, the output is named after the input
    shutil.copy(TOY, tmp_path / "toy.cool")
    run_cli("zoomify", tmp_path / "toy.cool", "--resolutions", "2,4")
    assert resolutions(tmp_path / "toy.mcool") == [2, 4]


def test_zoomify_resolution_sequences(run_cli, tmp_path):
    # 2 Mb bins: a 256 x 256 tile holds the genome at about 12 Mb
    clr = DATA / "hg19.GM12878-MboI.matrix.2000kb.cool"
    out = tmp_path / "hg19.mcool"
    for spec, expected in [
        ("n", [2, 4, 10]),
        ("b", [2, 4, 8]),
        ("2000000n", [2, 4, 10]),
        ("2000000b", [2, 4, 8]),
        ("2000000,6000000", [2, 6]),
    ]:
        run_cli("zoomify", clr, "--resolutions", spec, "-o", out)
        assert resolutions(out) == [r * 1_000_000 for r in expected], spec

    # 4DN resolutions start at 1 kb, which can't be made from 2 Mb bins
    run_cli("zoomify", clr, "--resolutions", "4dn", "-o", out, exit_code=1)


def test_zoomify_variable_bins(run_cli, tmp_path):
    out = tmp_path / "var.mcool"
    run_cli("zoomify", DATA / "toy.symm.upper.var.cool", "-r", "1,2", "-o", out)
    assert resolutions(out) == [1, 2]


def weights(stdout):
    # NaN weights are written as blank lines
    return np.array([float(x) if x.strip() else np.nan for x in stdout.splitlines()])


def test_balance(run_cli, tmp_path):
    n_bins = cooler.Cooler(str(TOY)).info["nbins"]
    result = run_cli("balance", TOY, *NOFILTER, "--nproc", 2, "--stdout")
    assert len(result.stdout.splitlines()) == n_bins
    assert not mp.active_children()  # worker processes have exited

    # Convergence policies, without iterating
    no_iters = [*NOFILTER, "--max-iters", 0, "--stdout"]
    result = run_cli("balance", TOY, *no_iters, "--convergence-policy", "store_final")
    assert len(result.stdout.splitlines()) == n_bins
    result = run_cli("balance", TOY, *no_iters, "--convergence-policy", "store_nan")
    assert np.isnan(weights(result.stdout)).all()
    assert len(weights(result.stdout)) == n_bins
    result = run_cli("balance", TOY, *no_iters, "--convergence-policy", "discard")
    assert result.stdout == ""
    run_cli("balance", TOY, *no_iters, "--convergence-policy", "error", exit_code=1)

    run_cli("balance", TOY, "--cis-only", "--trans-only", exit_code=2)

    # Storing the weights in the file
    path = tmp_path / "toy.cool"
    shutil.copy(TOY, path)
    run_cli("balance", path, "--check", exit_code=1)
    run_cli("balance", path, *NOFILTER)
    assert "is balanced" in run_cli("balance", path, "--check").stdout
    run_cli("balance", path, *NOFILTER, exit_code=1)  # already exists
    run_cli("balance", path, *NOFILTER, "--force", "--name", "weight")
    run_cli("balance", path, *NOFILTER, "--name", "other")
    assert {"weight", "other"} <= set(cooler.Cooler(str(path)).bins().columns)


def test_balance_blacklist(run_cli, tmp_path):
    blacklist = tmp_path / "blacklist.bed"
    pd.DataFrame({"chrom": ["chr1", "chr2"], "start": [5, 10], "end": [10, 20]}).to_csv(
        blacklist, sep="\t", index=False, header=False
    )
    result = run_cli("balance", TOY, *NOFILTER, "--blacklist", blacklist, "--stdout")
    w = weights(result.stdout)
    assert np.isnan(w[2:5]).all()  # chr1:5-10 at 2 bp
    assert np.isnan(w[21:26]).all()  # chr2:10-20

    # A single line, header/track lines, or no records
    for content in [
        "chr1\t5\t10\n",
        "track name=bl\nchrom\tstart\tend\nchr1\t5\t10\n",
        "#chrom\tstart\tend\nchr1\t5\t10\n",
        "",
    ]:
        blacklist.write_text(content)
        run_cli("balance", TOY, *NOFILTER, "--blacklist", blacklist, "--stdout")


def test_zoomify_does_not_overwrite_its_input(run_cli, tmp_path):
    # The default output name of an .mcool or .h5 input is the input itself
    mcool = tmp_path / "t.mcool"
    shutil.copy(DATA / "toy.symm.upper.2.mcool", mcool)
    before = resolutions(mcool)
    run_cli("zoomify", "-r", "4,8", f"{mcool}::resolutions/2", exit_code=2)
    h5 = tmp_path / "t.h5"
    shutil.copy(TOY, h5)
    run_cli("zoomify", "-r", "4", h5, "-o", h5, exit_code=2)
    assert resolutions(mcool) == before
    assert cooler.fileops.is_cooler(h5)

    with pytest.raises(ValueError, match="contains the input"):
        cooler.zoomify_cooler(str(h5), str(h5), [2, 4], 100)


def test_zoomify_bad_arguments(run_cli, tmp_path):
    out = tmp_path / "toy.mcool"
    run_cli("zoomify", TOY, "-r", "2,abc", "-o", out, exit_code=2)
    run_cli("zoomify", TOY, "-r", "2,4", "-o", f"{out}::group", exit_code=2)
    result = run_cli("zoomify", TOY, "-r", "3", "-o", out, exit_code=1)
    assert "cannot be derived" in result.stderr
