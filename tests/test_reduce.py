import os.path as op

import h5py
import numpy as np
import pytest
from _common import cooler_cmp, isolated_filesystem

import cooler
from cooler.reduce import coarsen_cooler, legacy_zoomify, merge_coolers, zoomify_cooler

testdir = op.realpath(op.dirname(__file__))
datadir = op.join(testdir, "data")


@pytest.mark.parametrize(
    "path1,path2",
    [
        (
            op.join(datadir, "hg19.GM12878-MboI.matrix.2000kb.cool"),
            op.join(datadir, "hg19.GM12878-MboI.matrix.2000kb.cool"),
        ),
        (
            op.join(datadir, "toy.asymm.2.cool"),
            op.join(datadir, "toy.asymm.2.cool")
        )
    ],
)
def test_merge(path1, path2):
    with isolated_filesystem():
        merge_coolers("test.cool", [path1, path2], mergebuf=int(15e6))
        single = cooler.Cooler(path1)
        merged = cooler.Cooler("test.cool")
        assert (
            merged.pixels()["count"][:].sum() == 2 * single.pixels()["count"][:].sum()
        )


def test_merge2():
    with isolated_filesystem():
        path1 = op.join(datadir, "toy.symm.upper.2.cool")
        path2 = op.join(datadir, "toy.symm.upper.2.cool")
        merge_coolers(
            "test.cool", [path1, path2], mergebuf=3, agg={'count': 'mean'}
        )
        single = cooler.Cooler(path1)
        merged = cooler.Cooler("test.cool")
        assert (
            merged.pixels()["count"][:].sum() == single.pixels()["count"][:].sum()
        )

        # different resolution
        path1 = op.join(datadir, "toy.symm.upper.2.cool")
        path2 = op.join(datadir, "toy.symm.upper.4.cool")
        with pytest.raises(ValueError):
            merge_coolers("test.cool", [path1, path2], mergebuf=3)

        # incompatible bins
        path1 = op.join(datadir, "toy.symm.upper.var.cool")
        path2 = op.join(datadir, "toy.symm.upper.2.cool")
        with pytest.raises(ValueError):
            merge_coolers("test.cool", [path1, path2], mergebuf=3)

        path2 = op.join(datadir, "toy.symm.upper.var.cool")
        path1 = op.join(datadir, "toy.symm.upper.2.cool")
        with pytest.raises(ValueError):
            merge_coolers("test.cool", [path1, path2], mergebuf=3)

        # incompatible symmetry
        path1 = op.join(datadir, "toy.symm.upper.2.cool")
        path2 = op.join(datadir, "toy.asymm.2.cool")
        with pytest.raises(ValueError):
            merge_coolers("test.cool", [path1, path2], mergebuf=3)

        # missing value column
        path1 = op.join(datadir, "toy.symm.upper.2.cool")
        path2 = op.join(datadir, "toy.symm.upper.2.cool")
        with pytest.raises(ValueError):
            merge_coolers(
                "test.cool", [path1, path2], mergebuf=3, columns=["missing"]
            )


@pytest.mark.parametrize(
    "input_uri,factor,ref_uri",
    [
        (
            op.join(datadir, "toy.symm.upper.2.cool"),
            2,
            op.join(datadir, "toy.symm.upper.4.cool"),
        ),
        (
            op.join(datadir, "toy.asymm.2.cool"),
            2,
            op.join(datadir, "toy.asymm.4.cool")
        ),
        (
            op.join(datadir, "toy.symm.upper.var.cool"),
            2,
            op.join(datadir, "toy.symm.upper.var2x.cool"),
        ),
    ],
)
def test_coarsen(input_uri, factor, ref_uri):

    with isolated_filesystem():
        kwargs = dict(
            chunksize=10, nproc=1, columns=None, dtypes=None, agg=None
        )
        coarsen_cooler(input_uri, "test.cool", factor, **kwargs)
        cooler_cmp("test.cool", ref_uri)

        # custom dtype
        kwargs = dict(
            chunksize=10, nproc=1, columns=None, dtypes={'count': np.float64}
        )
        coarsen_cooler(input_uri, "test.cool", factor, **kwargs)
        with h5py.File('test.cool', 'r') as f:
            assert f['pixels/count'].dtype.kind == 'f'

        # custom aggregator
        kwargs = dict(
            chunksize=10, nproc=1, columns=None, dtypes=None, agg={'count': 'mean'}
        )
        coarsen_cooler(input_uri, "test.cool", factor, **kwargs)

        # parallel
        kwargs = dict(
            chunksize=10, nproc=2, columns=None, dtypes=None, agg=None
        )
        coarsen_cooler(input_uri, "test.cool", factor, **kwargs)

        # raise on missing value column
        kwargs = dict(
            chunksize=10, nproc=2, columns=['missing'], dtypes=None, agg=None
        )
        with pytest.raises(ValueError):
            coarsen_cooler(input_uri, "test.cool", factor, **kwargs)


def test_coarsen_partitions_correctly():
    kwargs = dict(nproc=1, columns=None, dtypes=None, agg=None)
    with isolated_filesystem():
        f_ref = op.join(datadir, "odd.4.cool")
        f_in = op.join(datadir, "odd.1.cool")
        coarsen_cooler(f_in, "odd.1.coarsen_4.cool", factor=4, chunksize=2, **kwargs)
        pix1 = cooler.Cooler(f_ref).pixels()[:]
        pix2 = cooler.Cooler("odd.1.coarsen_4.cool").pixels()[:]
        assert len(pix1) == len(pix2)
        assert sum(pix2[["bin1_id", "bin2_id"]].duplicated()) == 0
        assert np.allclose(pix1, pix2)


def test_zoomify():
    kwargs = dict(chunksize=10, nproc=1, columns=None, dtypes=None, agg=None)
    with isolated_filesystem():
        zoomify_cooler(
            op.join(datadir, "toy.asymm.2.cool"),
            "test.2.mcool",
            resolutions=[4, 8, 16, 32],
            **kwargs
        )
        for res in [2, 4, 8, 16, 32]:
            cooler_cmp(
                f"test.2.mcool::resolutions/{res}",
                op.join(datadir, f"toy.asymm.{res}.cool"),
            )

        # include base resolution
        zoomify_cooler(
            op.join(datadir, "toy.asymm.2.cool"),
            "test.2.mcool",
            resolutions=[2, 4, 8, 16, 32],
            **kwargs
        )
        for res in [2, 4, 8, 16, 32]:
            cooler_cmp(
                f"test.2.mcool::resolutions/{res}",
                op.join(datadir, f"toy.asymm.{res}.cool"),
            )

        # impossible resolution to obtain
        with pytest.raises(ValueError):
            zoomify_cooler(
                op.join(datadir, "toy.asymm.2.cool"),
                "test.2.mcool",
                resolutions=[4, 5, 32],
                **kwargs
            )


def test_legacy_zoomify():
    infile = op.join(datadir, "hg19.GM12878-MboI.matrix.2000kb.cool")
    chunksize = int(10e6)
    # n_zooms = 2
    n_cpus = 1
    with isolated_filesystem():
        legacy_zoomify(infile, "test.multires.cool", n_cpus, chunksize)


def test_append_mode():
    # merge
    path1 = path2 = op.join(datadir, "toy.asymm.2.cool")
    out_path = "test.cool"
    for append in (True, False):
        with isolated_filesystem():
            with h5py.File("test.cool", "w") as f:
                f.attrs["xxxx"] = True
            merge_coolers(
                out_path,
                [path1, path2],
                mergebuf=int(15e6),
                mode="a" if append else "w"
            )
            with h5py.File(out_path, "r") as f:
                if append:
                    assert "xxxx" in f.attrs
                else:
                    assert "xxxx" not in f.attrs

    # coarsen
    input_path = op.join(datadir, "toy.symm.upper.2.cool")
    out_path = "test.cool"
    for append in (True, False):
        with isolated_filesystem():
            with h5py.File("test.cool", "w") as f:
                f.attrs["xxxx"] = True
            coarsen_cooler(
                input_path,
                "test.cool",
                2,
                chunksize=10,
                nproc=1,
                columns=None,
                dtypes=None,
                agg=None,
                mode="a" if append else "w"
            )
            with h5py.File(out_path, "r") as f:
                if append:
                    assert "xxxx" in f.attrs
                else:
                    assert "xxxx" not in f.attrs
