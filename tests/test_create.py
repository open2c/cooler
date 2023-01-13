import os
import os.path as op
import tempfile
from io import BytesIO

import h5py
import numpy as np
import pandas as pd
import pytest
from _common import isolated_filesystem

import cooler
import cooler.create

tmp = tempfile.gettempdir()
testdir = op.dirname(op.realpath(__file__))
datadir = op.join(testdir, "data")


@pytest.mark.parametrize(
    "fp", [op.join(datadir, "hg19.GM12878-MboI.matrix.2000kb.cool")]
)
def test_create_append(fp):
    import dask.dataframe as dd

    c = cooler.Cooler(fp)
    # chromsizes = c.chromsizes
    bins = c.bins()[:]
    pixels = c.pixels()[:]

    # create
    cooler.create.create(op.join(tmp, "test.df.2000kb.cool"), bins, pixels)
    cooler.create.create(
        op.join(tmp, "test.dict.2000kb.cool"),
        bins,
        {k: v for k, v in pixels.items()},
    )
    cooler.create.create(op.join(tmp, "test.iter_df.2000kb.cool"), bins, [pixels])
    cooler.create.create(
        op.join(tmp, "test.iter_dict.2000kb.cool"),
        bins,
        [{k: v for k, v in pixels.items()}],
    )
    ddf = dd.from_pandas(pixels, npartitions=3)
    cooler.create.create(op.join(tmp, "test.ddf.2000kb.cool"), bins, ddf)

    # Append
    cooler.create.append(
        op.join(tmp, "test.df.2000kb.cool"),
        "bins",
        {"start_1based": bins.apply(lambda x: x.start + 1, axis=1)},
    )
    cooler.create.append(op.join(tmp, "test.df.2000kb.cool"), "bins", {"ones": 1})
    series = ddf["count"] / ddf["count"].sum()
    series.name = "normed"
    cooler.create.append(op.join(tmp, "test.df.2000kb.cool"), "pixels", series)
    cooler.create.append(
        op.join(tmp, "test.df.2000kb.cool"), "pixels", series, force=True
    )
    cooler.create.append(
        op.join(tmp, "test.df.2000kb.cool"),
        "bins",
        {"twos": [np.ones(1000, dtype=int) * 2, np.ones(561, dtype=int) * 2]},
        chunked=True,
        force=True,
    )
    c2 = cooler.Cooler(op.join(tmp, "test.df.2000kb.cool"))
    assert len(c2.bins().columns) == 6
    assert len(c2.pixels().columns) == 4


@pytest.mark.parametrize(
    "f_hm,f_cool",
    [
        (
            op.join(datadir, "hg19.IMR90-MboI.matrix.2000kb.npy"),
            op.join(tmp, "test.cool"),
        )
    ],
)
def test_roundtrip(f_hm, f_cool):
    chromsizes = cooler.read_chromsizes(
        "http://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes",
        name_patterns=(r"^chr[0-9]+$", r"chrX$"),
    )
    binsize = 2000000
    bintable = cooler.binnify(chromsizes, binsize)

    heatmap = np.load(f_hm)
    reader = cooler.create.ArrayLoader(bintable, heatmap, 100000)
    cooler.create.create(f_cool, bintable, reader, assembly="hg19")

    h5 = h5py.File(f_cool, "r")
    new_chromtable = cooler.api.chroms(h5)
    assert np.all(chromsizes.index == new_chromtable["name"])

    new_bintable = cooler.api.bins(h5)
    assert np.all(bintable == new_bintable)

    info = cooler.api.info(h5)
    assert info["genome-assembly"] == "hg19"
    assert info["bin-type"] == "fixed"
    assert info["bin-size"] == binsize

    mat = cooler.api.matrix(h5, 0, 100, 0, 100, "count", balance=False)
    assert mat.shape == (100, 100)
    assert np.allclose(heatmap[:100, :100], mat)

    mat = cooler.Cooler(h5).matrix("count", balance=False)[:100, :100]
    assert mat.shape == (100, 100)
    assert np.allclose(heatmap[:100, :100], mat)

    mat = cooler.api.matrix(h5, 100, 200, 100, 200, "count", balance=False)
    assert mat.shape == (100, 100)
    assert np.allclose(heatmap[100:200, 100:200], mat)

    mat = cooler.Cooler(h5).matrix("count", balance=False)[100:200, 100:200]
    assert mat.shape == (100, 100)
    assert np.allclose(heatmap[100:200, 100:200], mat)

    try:
        os.remove(f_cool)
    except OSError:
        pass


def test_rename_chroms():
    from shutil import copyfile

    with isolated_filesystem():
        copyfile(op.join(datadir, "toy.asymm.4.cool"), "toy.asymm.4.cool")
        clr = cooler.Cooler("toy.asymm.4.cool")
        assert clr.chromnames == ["chr1", "chr2"]
        cooler.rename_chroms(clr, {"chr1": "1", "chr2": "2"})
        assert clr.chromnames == ["1", "2"]  # the Cooler object is refreshed


def test_create_custom_cols():

    with isolated_filesystem():
        df = pd.DataFrame(
            {
                "bin1_id": [0, 1, 1, 1, 2, 2, 3, 4, 5],
                "bin2_id": [1, 1, 3, 4, 5, 6, 7, 8, 9],
                "foo": [1, 1, 1, 1, 1, 2, 2, 2, 2],
                "bar": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
            },
            columns=["bin1_id", "bin2_id", "foo", "bar"],
        )
        bins = pd.DataFrame(
            {
                "chrom": ["chr1"] * 5 + ["chr2"] * 5,
                "start": list(range(5)) * 2,
                "end": list(range(1, 6)) * 2,
            }
        )
        # works in unordered mode
        cooler.create_cooler("test.cool", bins, df, columns=["foo", "bar"])
        clr = cooler.Cooler("test.cool")
        assert len(clr.pixels().columns) == 4
        assert np.allclose(df, clr.pixels()[["bin1_id", "bin2_id", "foo", "bar"]][:])

        # works in ordered mode
        cooler.create_cooler(
            "test.cool", bins, df, columns=["foo", "bar"], ordered=True
        )
        clr = cooler.Cooler("test.cool")
        assert len(clr.pixels().columns) == 4
        assert np.allclose(df, clr.pixels()[["bin1_id", "bin2_id", "foo", "bar"]][:])

        # raises if no custom columns specified and 'count' does not exist
        with pytest.raises(ValueError):
            cooler.create_cooler("test.cool", bins, df, columns=None, ordered=True)


def test_write():
    chroms = pd.DataFrame({
        'name': ['chr1', 'chr2', 'chr3'],
        'length': [32, 48, 56],
    }, columns=['name', 'length'])
    bins = pd.DataFrame({
        'chrom': [
            'chr1', 'chr1',
            'chr2', 'chr2', 'chr2',
            'chr3', 'chr3', 'chr3', 'chr3',
        ],
        'start': [0, 16, 0, 16, 32, 0, 16, 32, 48],
        'end': [16, 32, 16, 32, 48, 16, 32, 48, 56],
    }, columns=['chrom', 'start', 'end'])
    pixels = pd.DataFrame({
        'bin1_id': [0],
        'bin2_id': [5],
        'count': [1],
    }, columns=['bin1_id', 'bin2_id', 'count'])
    h5opts = {'compression': 'gzip', 'compression_opts': 6}

    # chroms
    b = BytesIO()
    f = h5py.File(b, 'r+')
    grp = f.create_group('chroms')
    cooler.create._create.write_chroms(grp, chroms, h5opts)
    f.flush()
    assert 'name' in f['chroms']
    assert 'length' in f['chroms']

    # chroms with extra column
    b = BytesIO()
    f = h5py.File(b, 'r+')
    grp = f.create_group('chroms')
    cooler.create._create.write_chroms(grp, chroms.assign(foo=42), h5opts)
    f.flush()
    assert 'foo' in f['chroms']

    # bins
    b = BytesIO()
    f = h5py.File(b, 'r+')
    grp = f.create_group('bins')
    cooler.create._create.write_bins(
        grp, bins.assign(foo=42), ['chr1', 'chr2', 'chr3'], h5opts
    )
    f.flush()
    assert 'chrom' in f['bins']
    assert 'start' in f['bins']
    assert 'end' in f['bins']
    assert 'foo' in f['bins']

    # bins
    b = BytesIO()
    f = h5py.File(b, 'r+')
    grp = f.create_group('pixels')
    cooler.create._create.prepare_pixels(
        grp, len(bins), 1000000, ['count', 'foo'], {'foo': float}, h5opts
    )
    f.close()

    # pixels
    from multiprocessing import Lock
    nnz, total = cooler.create._create.write_pixels(
        b,
        'pixels',
        ['count', 'foo'],
        (pixels.assign(foo=42.0),),
        h5opts,
        Lock(),
    )
    b.seek(0)
    f = h5py.File(b, 'r')
    assert 'bin1_id' in f['pixels']
    assert 'bin2_id' in f['pixels']
    assert 'count' in f['pixels']
    assert 'foo' in f['pixels']
    assert f['pixels/foo'].dtype.kind == 'f'
    assert total == 1
    assert nnz == 1


def test_many_contigs():
    chroms = pd.DataFrame({
        'name': [f'scaffold_{i:05}' for i in range(4000)],
        'length': np.full(4000, 20),
    }, columns=['name', 'length'])
    bins = cooler.util.binnify(chroms.set_index('name')['length'], 10)
    h5opts = {'compression': 'gzip', 'compression_opts': 6}

    # chroms
    b = BytesIO()
    f = h5py.File(b, 'w')
    cooler.create._create.write_chroms(
        f.create_group('chroms'), chroms, h5opts
    )
    cooler.create._create.write_bins(
        f.create_group('bins'), bins, chroms['name'].values, h5opts
    )
    f.flush()
    assert 'enum_path' in f['bins/chrom'].attrs

    # TODO: make this more robust
    # assert f['bins/chrom'].attrs['enum_path'] == '/chroms/name'

    cooler.create._create._rename_chroms(
        f,
        {
            f'scaffold_{i:05}' : f'contig_{i:05}'
            for i in range(4000)
        },
        h5opts
    )


def test_create_cooler():
    chromsizes = cooler.util.read_chromsizes(
        op.join(datadir, "toy.chrom.sizes")
    )
    bins = cooler.util.binnify(chromsizes, 1)
    pixels = pd.read_csv(
        op.join(datadir, "toy.symm.upper.1.zb.coo"),
        sep='\t',
        names=['bin1_id', 'bin2_id', 'count']
    )
    pixels['foo'] = 42.0

    with isolated_filesystem():
        cooler.create.create_cooler(
            "test.cool",
            bins,
            pixels,
            assembly='toy',
            metadata={'hello': 'world', 'list': [1, 2, 3]},
        )

        cooler.create.create_cooler(
            "test.cool::foo/bar",
            bins,
            pixels,
        )

        cooler.create.create_cooler(
            "test.cool",
            bins,
            pixels,
            symmetric_upper=False
        )

        cooler.create.create_cooler(
            "test.cool",
            bins,
            pixels,
            columns=['count', 'foo'],
            dtypes={'foo': np.float64}
        )

        cooler.create.create_cooler(
            "test.cool",
            bins,
            pixels.to_dict(orient='series'),
        )

        cooler.create.create_cooler(
            "test.cool",
            bins,
            (pixels,),
        )

        cooler.create.create_cooler(
            "test.cool",
            bins,
            (pixels.to_dict(orient='series'),),
        )

        two_piece = (
            pixels.iloc[:len(pixels) // 2], pixels.iloc[len(pixels) // 2:]
        )
        cooler.create.create_cooler(
            "test.cool",
            bins,
            two_piece,
            ordered=True
        )
        cooler.create.create_cooler(
            "test.cool",
            bins,
            two_piece[::-1],
            ordered=False
        )

        many_piece = tuple(
            pixels.iloc[lo:hi] for lo, hi in
            cooler.util.partition(0, len(pixels), 5)
        )[::-1]
        cooler.create.create_cooler(
            "test.cool",
            bins,
            many_piece,
            ordered=False,
            max_merge=10
        )

        with pytest.raises(ValueError):
            cooler.create.create_cooler(
                "test.cool",
                bins,
                pixels,
                columns=['count', 'missing'],
            )

        with pytest.raises(ValueError):
            cooler.create.create_cooler(
                "test.cool",
                bins[['start', 'end']],
                pixels,
                columns=['count', 'missing'],
            )

        with pytest.raises(ValueError):
            cooler.create.create_cooler(
                "test.cool",
                bins[['start', 'end']],
                pixels,
                h5opts={'shuffuffle': 'boing'}
            )


def test_create_cooler_from_dask():
    dd = pytest.importorskip("dask.dataframe")

    chromsizes = cooler.util.read_chromsizes(
        op.join(datadir, "toy.chrom.sizes")
    )
    bins = cooler.util.binnify(chromsizes, 1)
    pixels = pd.read_csv(
        op.join(datadir, "toy.symm.upper.1.zb.coo"),
        sep='\t',
        names=['bin1_id', 'bin2_id', 'count']
    )
    pixels = dd.from_pandas(pixels, npartitions=10)

    with isolated_filesystem():
        cooler.create.create_cooler(
            "test.cool",
            bins,
            pixels,
            ordered=True
        )

        # TODO: unordered with dask is broken...
        # cooler.create.create_cooler(
        #     "test.cool",
        #     bins,
        #     pixels,
        #     ordered=False
        # )


@pytest.mark.parametrize(
    "fp", [op.join(datadir, "hg19.GM12878-MboI.matrix.2000kb.cool")]
)
def test_create_scool(fp):
    c = cooler.Cooler(fp)
    # chromsizes = c.chromsizes
    bins = c.bins()[:]
    pixels = c.pixels()[:]

    # random and different content to prove only chrom, start, end is linked and the rest is independent for each cell
    from copy import deepcopy
    bins_cell1 = deepcopy(bins)
    bins_cell2 = deepcopy(bins)
    bins_cell3 = deepcopy(bins)
    bins_cell1['weight'] = np.array([0] * len(bins_cell1["start"]))
    bins_cell2['weight'] = np.array([1] * len(bins_cell1["start"]))
    bins_cell3['weight'] = np.array([2] * len(bins_cell1["start"]))

    bins_cell1['KR'] = np.array([3] * len(bins_cell1["start"]))
    bins_cell2['KR'] = np.array([4] * len(bins_cell1["start"]))
    bins_cell3['KR'] = np.array([5] * len(bins_cell1["start"]))

    name_pixel_dict = {'cell1': pixels, 'cell2': pixels, 'cell3': pixels}
    name_bins_dict = {'cell1': bins_cell1, 'cell2': bins_cell2, 'cell3': bins_cell3}

    with isolated_filesystem():
        cooler.create_scool('outfile_test.scool', name_bins_dict, name_pixel_dict)
        content_of_scool = cooler.fileops.list_scool_cells('outfile_test.scool')
        content_expected = ['/cells/cell1', '/cells/cell2', '/cells/cell3']
        for content in content_expected:
            assert content in content_of_scool

        cooler.create_scool('outfile_test.scool', bins, name_pixel_dict)
        content_of_scool = cooler.fileops.list_scool_cells('outfile_test.scool')
        content_expected = ['/cells/cell1', '/cells/cell2', '/cells/cell3']
        for content in content_expected:
            assert content in content_of_scool
