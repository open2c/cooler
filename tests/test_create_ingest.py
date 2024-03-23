# import filecmp
import importlib.util
import os
import os.path as op
import tempfile

import h5py
import numpy as np
import pandas as pd
import pytest
from pandas.api import types

import cooler
from cooler.cli.cload import pairix as cload_pairix
from cooler.cli.cload import pairs as cload_pairs
from cooler.cli.cload import tabix as cload_tabix
from cooler.cli.load import load

pysam_missing = importlib.util.find_spec("pysam") is None
pairix_missing = importlib.util.find_spec("pypairix") is None
_pandas_major_version = int(pd.__version__.split(".")[0])

tmp = tempfile.gettempdir()
testdir = op.realpath(op.dirname(__file__))
testcool_path = op.join(tmp, "test.cool")


def _alternating_bins(chromsizes, steps):
    def _each(chrom):
        clen = chromsizes[chrom]
        edges = [0]
        total = 0
        i = 0
        while edges[i] < clen:
            step = steps[i % len(steps)]
            total += step
            i += 1
            edges.append(total)
            print(edges)
        edges[-1] = clen
        return pd.DataFrame(
            {"chrom": chrom, "start": edges[:-1], "end": edges[1:]},
            columns=["chrom", "start", "end"],
        )

    return pd.concat(map(_each, chromsizes.keys()), axis=0, ignore_index=True)


def test_from_hdf5_pairs():
    def should_not_depend_on_chunksize(chromsizes, bintable, mock_pairs):
        # try different chunk sizes
        binner = cooler.create.HDF5Aggregator(
            mock_pairs, chromsizes, bintable, chunksize=66
        )
        cooler.create.create(testcool_path, bintable, binner)
        with h5py.File(testcool_path, "r") as h5:
            oc1 = h5["indexes"]["chrom_offset"][:]
            ob1 = h5["indexes"]["bin1_offset"][:]
            p1 = cooler.api.pixels(h5, join=False)

        binner = cooler.create.HDF5Aggregator(
            mock_pairs, chromsizes, bintable, chunksize=666
        )
        cooler.create.create(testcool_path, bintable, binner)
        with h5py.File(testcool_path, "r") as h5:
            oc2 = h5["indexes"]["chrom_offset"][:]
            ob2 = h5["indexes"]["bin1_offset"][:]
            p2 = cooler.api.pixels(h5, join=False)

        assert np.all(oc1 == oc2)
        assert np.all(ob1 == ob2)
        assert np.all(p1.values == p2.values)

    def should_raise_if_input_not_sorted(chromsizes, bintable, mock_pairs):
        # not sorted by chrm1
        # with h5py.File(testcool_path, 'w') as h5:
        bad_reads = {
            "chrms1": mock_pairs["chrms2"],
            "cuts1": mock_pairs["cuts2"],
            "chrms2": mock_pairs["chrms1"],
            "cuts2": mock_pairs["cuts1"],
        }
        with pytest.raises(ValueError):
            cooler.create.HDF5Aggregator(bad_reads, chromsizes, bintable, chunksize=66)

        # not triu
        bad_reads = {
            "chrms1": mock_pairs["chrms1"].copy(),
            "cuts1": mock_pairs["cuts1"].copy(),
            "chrms2": mock_pairs["chrms2"].copy(),
            "cuts2": mock_pairs["cuts2"].copy(),
        }
        bad_reads["chrms1"][0] = 0
        bad_reads["chrms2"][0] = 0
        bad_reads["cuts1"][0] = 10
        bad_reads["cuts2"][0] = 9
        binner = cooler.create.HDF5Aggregator(
            bad_reads, chromsizes, bintable, chunksize=66
        )
        with pytest.raises(ValueError):
            cooler.create.create(testcool_path, bintable, binner)

    def should_work_with_int32_cols(chromsizes, bintable, mock_pairs):
        # int64
        binner = cooler.create.HDF5Aggregator(
            mock_pairs, chromsizes, bintable, chunksize=66
        )
        cooler.create.create(testcool_path, bintable, binner)
        with h5py.File(testcool_path, "r") as h5:
            oc1 = h5["indexes"]["chrom_offset"][:]
            ob1 = h5["indexes"]["bin1_offset"][:]
            p1 = cooler.api.pixels(h5, join=False)

        # int32
        mock_pairs32 = {
            "chrms1": mock_pairs["chrms1"].astype(np.int32),
            "cuts1": mock_pairs["cuts1"].astype(np.int32),
            "chrms2": mock_pairs["chrms2"].astype(np.int32),
            "cuts2": mock_pairs["cuts2"].astype(np.int32),
        }
        binner = cooler.create.HDF5Aggregator(
            mock_pairs32, chromsizes, bintable, chunksize=66
        )
        cooler.create.create(testcool_path, bintable, binner)
        with h5py.File(testcool_path, "r") as h5:
            oc2 = h5["indexes"]["chrom_offset"][:]
            ob2 = h5["indexes"]["bin1_offset"][:]
            p2 = cooler.api.pixels(h5, join=False)

        assert np.all(oc1 == oc2)
        assert np.all(ob1 == ob2)
        assert np.all(p1.values == p2.values)

    def _mock_hdf5_pairs():
        np.random.seed(1)
        chrms = np.random.randint(0, n_chroms, n_records * 2)
        cuts = np.random.randint(0, clen, n_records * 2)
        abs_cuts = np.array([clen * chrm + cut for chrm, cut in zip(chrms, cuts)])
        abs_cuts1, abs_cuts2 = abs_cuts[:n_records], abs_cuts[n_records:]
        mock_pairs = {
            "chrms1": chrms[:n_records],
            "cuts1": cuts[:n_records],
            "chrms2": chrms[n_records:],
            "cuts2": cuts[n_records:],
        }
        # Triu-sort
        mask = abs_cuts1 > abs_cuts2
        mock_pairs["chrms1"][mask], mock_pairs["chrms2"][mask] = (
            mock_pairs["chrms2"][mask],
            mock_pairs["chrms1"][mask],
        )
        mock_pairs["cuts1"][mask], mock_pairs["cuts2"][mask] = (
            mock_pairs["cuts2"][mask],
            mock_pairs["cuts1"][mask],
        )
        abs_cuts1[mask], abs_cuts2[mask] = abs_cuts2[mask], abs_cuts1[mask]
        idx = np.lexsort([abs_cuts2, abs_cuts1])
        for key in mock_pairs:
            mock_pairs[key] = mock_pairs[key][idx]
        return mock_pairs

    n_chroms = 2
    clen = 2000
    n_records = 3000
    chromsizes = pd.Series(index=["chr1", "chr2"], data=[clen, clen])
    mock_pairs = _mock_hdf5_pairs()

    # uniform bins
    bintable = cooler.binnify(chromsizes, 100)
    should_not_depend_on_chunksize(chromsizes, bintable, mock_pairs)
    should_raise_if_input_not_sorted(chromsizes, bintable, mock_pairs)
    should_work_with_int32_cols(chromsizes, bintable, mock_pairs)

    # non-uniform bins
    bintable = _alternating_bins(chromsizes, [10, 100])
    should_not_depend_on_chunksize(chromsizes, bintable, mock_pairs)
    should_raise_if_input_not_sorted(chromsizes, bintable, mock_pairs)
    should_work_with_int32_cols(chromsizes, bintable, mock_pairs)


@pytest.mark.skipif(pysam_missing, reason="pysam not installed")
@pytest.mark.filterwarnings("ignore")
@pytest.mark.parametrize(
    "bins_path,pairs_path,ref_path,nproc",
    [
        (
            op.join(testdir, "data", "hg19.bins.2000kb.bed.gz"),
            op.join(testdir, "data", "hg19.GM12878-MboI.pairs.subsample.sorted.txt.gz"),
            op.join(testdir, "data", "hg19.GM12878-MboI.matrix.2000kb.cool"),
            8,
        ),
        (
            op.join(testdir, "data", "hg19.bins.2000kb.bed.gz"),
            op.join(testdir, "data", "hg19.GM12878-MboI.pairs.subsample.sorted.txt.gz"),
            op.join(testdir, "data", "hg19.GM12878-MboI.matrix.2000kb.cool"),
            1,
        ),
    ],
)
def test_cload_tabix(bins_path, pairs_path, ref_path, nproc):
    cload_tabix.callback(
        bins_path,
        pairs_path,
        testcool_path,
        metadata=None,
        assembly="hg19",
        nproc=nproc,
        zero_based=False,
        max_split=2,
    )
    with h5py.File(testcool_path, "r") as f1, h5py.File(ref_path, "r") as f2:
        assert np.all(f1["pixels/bin1_id"][:] == f2["pixels/bin1_id"][:])
        assert np.all(f1["pixels/bin2_id"][:] == f2["pixels/bin2_id"][:])
        assert np.all(f1["pixels/count"][:] == f2["pixels/count"][:])
    try:
        os.remove(testcool_path)
    except OSError:
        pass


@pytest.mark.skipif(pairix_missing, reason="pairix not installed")
@pytest.mark.parametrize(
    "bins_path,pairs_path,ref_path,nproc",
    [
        (
            op.join(testdir, "data", "hg19.bins.2000kb.bed.gz"),
            op.join(testdir, "data", "hg19.GM12878-MboI.pairs.subsample.blksrt.txt.gz"),
            op.join(testdir, "data", "hg19.GM12878-MboI.matrix.2000kb.cool"),
            8,
        ),
        (
            op.join(testdir, "data", "hg19.bins.2000kb.bed.gz"),
            op.join(testdir, "data", "hg19.GM12878-MboI.pairs.subsample.blksrt.txt.gz"),
            op.join(testdir, "data", "hg19.GM12878-MboI.matrix.2000kb.cool"),
            1,
        ),
    ],
)
def test_cload_pairix(bins_path, pairs_path, ref_path, nproc):
    cload_pairix.callback(
        bins_path,
        pairs_path,
        testcool_path,
        metadata=None,
        assembly="hg19",
        nproc=nproc,
        zero_based=False,
        max_split=2,
        block_char="|",
    )
    with h5py.File(testcool_path, "r") as f1, h5py.File(ref_path, "r") as f2:
        assert np.all(f1["pixels/bin1_id"][:] == f2["pixels/bin1_id"][:])
        assert np.all(f1["pixels/bin2_id"][:] == f2["pixels/bin2_id"][:])
        assert np.all(f1["pixels/count"][:] == f2["pixels/count"][:])
    try:
        os.remove(testcool_path)
    except OSError:
        pass


@pytest.mark.skipif(pairix_missing, reason="pairix not installed")
def test_cload_pairix_wrong_block_char():
    with pytest.raises(ValueError):
        cload_pairix.callback(
            op.join(testdir, "data", "hg19.bins.2000kb.bed.gz"),
            op.join(testdir, "data", "hg19.GM12878-MboI.pairs.subsample.blksrt.txt.gz"),
            testcool_path,
            metadata=None,
            assembly="hg19",
            nproc=1,
            zero_based=False,
            max_split=2,
            block_char="?",
        )


@pytest.mark.skipif(
    _pandas_major_version < 1, reason="hash fix only works with pandas >= 1.0"
)
@pytest.mark.parametrize(
    "bins_path,pairs_path,ref_path",
    [
        (
            op.join(testdir, "data", "toy.bins.var.bed"),
            op.join(testdir, "data", "toy.pairs"),
            op.join(testdir, "data", "toy.symm.upper.var.cool"),
        ),
        (
            op.join(testdir, "data", "toy.bins.var.bed"),
            op.join(testdir, "data", "toy_hash.pairs"),
            op.join(testdir, "data", "toy.symm.upper.var.cool"),
        ),
        (
            op.join(testdir, "data", "toy.bins.var.bed"),
            op.join(testdir, "data", "toy_hash.pairs.gz"),
            op.join(testdir, "data", "toy.symm.upper.var.cool"),
        ),
    ],
)
def test_cload_pairs(bins_path, pairs_path, ref_path):
    kwargs = dict(
        metadata=None,
        assembly="hg19",
        zero_based=False,
        comment_char="#",
        input_copy_status="unique",
        no_symmetric_upper=False,
        field=(),
        chunksize=int(15e6),
        mergebuf=int(20e6),
        max_merge=200,
        temp_dir=None,
        no_delete_temp=False,
        storage_options=None,
        no_count=False,
        chrom1=2,
        pos1=3,
        chrom2=4,
        pos2=5,
        append=False,
    )
    cload_pairs.callback(bins_path, pairs_path, testcool_path, **kwargs)
    with h5py.File(testcool_path, "r") as f1, h5py.File(ref_path, "r") as f2:
        assert np.all(f1["pixels/bin1_id"][:] == f2["pixels/bin1_id"][:])
        assert np.all(f1["pixels/bin2_id"][:] == f2["pixels/bin2_id"][:])
        assert np.all(f1["pixels/count"][:] == f2["pixels/count"][:])
    try:
        os.remove(testcool_path)
    except OSError:
        pass


@pytest.mark.parametrize(
    "bins_path,pairs_path",
    [
        (
            op.join(testdir, "data", "toy.chrom.sizes") + ":2",
            op.join(testdir, "data", "toy.pairs"),
        )
    ],
)
def test_cload_field(bins_path, pairs_path):
    kwargs = dict(
        metadata=None,
        assembly="toy",
        zero_based=False,
        comment_char="#",
        input_copy_status="unique",
        no_symmetric_upper=False,
        chunksize=10,
        mergebuf=10,
        max_merge=200,
        temp_dir=None,
        no_delete_temp=False,
        storage_options=None,
        no_count=True,
        chrom1=2,
        pos1=3,
        chrom2=4,
        pos2=5,
        append=False,
    )
    cload_pairs.callback(
        bins_path, pairs_path, testcool_path, field=("score=8:dtype=float",), **kwargs
    )
    pixels = cooler.Cooler(testcool_path).pixels()[:]
    assert "count" in pixels.columns and types.is_integer_dtype(pixels.dtypes["count"])
    assert "score" in pixels.columns and types.is_float_dtype(pixels.dtypes["score"])


@pytest.mark.parametrize(
    "bins_path,pairs_path",
    [
        (
            op.join(testdir, "data", "toy.chrom.sizes") + ":2",
            op.join(testdir, "data", "toy.pairs"),
        )
    ],
)
def test_cload_custom_tempdir(bins_path, pairs_path):
    for temp_dir in [op.join(testdir, "data"), "-"]:
        cload_pairs.callback(
            bins_path,
            pairs_path,
            testcool_path,
            metadata=None,
            assembly="toy",
            zero_based=False,
            comment_char="#",
            input_copy_status="unique",
            no_symmetric_upper=False,
            field=(),
            chunksize=10,
            mergebuf=10,
            temp_dir=temp_dir,
            no_delete_temp=False,
            max_merge=200,
            storage_options=None,
            no_count=True,
            chrom1=2,
            pos1=3,
            chrom2=4,
            pos2=5,
            append=False,
        )
        pixels = cooler.Cooler(testcool_path).pixels()[:]
        assert "count" in pixels.columns and types.is_integer_dtype(
            pixels.dtypes["count"]
        )


def test_load_bg2_vs_coo():
    kwargs = dict(
        metadata=None,
        assembly="hg19",
        chunksize=int(20e6),
        mergebuf=int(15e6),
        max_merge=200,
        field=(),
        count_as_float=False,
        one_based=False,
        comment_char="#",
        input_copy_status="unique",
        no_symmetric_upper=False,
        temp_dir=None,
        no_delete_temp=False,
        storage_options=None,
        append=False,
    )

    out_path1 = op.join(tmp, "test1.cool")
    out_path2 = op.join(tmp, "test2.cool")

    load.callback(
        op.join(testdir, "data", "hg19.bins.2000kb.bed.gz"),
        op.join(testdir, "data", "hg19.GM12878-MboI.matrix.2000kb.bg2.gz"),
        out_path1,
        format="bg2",
        **kwargs,
    )
    load.callback(
        op.join(testdir, "data", "hg19.bins.2000kb.bed.gz"),
        op.join(testdir, "data", "hg19.GM12878-MboI.matrix.2000kb.coo.txt"),
        out_path2,
        format="coo",
        **kwargs,
    )

    with h5py.File(out_path1, "r") as f1, h5py.File(out_path2, "r") as f2:
        for col in ["bin1_id", "bin2_id", "count"]:
            assert np.all(f1["pixels"][col][:] == f2["pixels"][col][:])

    for fp in [out_path1, out_path2]:
        try:
            os.remove(fp)
        except OSError:
            pass


def test_load_zero_one_based_bg2():
    kwargs = dict(
        format="bg2",
        metadata=None,
        assembly="toy",
        chunksize=10,
        mergebuf=10,
        max_merge=200,
        field=(),
        count_as_float=False,
        comment_char="#",
        input_copy_status="unique",
        no_symmetric_upper=False,
        temp_dir=None,
        no_delete_temp=False,
        storage_options=None,
        append=False,
    )
    # 1-based-start BG2 input
    ref = "toy.symm.upper.1.ob.bg2"
    bins_path = op.join(testdir, "data", "toy.chrom.sizes") + ":1"
    pixels_path = op.join(testdir, "data", ref)
    load.callback(bins_path, pixels_path, testcool_path, one_based=True, **kwargs)
    # reference, 1-based starts
    ref_df = pd.read_csv(
        pixels_path,
        sep="\t",
        names=["chrom1", "start1", "end1", "chrom2", "start2", "end2", "count"],
    )
    # output
    out_df = cooler.Cooler(testcool_path).pixels(join=True)[:]
    out_df["start1"] += 1
    out_df["start2"] += 1
    assert np.all(out_df == ref_df)

    # 0-based-start BG2 input
    ref = "toy.symm.upper.1.zb.bg2"
    bins_path = op.join(testdir, "data", "toy.chrom.sizes") + ":1"
    pixels_path = op.join(testdir, "data", ref)
    load.callback(bins_path, pixels_path, testcool_path, one_based=False, **kwargs)
    # reference, 0-based starts
    ref_df = pd.read_csv(
        pixels_path,
        sep="\t",
        names=["chrom1", "start1", "end1", "chrom2", "start2", "end2", "count"],
    )
    # output
    out_df = cooler.Cooler(testcool_path).pixels(join=True)[:]
    assert np.all(out_df == ref_df)


def test_load_zero_one_based_coo():
    kwargs = dict(
        format="coo",
        metadata=None,
        assembly="toy",
        chunksize=10,
        mergebuf=10,
        max_merge=200,
        field=(),
        count_as_float=False,
        comment_char="#",
        input_copy_status="unique",
        no_symmetric_upper=False,
        temp_dir=None,
        no_delete_temp=False,
        storage_options=None,
        append=False,
    )
    # 1-based-start COO input
    ref = "toy.symm.upper.1.ob.coo"
    bins_path = op.join(testdir, "data", "toy.chrom.sizes") + ":1"
    pixels_path = op.join(testdir, "data", ref)
    load.callback(bins_path, pixels_path, testcool_path, one_based=True, **kwargs)
    # reference, 1-based starts
    ref_df = pd.read_csv(pixels_path, sep="\t", names=["bin1_id", "bin2_id", "count"])
    # output
    out_df = cooler.Cooler(testcool_path).pixels()[:]
    out_df["bin1_id"] += 1
    out_df["bin2_id"] += 1
    assert np.all(out_df == ref_df)

    # 0-based-start COO input
    ref = "toy.symm.upper.1.zb.coo"
    bins_path = op.join(testdir, "data", "toy.chrom.sizes") + ":1"
    pixels_path = op.join(testdir, "data", ref)
    load.callback(bins_path, pixels_path, testcool_path, one_based=False, **kwargs)
    # reference, 0-based starts
    ref_df = pd.read_csv(pixels_path, sep="\t", names=["bin1_id", "bin2_id", "count"])
    # output
    out_df = cooler.Cooler(testcool_path).pixels()[:]
    assert np.all(out_df == ref_df)


def test_array_loader():
    chromsizes = cooler.util.read_chromsizes(
        op.join(testdir, "data", "toy.chrom.sizes")
    )
    bins = cooler.util.binnify(chromsizes, 10)
    n = len(bins)
    array = np.ones((n, n))
    iterator = cooler.create.ArrayLoader(bins, array, chunksize=100)
    list(iterator)

    array = np.ones((n + 1, n + 1))
    with pytest.raises(ValueError):
        cooler.create.ArrayLoader(bins, array, chunksize=100)


def test_cload_append_mode():
    for append in (True, False):
        out_path = op.join(tmp, "test.cool")
        with h5py.File(out_path, "w") as f:
            f.attrs["xxxx"] = True
        kwargs = dict(
            metadata=None,
            assembly="hg19",
            chunksize=int(15e6),
            mergebuf=int(15e6),
            zero_based=False,
            comment_char="#",
            input_copy_status="unique",
            no_symmetric_upper=False,
            field=(),
            temp_dir=None,
            no_delete_temp=False,
            storage_options=None,
            no_count=False,
            max_merge=200,
            chrom1=2,
            pos1=3,
            chrom2=4,
            pos2=5,
            append=append,
        )
        cload_pairs.callback(
            op.join(testdir, "data", "toy.bins.var.bed"),
            op.join(testdir, "data", "toy.pairs"),
            out_path,
            **kwargs,
        )
        with h5py.File(out_path, "r") as f:
            if append:
                assert "xxxx" in f.attrs
            else:
                assert "xxxx" not in f.attrs
        try:
            os.remove(out_path)
        except OSError:
            pass


def test_load_append_mode():
    for append in (True, False):
        out_path = op.join(tmp, "test.cool")
        with h5py.File(out_path, "w") as f:
            f.attrs["xxxx"] = True
        kwargs = dict(
            metadata=None,
            assembly="hg19",
            chunksize=int(20e6),
            mergebuf=int(20e6),
            max_merge=200,
            field=(),
            count_as_float=False,
            one_based=False,
            comment_char="#",
            input_copy_status="unique",
            no_symmetric_upper=False,
            temp_dir=None,
            no_delete_temp=False,
            storage_options=None,
            append=append,
        )
        load.callback(
            op.join(testdir, "data", "hg19.bins.2000kb.bed.gz"),
            op.join(testdir, "data", "hg19.GM12878-MboI.matrix.2000kb.coo.txt"),
            out_path,
            format="coo",
            **kwargs,
        )
        with h5py.File(out_path, "r") as f:
            if append:
                assert "xxxx" in f.attrs
            else:
                assert "xxxx" not in f.attrs
        try:
            os.remove(out_path)
        except OSError:
            pass


def test_load_custom_tempdir():
    for temp_dir in [op.join(testdir, "data"), "-"]:
        load.callback(
            op.join(testdir, "data", "hg19.bins.2000kb.bed.gz"),
            op.join(testdir, "data", "hg19.GM12878-MboI.matrix.2000kb.coo.txt"),
            testcool_path,
            format="coo",
            metadata=None,
            assembly="hg19",
            chunksize=int(20e6),
            mergebuf=int(20e6),
            max_merge=200,
            field=(),
            count_as_float=False,
            one_based=False,
            comment_char="#",
            input_copy_status="unique",
            no_symmetric_upper=False,
            temp_dir=temp_dir,
            no_delete_temp=False,
            storage_options=None,
            append=False,
        )
        pixels = cooler.Cooler(testcool_path).pixels()[:]
        assert "count" in pixels.columns and types.is_integer_dtype(
            pixels.dtypes["count"]
        )
