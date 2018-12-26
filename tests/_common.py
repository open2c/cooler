import contextlib
import tempfile
import shutil
import os

from pandas.api.types import is_numeric_dtype
import numpy as np
import h5py
import cooler


@contextlib.contextmanager
def isolated_filesystem():
    """A context manager that creates a temporary folder and changes
    the current working directory to it for isolated filesystem tests.
    """
    cwd = os.getcwd()
    t = tempfile.mkdtemp()
    os.chdir(t)
    try:
        yield t
    finally:
        os.chdir(cwd)
        try:
            shutil.rmtree(t)
        except (OSError, IOError):
            pass


def cooler_cmp(uri1, uri2):
    c1 = cooler.Cooler(uri1)
    c2 = cooler.Cooler(uri2)
    with c1.open('r') as f1, \
         c2.open('r') as f2:
         for path in (
                'chroms/name', 'chroms/length',
                'bins/chrom', 'bins/start', 'bins/end',
                'pixels/bin1_id', 'pixels/bin2_id', 'pixels/count'):
            dset1, dset2 = f1[path], f2[path]
            dtype = dset1.dtype
            assert dtype == dset2.dtype
            if is_numeric_dtype(dtype):
                assert np.allclose(dset1[:], dset2[:])
            else:
                assert np.all(dset1[:] == dset2[:])
