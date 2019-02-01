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

            dtype1 = dset1.dtype
            dtype2 = dset2.dtype

            if dtype1.kind == 'S':
                # Null padding of ascii arrays is not guaranteed to be
                # preserved so we only check the kind.
                assert dtype2.kind == 'S'
            else:
                assert dtype1 == dtype2

            if is_numeric_dtype(dtype1):
                assert np.allclose(dset1[:], dset2[:])
            else:
                assert np.all(dset1[:] == dset2[:])
