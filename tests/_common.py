import contextlib
import tempfile
import shutil
import os

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
        assert np.allclose(f1['pixels/bin1_id'][:], f2['pixels/bin1_id'][:])
        assert np.allclose(f1['pixels/bin2_id'][:], f2['pixels/bin2_id'][:])
        assert np.allclose(f1['pixels/count'][:], f2['pixels/count'][:])
