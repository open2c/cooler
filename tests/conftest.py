from scipy import sparse
import numpy as np

import pytest
import mock


class MockGroup(dict):
    file = mock.Mock(['mode'])
    attrs = mock.MagicMock(dict)

    def __getitem__(self, path):
        current_item = self
        if set(path) == {'/'}:
            return self
        for item in path.split('/'):
            current_item = dict.__getitem__(current_item, item)
        return current_item

    def __setitem__(self, path, value):
        if isinstance(value, dict):
            value = MockGroup()
        # else:
        #     value = mock.Mock(h5py.Dataset)
        current_item = self
        previous_item = None
        for item in path.split('/'):
            previous_item = current_item
            try:
                current_item = dict.__getitem__(current_item, item)
            except KeyError:
                dict.__setitem__(current_item, item, MockGroup())
                current_item = dict.__getitem__(current_item, item)
        dict.__setitem__(previous_item, item, value)

    def __delitem__(self, path):
        current_item = self
        previous_item = None
        for item in path.split('/'):
            previous_item = current_item
            current_item = dict.__getitem__(current_item, item)
        dict.__delitem__(previous_item, item)


class MockCooler(MockGroup):

    @classmethod
    def make_random(cls, chrom_offsets, binsize, density):
        chrom_nbins = np.diff(chrom_offsets)
        assert chrom_offsets[0] == 0 and np.all(np.diff(chrom_offsets) >= 0)
        n_chroms = len(chrom_offsets) - 1
        n_bins = chrom_offsets[-1]
        chroms = {
            'name': np.array(
                ['chr' + str(i) for i in range(1, n_chroms + 1)],
                dtype='S'),
            'length': np.array(
                [chrom_nbins[i] * binsize for i in range(n_chroms)],
                dtype=np.int32)
        }
        bins = {
            'chrom': np.concatenate(
                [[i] * chrom_nbins[i]
                    for i in range(n_chroms)]),
            'start': np.concatenate([
                np.arange(0, chrom_nbins[i] * binsize, binsize)
                    for i in range(n_chroms)]),
            'end': np.concatenate([
                np.arange(binsize, chrom_nbins[i] * (binsize + 1), binsize)
                    for i in range(n_chroms)])
        }
        r = sparse.random(n_bins, n_bins, density=density, random_state=1)
        r = sparse.triu(r, k=1).tocsr()
        pixels = {
            'bin1_id':  r.tocoo().row,
            'bin2_id':  r.indices,
            'count':    r.data,
        }
        indexes = {
            'chrom_offset': np.array(chrom_offsets, dtype=np.int32),
            'bin1_offset': r.indptr,
        }
        return cls(chroms, bins, pixels, indexes, binsize)

    def __init__(self, chroms, bins, pixels, indexes, binsize):
        self.file = self
        self.file.mode = 'r'
        self.file.filename = 'mock.cool'
        self.name = '/'
        self.attrs = {
            'bin-size': binsize,
            'bin-type': 'fixed',
            'symmetric-storage-mode': 'upper',
            'nchroms': len(chroms),
            'nbins': len(bins['chrom']),
            'nnz': len(pixels['bin1_id']),
            'metadata': '{}',
        }
        super(MockCooler, self).__init__({
            'chroms':  chroms,
            'bins': bins,
            'pixels': pixels,
            'indexes': indexes,
        })



@pytest.fixture(params=[([0, 10, 20], 100, 1)])
def mock_cooler(request):
    chrom_offsets, binsize, density = request.param
    return MockCooler.make_random(chrom_offsets, binsize, density)
