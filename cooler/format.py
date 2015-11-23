# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals
from datetime import datetime
import warnings
import six

from scipy.sparse import coo_matrix
import numpy as np
import pandas
import h5py

from .genome import Region


def from_dense(h5, chrom_table, bin_table, heatmap, metadata=None, 
               bintype='fixed', h5opts=None):
    metadata = {} if metadata is None else metadata
    h5opts = {} if h5opts is None else h5opts

    n_chroms = len(chrom_table)
    n_bins = len(bin_table)
    if len(heatmap) != n_bins:
        raise ValueError(
            "length mismatch:" + 
            " heatmap length is {0}, bin table length is {1}".format(
                len(heatmap), n_bins))

    bintype = bintype.lower()
    if bintype not in ('fixed', 'variable'):
        raise ValueError("Bin type must be one of 'fixed' or 'variable'")

    # Sparsify matrix
    i, j = np.nonzero(heatmap)
    mask = i >= j
    tril_i = i[mask]
    tril_j = j[mask]
    values = heatmap[tril_i, tril_j]
    nnz = len(values)

    # Attributes
    print('metadata')
    h5.attrs['id'] = metadata.get('id', "No ID")

    h5.attrs['bin-type'] = bintype
    if bintype == 'fixed':
        binsize = bin_table['end'].iat[0] - bin_table['start'].iat[0]
    else:
        binsize = None
    h5.attrs['bin-size'] = binsize
    h5.attrs['genome-assembly'] = metadata.get('genome-assembly', 'unknown')
    h5.attrs['format-url'] = "https://bitbucket.org/nvictus/cooler"
    h5.attrs['format-version'] = (0, 1)
    h5.attrs['generated-by'] = metadata.get('generated-by', "cooler")
    h5.attrs['creation-date'] = datetime.now().isoformat()
    h5.attrs['shape'] = heatmap.shape
    h5.attrs['nnz'] = nnz

    # Data groups
    print('sequence assemblies')
    grp = h5.create_group('scaffolds')
    grp.create_dataset('name', 
                        shape=(len(chrom_table),), dtype='S32',
                        data=np.array(chrom_table['name'], dtype='S32'),
                        **h5opts)
    grp.create_dataset('length', 
                        shape=(len(chrom_table),), dtype=np.int64,
                        data=chrom_table['length'],
                        **h5opts)

    print('bins')
    grp = h5.create_group('bins')
    if 'id' not in chrom_table.columns:
        chrom_table['id'] = np.arange(len(chrom_table))
    chrom_ids = chrom_table['id'].loc[bin_table['chrom']]
    grp.create_dataset('chrom_id',
                        shape=(len(bin_table),), dtype=np.int32,
                        data=chrom_ids,
                        **h5opts)
    grp.create_dataset('start',
                        shape=(len(bin_table),), dtype=np.int64,
                        data=bin_table['start'], 
                        **h5opts)
    grp.create_dataset('end',
                        shape=(len(bin_table),), dtype=np.int64,
                        data=bin_table['end'], 
                        **h5opts)

    print('matrix')
    grp = h5.create_group('matrix')
    grp.create_dataset('bin1_id',
                        shape=(len(values),), dtype=np.int32,
                        data=tril_i, **h5opts)
    grp.create_dataset('bin2_id',
                        shape=(len(values),), dtype=np.int32,
                        data=tril_j, **h5opts)
    grp.create_dataset('count',
                        shape=(len(values),), dtype=np.int32,
                        data=values, **h5opts)

    # Indexes
    print('indexes')
    chrom_binedges = np.r_[0, np.cumsum(np.ceil(chrom_table['length']/binsize))]
    mat_lo = np.searchsorted(tril_i, np.arange(n_bins), side='left')
    mat_hi = np.r_[mat_lo[1:], nnz]

    grp = h5.create_group('indexes')
    idx1 = grp.create_group('chrom_to_bin')
    idx1.create_dataset('bin_lo',
                         shape=(n_chroms,), dtype=np.int64,
                         data=chrom_binedges[:-1], **h5opts)
    idx1.create_dataset('bin_hi',
                         shape=(n_chroms,), dtype=np.int64,
                         data=chrom_binedges[1:], **h5opts)
    idx2 = grp.create_group('bin_to_matrix')
    idx2.create_dataset('mat_lo',
                         shape=(n_bins,), dtype=np.int32,
                         data=mat_lo, **h5opts)
    idx2.create_dataset('mat_hi',
                         shape=(n_bins,), dtype=np.int32,
                         data=mat_hi, **h5opts)


def rebuild_indexes(h5):
    pass


def _get_region_extent(h5, chrom_id, region, binsize=None):
    chrom, start, end = region
    if binsize is not None:
        bin_chrom_lo = h5['indexes']['chrom_to_bin']['bin_lo'][chrom_id]
        bin_lo = bin_chrom_lo + int(np.floor(start/binsize))
        bin_hi = bin_chrom_lo + int(np.ceil(end/binsize))
    else:
        bin_chrom_lo = h5['indexes']['chrom_to_bin']['bin_lo'][chrom_id]
        bin_chrom_hi = h5['indexes']['chrom_to_bin']['bin_hi'][chrom_id]
        bin_chrom = h5['bins']['start'][bin_chrom_lo:bin_chrom_hi]
        bin_lo = bin_chrom_lo + np.searchsorted(bin_chrom, start, 'left')
        bin_hi = bin_chrom_lo + np.searchsorted(bin_chrom, end, 'right')
    return bin_lo, bin_hi


def _get_slice_from_bins(h5, bin1_lo, bin1_hi, bin2_lo, bin2_hi):    
    bin1_ids = np.arange(bin1_lo, bin1_hi)
    hm1_lo = h5['indexes']['bin_to_matrix']['mat_lo'][bin1_lo:bin1_hi]
    hm1_hi = h5['indexes']['bin_to_matrix']['mat_hi'][bin1_lo:bin1_hi]
    if bin1_hi - bin1_lo != 0 or bin2_hi - bin2_lo != 0:
        i, j, v = [], [], []
        for bin1_id, lo1, hi1 in zip(bin1_ids, hm1_lo, hm1_hi):
            hm_bin2 = h5['matrix']['bin2_id'][lo1:hi1]
            lo = lo1 + np.searchsorted(hm_bin2, bin2_lo)
            hi = lo1 + np.searchsorted(hm_bin2, bin2_hi)
            i.append(np.zeros(hi-lo) + bin1_id)
            j.append(h5['matrix']['bin2_id'][lo:hi])
            v.append(h5['matrix']['count'][lo:hi])
        i = np.concatenate(i, axis=0).astype(int)
        j = np.concatenate(j, axis=0).astype(int)
        v = np.concatenate(v, axis=0)
        i -= bin1_lo
        j -= bin2_lo
    else:
        i, j, v = np.array([], dtype=int), np.array([], dtype=int), np.array([])
    m = bin1_hi - bin1_lo
    n = bin2_hi - bin2_lo
    return (i, j, v), (m, n)


def _get_slice_overlap(h5, region1, region2):
    # Two corner cases: the regions partially overlap or one is inside the other
    transpose = False
    if region1.start < region2.start:
        region2, region1 = region1, region2
        transpose = True

    chrom = region1.chrom
    if region2.comes_before(region1):    
        b1a1 = Region((chrom, region2.start, region1.start))
        a1b2 = Region((chrom, region1.start, region2.end))
        b2a2 = Region((chrom, region2.end, region1.end))
        (i1, j1, v1), (m1, n1) = get_slice(h5, region1, b1a1)
        (i2, j2, v2), (m2, n2) = get_slice(h5, a1b2)
        (i3, j3, v3), (m3, n3) = get_slice(h5, b2a2, a1b2)
        j2 += n1
        i3 += m2
        j3 += n1
        i, j, v = np.r_[i1, i2, i3], np.r_[j1, j2, j3], np.r_[v1, v2, v3]
        m, n = m1, n1 + n2
    elif region2.contains(region1):
        b1a1 = Region((chrom, region2.start, region1.start))
        a1a2 = Region((chrom, region1.start, region1.end))
        a2b2 = Region((chrom, region1.end, region2.end))
        (i1, j1, v1), (m1, n1) = get_slice(h5, a1a2, b1a1)
        (i2, j2, v2), (m2, n2) = get_slice(h5, a1a2)
        (j3, i3, v3), (n3, m3) = get_slice(h5, a2b2, a1a2)
        j2 += n1
        j3 += n1 + n2
        i, j, v = np.r_[i1, i2, i3], np.r_[j1, j2, j3], np.r_[v1, v2, v3]
        m, n = m1, n1 + n2 + n3
    else:
        raise RuntimeError("This shouldn't happen")

    if transpose:
        i, j = j, i
        m, n = n, m
    return (i, j, v), (m, n)


def get_slice(h5, region=None, region2=None):
    # If no regions specified, return everything
    if region is None:
        i = h5['matrix']['bin1_id'][:]
        j = h5['matrix']['bin2_id'][:]
        v = h5['matrix']['count'][:]
        m, n = h5.attrs['shape']
        return (i, j, v), (m, n)

    # Parse the region queries
    if region2 is None: region2 = region
    chrom_table = get_scaffolds(h5)
    binsize = h5.attrs.get('bin-size', None)
    region1 = Region(region, chrom_table['length'])
    region2 = Region(region2, chrom_table['length'])
    chrom1_id = chrom_table['id'].at[region1[0]]
    chrom2_id = chrom_table['id'].at[region2[0]]
    
    # Three query cases: same, different but overlapping, non-overlapping
    if region1 == region2:
        bin1_lo, bin1_hi = _get_region_extent(h5, chrom1_id, region1, binsize)
        bin2_lo, bin2_hi = bin1_lo, bin1_hi
        (i, j, v), (m, n) = _get_slice_from_bins(
            h5, bin1_lo, bin1_hi, bin2_lo, bin2_hi)
        i, j, v = np.r_[i, j], np.r_[j, i], np.r_[v, v]
        return (i, j, v), (m, n)
    elif region1.overlaps(region2):
        return _get_slice_overlap(h5, region1, region2)
    else:
        # Since the data is lower triangular, the region along axis1 region must 
        # come strictly before the region along axis0, otherwise swap and 
        # transpose the result.
        transpose = False
        if chrom1_id < chrom2_id or region1.start < region2.start:
            chrom2_id, chrom1_id = chrom1_id, chrom2_id
            region2, region1 = region1, region2
            transpose = True
        bin1_lo, bin1_hi = _get_region_extent(h5, chrom1_id, region1, binsize)
        bin2_lo, bin2_hi = _get_region_extent(h5, chrom2_id, region2, binsize)
        (i, j, v), (m, n) = _get_slice_from_bins(
            h5, bin1_lo, bin1_hi, bin2_lo, bin2_hi)
        if transpose:
            i, j, m, n = j, i, n, m
        return (i, j, v), (m, n)


def get_matrix(h5, region=None, region2=None, dense=False):
    (i, j, v), (m, n) = get_slice(h5, region, region2)
    A = coo_matrix((v, (i, j)), (m, n))
    if dense:
        A = A.toarray()
    return A


def get_scaffolds(h5):
    names = h5['scaffolds']['name'][:].astype('U')
    lengths = h5['scaffolds']['length'][:]
    return pandas.DataFrame(
        {'name': names, 'id': np.arange(len(names)), 'length': lengths}, 
        columns=['name', 'id', 'length'], index=names)


def get_bins(h5):
    names = h5['scaffolds']['name'][:].astype('U')
    chrom_ids = h5['bins']['chrom_id'][:]
    chroms = names[chrom_ids]
    starts = h5['bins']['start'][:]
    ends = h5['bins']['end'][:]
    return pandas.DataFrame(
        {'chrom': chroms, 'start': starts, 'end': ends},
        columns=['chrom', 'start', 'end'])

def get_metadata(h5):
    return dict(h5.attrs.items())
