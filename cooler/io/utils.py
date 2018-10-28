# -*- coding: utf-8 -*-
from __future__ import division, print_function
import posixpath
import tempfile
import warnings
from six.moves import map
import six

from pandas.api.types import is_categorical, is_integer
import pandas as pd
import numpy as np
import h5py

from ..core import put
from ..util import get_binsize, get_chromsizes, infer_meta, get_meta, parse_region
from .. import get_logger


logger = get_logger()


def parse_cooler_uri(s):
    """
    Parse a Cooler URI string

    e.g. /path/to/mycoolers.cool::/path/to/cooler

    """
    parts = s.split('::')
    if len(parts) == 1:
        file_path, group_path = parts[0], '/'
    elif len(parts) == 2:
        file_path, group_path = parts
        if not group_path.startswith('/'):
            group_path = '/' + group_path
    else:
        raise ValueError("Invalid Cooler URI string")
    return file_path, group_path


def is_cooler(filepath, group=None):
    """
    Determine if a file contains a valid Cooler data hierarchy.

    Parameters
    ----------
    filepath : str
    group : str, optional
        Path to specific group to check. Otherwise returns True
        if any Cooler paths are detected.

    Returns
    -------
    bool
    
    """
    if not h5py.is_hdf5(filepath):
        return False
    if group is None:
        return len(ls(filepath)) > 0
    if not group.startswith('/'):
        group = '/' + group
    return group in ls(filepath)


def is_multires_cooler():
    pass


def check_bins(bins, chromsizes):
    is_cat = pd.api.types.is_categorical(bins['chrom'])
    bins = bins.copy()
    if not is_cat:
        bins['chrom'] = pd.Categorical(
            bins.chrom, 
            categories=list(chromsizes.index), 
            ordered=True)
    else:
        assert (bins['chrom'].cat.categories == chromsizes.index).all()

    return bins


def balanced_partition(gs, n_chunk_max, file_contigs, loadings=None):
    n_bins = len(gs.bins)
    grouped = gs._bins_grouped
    
    chrom_nbins = grouped.size()
    if loadings is None:
        loadings = chrom_nbins
    chrmax = loadings.idxmax()
    loadings = loadings / loadings.loc[chrmax]
    const = chrom_nbins.loc[chrmax] / n_chunk_max

    granges = []
    for chrom, group in grouped:
        if chrom not in file_contigs:
            continue
        clen = gs.chromsizes[chrom]
        step = int(np.ceil(const / loadings.loc[chrom]))        
        anchors = group.start.values[::step]
        if anchors[-1] != clen:
            anchors = np.r_[anchors, clen]
        granges.extend( (chrom, start, end) 
            for start, end in zip(anchors[:-1], anchors[1:]))
    return granges


class GenomeSegmentation(object):
    def __init__(self, chromsizes, bins):
        bins = check_bins(bins, chromsizes)
        self._bins_grouped = bins.groupby('chrom', sort=False)
        nbins_per_chrom = self._bins_grouped.size().values

        self.chromsizes = chromsizes
        self.binsize = get_binsize(bins)
        self.contigs = list(chromsizes.keys())
        self.bins = bins
        self.idmap = pd.Series(
            index=chromsizes.keys(), 
            data=range(len(chromsizes)))
        self.chrom_binoffset = np.r_[0, np.cumsum(nbins_per_chrom)]
        self.chrom_abspos = np.r_[0, np.cumsum(chromsizes.values)]
        self.start_abspos = (self.chrom_abspos[bins['chrom'].cat.codes] + 
                             bins['start'].values)
    
    def fetch(self, region):
        chrom, start, end = parse_region(region, self.chromsizes)
        result = self._bins_grouped.get_group(chrom)
        if start > 0 or end < self.chromsizes[chrom]:
            lo = result['end'].values.searchsorted(start, side='right')
            hi = lo + result['start'].values[lo:].searchsorted(end, side='left')
            result = result.iloc[lo:hi]
        return result


def rename_scaffolds(clr, scaffolds, h5opts=None):
    
    from .creation import CHROM_DTYPE, CHROMID_DTYPE

    if h5opts is None:
        h5opts = dict(compression='gzip', compression_opts=6)
    
    with clr.open('r+') as h5:
        chroms = cooler.core.get(h5['chroms']).set_index('name')
        n_chroms = len(chroms)
        new_names = np.array(chroms.rename(scaffolds).index.values, 
                             dtype=CHROM_DTYPE)  # auto-adjusts char length
        
        del h5['chroms/name']
        h5['chroms'].create_dataset('name',
                           shape=(n_chroms,),
                           dtype=new_names.dtype,
                           data=new_names,
                           **h5opts)
        
        bins = cooler.core.get(h5['bins'])
        n_bins = len(bins)
        idmap = dict(zip(new_names, range(n_chroms)))
        if is_categorical(bins['chrom']) or is_integer(bins['chrom']):            
            chrom_ids = bins['chrom'].cat.codes
            chrom_dtype = h5py.special_dtype(enum=(CHROMID_DTYPE, idmap))
            del h5['bins/chrom']
            try:
                chrom_dset = h5['bins'].create_dataset('chrom',
                                   shape=(n_bins,),
                                   dtype=chrom_dtype,
                                   data=chrom_ids,
                                   **h5opts)
            except ValueError:
                # If HDF5 enum header would be too large, 
                # try storing chrom IDs as raw int instead
                chrom_dtype = CHROMID_DTYPE
                chrom_dset = h5['bins'].create_dataset('chrom',
                                   shape=(n_bins,),
                                   dtype=chrom_dtype,
                                   data=chrom_ids,
                                   **h5opts)


def ls(filepath):
    """
    Traverse a file's data hierarchy and list all Cooler nodes.

    Parameters
    ----------
    filepath : str

    Returns
    -------
    list of Cooler group paths in the file
    
    """
    from .creation import MAGIC, URL
    listing = []
    keys = ['chroms', 'bins', 'pixels', 'indexes']
    def _check_group(pth, grp):
        fmt = grp.attrs.get('format', None)
        url = grp.attrs.get('format-url', None)
        if fmt == MAGIC or url == URL:
            if not all(name in grp.keys() for name in keys):
                warnings.warn(
                    'Cooler path /{} appears to be corrupt'.format(pth))
            listing.append('/' + pth if not pth.startswith('/') else pth)

    with h5py.File(filepath, 'r') as f:
        _check_group('/', f)
        f.visititems(_check_group)

    return listing


def _copy(src_uri, dst_uri, overwrite, link, rename, soft_link):
    """
    Copy a Cooler from one file to another or within the same file.

    See also: h5copy, h5repack tools from HDF5 suite

    \b\bArguments:

    SRC_URI : Path to source file or URI to source Cooler group

    DST_URI : Path to destination file or URI to destination Cooler group

    """
    src_path, src_group = parse_cooler_uri(src_uri)
    dst_path, dst_group = parse_cooler_uri(dst_uri)

    if sum([link, rename, soft_link]) > 1:
        raise ValueError(
            'Must provide at most one of: "link", "rename", "soft_link"')

    if not os.path.isfile(dst_path) or overwrite:
        write_mode = 'w'
    else:
        write_mode = 'r+'

    with h5py.File(src_path, 'r+') as src, \
         h5py.File(dst_path, write_mode) as dst:

        # if dst_group in dst and dst_group != '/':
        #     click.confirm(
        #         "A group named '{}' already exists in '{}'. Overwrite?".format(
        #             dst_group, dst_path), 
        #         abort=True)
        #     del dst[dst_group]

        if src_path == dst_path:
            if link or rename:
                src[dst_group] = src[src_group]
                if rename:
                    del src[src_group]
            elif soft_link:
                src[dst_group] = h5py.SoftLink(src_group)
        else:
            if link:
                raise OSError("Can't hard link between two different files.")
            elif soft_link:
                dst[dst_group] = h5py.ExternalLink(src_path, src_group)
            else:
                if dst_group == '/':
                    for subgrp in src[src_group].keys():
                        src.copy(src_group + '/' + subgrp, dst, subgrp)
                    dst[dst_group].attrs.update(src[src_group].attrs)
                else:
                    src.copy(
                        src_group, dst, 
                        dst_group if dst_group != '/' else None)


def cp(src_uri, dst_uri, overwrite=False):
    _copy(src_uri, dst_uri, overwrite, link=False, rename=False, soft_link=False)


def mv(src_uri, dst_uri, overwrite=False):
    _copy(src_uri, dst_uri, overwrite, link=False, rename=True, soft_link=False)


def ln(src_uri, dst_uri, soft=False, overwrite=False):
    _copy(src_uri, dst_uri, overwrite, link=not soft, rename=False, soft_link=soft)
