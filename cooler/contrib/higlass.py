from __future__ import division, print_function, absolute_import
 
import numpy as np
import pandas as pd
import cooler
import h5py
import logging
 
logger = logging.getLogger(__name__)
 
TILE_SIZE = 256


def abs_coord_2_bin(c, abs_pos, chroms, chrom_cum_lengths, chrom_sizes):
    """Get bin ID from absolute coordinates.
 
    Args:
        c (Cooler): Cooler instance of a .cool file.
        abs_pos (int): Absolute coordinate to be translated.
 
    Returns:
        int: Bin number.
    """
 
    try:
        chr_id = np.flatnonzero(chrom_cum_lengths > abs_pos)[0] - 1
    except IndexError:
        return c.info['nbins']
 
    chrom = chroms[chr_id]
    rel_pos = abs_pos - chrom_cum_lengths[chr_id]
 
    return c.offset((chrom, rel_pos, chrom_sizes[chrom]))


def get_chromosome_names_cumul_lengths(c):
    '''
    Get the chromosome names and cumulative lengths:
 
    Args:
 
    c (Cooler): A cooler file
 
    Return:
 
    (names, sizes, lengths) -> (list(string), dict, np.array(int))
    '''
    chrom_names = c.chromnames
    chrom_sizes = dict(c.chromsizes)
    chrom_cum_lengths = np.r_[0, np.cumsum(c.chromsizes.values)]
    return chrom_names, chrom_sizes, chrom_cum_lengths
 
 
def get_data(f, start_pos_1, end_pos_1, start_pos_2, end_pos_2, transform='default'):
    """Get balanced pixel data.
 
    Args:
        f: h5py.File
            An HDF5 Group that contains the cooler for this resolution
        start_pos_1 (int): Test.
        end_pos_1 (int): Test.
        start_pos_2 (int): Test.
        end_pos_2 (int): Test.
 
    Returns:
        DataFrame: Annotated cooler pixels.
    """
 
    c = cooler.Cooler(f)
 
    (chroms, chrom_sizes, chrom_cum_lengths) = get_chromosome_names_cumul_lengths(c)
 
    i0 = abs_coord_2_bin(c, start_pos_1, chroms, chrom_cum_lengths, chrom_sizes)
    i1 = abs_coord_2_bin(c, end_pos_1, chroms, chrom_cum_lengths, chrom_sizes)

    j0 = abs_coord_2_bin(c, start_pos_2, chroms, chrom_cum_lengths, chrom_sizes)
    j1 = abs_coord_2_bin(c, end_pos_2, chroms, chrom_cum_lengths, chrom_sizes)

    '''
    print('i', i0, i1)
    print('j', j0, j1)
    '''
    matrix = c.matrix(as_pixels=True, balance=False, max_chunk=np.inf)

    if i0 >= matrix.shape[0] or j0 >= matrix.shape[1]:
        # query beyond the bounds of the matrix
        # return an empty matrix
        i0,i1,i1,j1 = 0,0,0,0
    else:
        # limit the range of the query to be within bounds
        i1 = min(i1, matrix.shape[0]-1)
        j1 = min(j1, matrix.shape[1]-1)

    #print("size", matrix.shape)

    pixels = matrix[i0:i1+1, j0:j1+1]

    if not len(pixels):
        return pd.DataFrame(columns=['genome_start1', 'genome_start2', 'balanced'])
 
    # select bin columns to extract
    cols = ['chrom', 'start', 'end']
    if (transform == 'default' and 'weight' in c.bins()) or transform == 'weight':
        cols.append('weight')
    elif transform in ('KR', 'VC', 'VC_SQRT'):
        cols.append(transform)

    bins = c.bins(convert_enum=False)[cols]    
    pixels = cooler.annotate(pixels, bins)
    pixels['genome_start1'] = chrom_cum_lengths[pixels['chrom1']] + pixels['start1']
    pixels['genome_start2'] = chrom_cum_lengths[pixels['chrom2']] + pixels['start2']

    # apply transform
    if (transform == 'default' and 'weight' in c.bins()) or transform == 'weight':
        pixels['balanced'] = (
            pixels['count'] * pixels['weight1'] * pixels['weight2']
        )
        return pixels[['genome_start1', 'genome_start2', 'balanced']]
    elif transform in ('KR', 'VC', 'VC_SQRT'):
        pixels['balanced'] = (
            pixels['count'] / pixels[transform+'1'] / pixels[transform+'2']
        )
        return pixels[['genome_start1', 'genome_start2', 'balanced']]
    else:
        return pixels[['genome_start1', 'genome_start2', 'count']]

 
def get_info(file_path):
    """Get information of a cooler file.
 
    Args:
        file_path (str): Path to a cooler file.
 
    Returns:
        dict: Dictionary containing basic information about the cooler file.
    """
 
    with h5py.File(file_path, 'r') as f:
        max_zoom = f.attrs.get('max-zoom')
 
        if max_zoom is None:
            logger.info('no zoom found')
            raise ValueError(
                'The `max_zoom` attribute is missing.'
            )
 
        c = cooler.Cooler(f["0"])
 
        (chroms, chrom_sizes, chrom_cum_lengths) = get_chromosome_names_cumul_lengths(c)
 
        total_length = int(chrom_cum_lengths[-1])
        max_zoom = f.attrs['max-zoom']
        bin_size = int(f[str(max_zoom)].attrs['bin-size'])
 
        max_width = bin_size * TILE_SIZE * 2**max_zoom
        
        # the list of available data transforms
        transforms = {}
        
        for i in range(max_zoom):
            f_for_zoom = f[str(i)]['bins']

            if 'weight' in f_for_zoom:
                transforms['weight'] = {'name': 'ICE', 'value': 'weight'}
            if 'KR' in f_for_zoom:
                transforms['KR'] = {'name': 'KR', 'value': 'KR'}
            if 'VC' in f_for_zoom:
                transforms['VC'] = {'name': 'VC', 'value': 'VC'}
            if 'VC_SQRT' in f_for_zoom:
                transforms['VC_SQRT'] = {'name': 'VC_SQRT', 'value': 'VC_SQRT'}

        info = {
            'min_pos': [0.0, 0.0],
            'max_pos': [total_length, total_length],
            'max_zoom': max_zoom,
            'max_width': max_width,
            'bins_per_dimension': TILE_SIZE,
            'transforms': transforms.values()
        }
 
    return info


def get_quadtree_depth(chromsizes, binsize):
    """
    Depth of quad tree necessary to tesselate the concatenated genome with quad
    tiles such that linear dimension of the tiles is a preset multiple of the 
    genomic resolution.

    """
    tile_size_bp = TILE_SIZE * binsize
    min_tile_cover = np.ceil(sum(chromsizes) / tile_size_bp)
    return int(np.ceil(np.log2(min_tile_cover)))


def get_zoom_resolutions(chromsizes, base_res):
    return [base_res * 2**x for x in range(get_quadtree_depth(chromsizes, base_res)+1)]


def print_zoom_resolutions(chromsizes_file, base_res):
    """
    Print comma-separated list of zoom resolutions for a given genome
    and base resolution.

    """
    chromsizes = cooler.util.read_chromsizes(chromsizes_file, all_names=True)
    resolutions = get_zoom_resolutions(chromsizes, base_res)
    print(','.join(str(res) for res in resolutions))
