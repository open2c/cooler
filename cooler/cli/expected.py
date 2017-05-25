# # -*- coding: utf-8 -*-
# from __future__ import division, print_function
# from multiprocess import Pool
# from functools import partial
# from itertools import chain
# import itertools
# import argparse
# import sys

# from scipy.signal import fftconvolve
# import numpy as np
# import pandas

# import cooler
# from cooler.tools import split, partition

# where = np.flatnonzero
# concat = chain.from_iterable


# def bedslice(grouped, chrom, start, end):
#     """Assumes no proper nesting of intervals"""
#     chromdf = grouped.get_group(chrom)
#     lo = chromdf['end'].values.searchsorted(start, side='right')
#     hi = lo + chromdf['start'].values[lo:].searchsorted(end, side='left')
#     return chromdf.iloc[lo:hi]


# def apply_by_chrom(df, func, *args, chrom_col='chrom', **kwargs):
#     grouped = df.groupby(chrom_col, sort=False)
#     results = []
#     for name, group in grouped:
#         result = func(name, group, *args, **kwargs)
#         if np.isscalar(result):
#             result = pandas.DataFrame([result] * len(group), index=group.index)
#         elif not isinstance(result, pandas.DataFrame):
#             result = pandas.DataFrame(result, index=group.index)
#         results.append(result)
#     return pandas.concat(results, axis=0)


# def geomprog(factor, start=1):
#     yield start
#     while True:
#         start *= factor
#         yield start


# def logsmooth(x, factor):
#     y = np.zeros_like(x)
#     steps = []
#     pos = 0
#     for i, step in enumerate(geomprog(factor)):
#         if step > len(x) - i:
#             break 
#         istep = min(len(x) - i, int(step))
#         y[i] = x[i:i+istep].mean()
#     return y


# def pdist1dcount_fft(n, locs):
#     """
#     See <http://stackoverflow.com/questions/42423823/distribution-of-pairwise-distances-between-many-integers>.
    
#     """
#     x = np.zeros(n)
#     x[locs] = 1
#     return np.round(fftconvolve(x, x[::-1], mode='full')).astype(int)[-n:]


# def count_all_pixels_per_diag(n):
#     return np.arange(n, 0, -1)


# def count_bad_pixels_per_diag(n, bad_bins):
#     """
#     Efficiently count the number of bad pixels on each positive diagonal of a 
#     matrix assuming a sequence of bad bins forms a "grid" of invalid pixels.
    
#     Each bad bin bifurcates into two a row and column of bad pixels, so the 
#     maximum count is 2*k, where k is the number of bad bins. For a given 
#     diagonal, one needs to subtract the contribution from "out-of-bounds" 
#     regions and the contribution of the intersection points of bad rows with 
#     bad columns.
    
#     Parameters
#     ----------
#     n : int
#         total number of bins
#     bad_bins : 1D array of int
#         sorted array of bad bin indexes
#     Returns
#     -------
#     dcount : 1D array of length n
#         dcount[d] == number of bad bins on diagonal d
    
#     """
#     k = len(bad_bins)
        
#     # Store all intersection pixels in a separate array
#     # ~O(n log n) with fft
#     ixn = pdist1dcount_fft(n, bad_bins)
    
#     # Keep track of out-of-bounds pixels
#     # ~O(n + 2k), I think
#     dcount = np.zeros(n, dtype=int)
#     dcount[0] = ixn[0]
#     pl = 0
#     pr = k
#     for d in range(1, n):
#         if pl < k:
#             bl = bad_bins[pl] - d
#             while bl < 0:
#                 pl += 1
#                 if pl == k:
#                     break
#                 bl = bad_bins[pl] - d
#         if pr > 0:
#             br = bad_bins[pr-1] + d
#             while br >= n:
#                 pr -= 1
#                 if pr == 0:
#                     break
#                 br = bad_bins[pr-1] + d
#         dcount[d] = 2*k - ixn[d] - pl - (k - pr)
#     return dcount


# def _accum_by_cisdiag(c, bins, span):
#     """Sum properties along the diagonals of the intrachromosomal matrices"""
#     lo, hi = span
#     pixels = c.pixels()[lo:hi]
    
#     # assign chroms and filter for cis records
#     pixels = cooler.annotate(pixels, bins[['chrom', 'weight']], replace=False)
#     pixels = pixels[pixels.chrom1 == pixels.chrom2].copy()
    
#     # assign diagonal indices
#     pixels = pixels.rename(columns={'chrom1': 'chrom'})
#     pixels['diag'] = pixels['bin2_id'] - pixels['bin1_id']

#     # balance
#     pixels['balanced'] = pixels['count'] * pixels['weight1'] * pixels['weight2']
#     pixels['balanced2'] = pixels['balanced'] * pixels['balanced']

#     # group by diagonal and accumulate
#     grouped = pixels.groupby(['chrom', 'diag'], sort=False)
#     agg = grouped.aggregate({
#             'balanced': np.sum,
#             'balanced2': np.sum,
#         })
#     return agg.reset_index()


# def compute_expected(c, binsize, drop_diags, chunksize, map_impl=map, 
#                      regions=None, smooth_factor=None):
#     bins = c.bins()[:]
    
#     if regions is None:
#         #names = [item[0] for item in bins.groupby('chrom', sort=False)]
#         groups = [item[1] for item in bins.groupby('chrom', sort=False)]
#     else:
#         groups = []
#         g = bins.groupby('chrom', sort=False)
#         for _, region in regions.iterrows():
#             #names.extend([region['name']] * len(g))
#             groups.append(
#                 bedslice(g, region['chrom'], region['start'], region['end']))
    
#     n_bins_per_group = [len(g) for g in groups]
#     bad_bins_per_group = [where(np.isnan(g['weight'].values)) for g in groups]
    
#     # initialize
#     ex = bins[['chrom']].copy()
#     #ex['name'] = bins
#     ex['diag'] = bins['start'] // binsize
#     ex['balanced'] = 0
#     ex['bad']   = list(concat(map_impl(count_bad_pixels_per_diag, 
#                                        n_bins_per_group, 
#                                        bad_bins_per_group)))
#     ex['total'] = list(concat(map_impl(count_all_pixels_per_diag, 
#                                        n_bins_per_group)))

#     # split records into chunks
#     args = partition(0, len(c.pixels()), chunksize)

#     # apply + combine
#     combined = pandas.concat(
#         map_impl(partial(_accum_by_cisdiag, c, bins), args),
#         axis=0,
#         ignore_index=True)
#     combined = combined.groupby(['chrom', 'diag']).sum()
    
#     ex = ex.set_index(['chrom', 'diag'])
#     ex = ex.add(combined, fill_value=0)
#     ex = ex.reset_index()

#     if smooth_factor is not None:
#         ex['balanced'] = apply_by_chrom(
#             ex, lambda c,g: logsmooth(g['balanced'], smooth_factor))
#         ex['balanced2'] = apply_by_chrom(
#             ex, lambda c,g: logsmooth(g['balanced2'], smooth_factor))

#     # average over valid elements
#     n = ex['total'] - ex['bad']
#     ex['average'] = ex['balanced'] / n
#     ex['std'] = np.sqrt(
#         ex['balanced2'] / n - (ex['balanced'] / n)**2
#     )
    
#     # mask out bad diagonals
#     ex.loc[ex['diag'] < drop_diags, 'average'] = np.nan
#     ex.loc[ex['diag'] < drop_diags, 'std'] = np.nan

#     return ex


# def variadic_poolmap(pool):
#     def decorated(job, *args, **kwargs):
#         def job_wrapper(a):
#             return job(*a)
#         return pool.map(job_wrapper, zip(*args), **kwargs)
#     return decorated


# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(
#         description="Compute expected")
#     parser.add_argument(
#         "cool_path",
#         help="Cooler file",
#         metavar="COOL_PATH")
#     parser.add_argument(
#         "--drop-diags",
#         type=int,
#         default=2)
#     parser.add_argument(
#         "--chunksize", "-c",
#         type=int,
#         default=int(10e6))
#     parser.add_argument(
#         "--nproc", "-p",
#         type=int,
#         default=8)
#     parser.add_argument(
#         "--out", "-o",
#         type=str)

#     args = vars(parser.parse_args())

#     c = cooler.Cooler(args['cool_path'])
#     binsize = c.info['bin-size']
#     drop_diags = args['drop_diags']
    
#     with Pool(args['nproc']) as pool:
#         ex = compute_expected(
#             c, 
#             binsize, 
#             args['drop_diags'], 
#             args['chunksize'], 
#             map_impl=variadic_poolmap(pool))

#     try:
#         out = args['out']
#         if out is None:
#             f = sys.stdout
#         else:
#             f = open(out, 'wt')
#         ex[['chrom', 'diag', 'average', 'std']].to_csv(
#             f, sep='\t', index=False, header=True)
#     except (IOError, OSError) as e:
#         if e.errno == 32:  # broken pipe
#             try:
#                 f.close()
#             except OSError:
#                 pass
#         else:
#             raise
#     else:
#         f.close()
