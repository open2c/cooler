# import multiprocess as mp
# import click

# from ..util import parse_cooler_uri


# def parse_bins_arg(arg):
#     # Provided chromsizes and binsize
#     if ":" in arg:
#         chromsizes_file, binsize = arg.split(":")
#         if not op.exists(chromsizes_file):
#             raise ValueError('File "{}" not found'.format(chromsizes_file))
#         try:
#             binsize = int(binsize)
#         except ValueError:
#             raise ValueError(
#                 'Expected integer binsize argument (bp), got "{}"'.format(binsize))
#         chromsizes = util.read_chromsizes(chromsizes_file, all_names=True)
#         bins = util.binnify(chromsizes, binsize)

#     # Provided bins
#     elif op.exists(arg):
#         try:
#             bins = pd.read_csv(
#                 arg,
#                 sep='\t',
#                 names=['chrom', 'start', 'end'],
#                 usecols=[0, 1, 2],
#                 dtype={'chrom': str})
#         except pd.parser.CParserError as e:
#             raise ValueError(
#                 'Failed to parse bins file "{}": {}'.format(arg, str(e)))

#         chromtable = (
#             bins.drop_duplicates(['chrom'], keep='last')[['chrom', 'end']]
#                 .reset_index(drop=True)
#                 .rename(columns={'chrom': 'name', 'end': 'length'})
#         )
#         chroms, lengths = list(chromtable['name']), list(chromtable['length'])
#         chromsizes = pd.Series(index=chroms, data=lengths)
        
#     else:
#         raise ValueError(
#             'Expected BINS to be either <Path to bins file> or '
#             '<Path to chromsizes file>:<binsize in bp>.')

#     return chromsizes, bins


# def check_ncpus(arg_value):
#     arg_value = int(arg_value)

#     if arg_value <= 0:
#         raise argparse.ArgumentTypeError("n_cpus must be >= 1")
#     else:
#         return min(arg_value, mp.cpu_count())


# class CoolerURI(click.Path):
# 	def convert(self, value, param, ctx):
# 		file_path, group_path = parse_cooler_uri(value)
# 		file_path = super(CoolerURI, self).convert(file_path, param, ctx)
# 		return file_path + '::' + group_path
