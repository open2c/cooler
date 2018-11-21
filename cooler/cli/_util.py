# import multiprocess as mp
import os.path as op
import pandas as pd
import numpy as np
from .. import util
import click


def _parse_field_params(args):
    extra_fields = []
    bad_param = False
    for arg in args:

        parts = arg.split(',')
        if len(parts) == 1 or len(parts) > 3:
            bad_param = True
        else:
            name = parts[0]
            try:
                number = int(parts[1]) - 1
            except ValueError:
                bad_param = True

            if number < 0:
                raise click.BadParameter(
                    "Field numbers are assumed to be 1-based.")

            if len(parts) == 3:
                dtype = np.dtype(parts[2])
            else:
                dtype = None

        if bad_param:
            raise click.BadParameter(
                "Expected '--field {{name}},{{number}}' "
                "or '--field {{name}},{{number}},{{dtype}}'; "
                "got '{}'".format(arg))
        extra_fields.append((name, number, dtype))

    return extra_fields


def _parse_bins(arg):
    # Provided chromsizes and binsize
    if ":" in arg:
        chromsizes_file, binsize = arg.split(":")
        if not op.exists(chromsizes_file):
            raise ValueError('File "{}" not found'.format(chromsizes_file))
        try:
            binsize = int(binsize)
        except ValueError:
            raise ValueError(
                'Expected integer binsize argument (bp), got "{}"'.format(binsize))
        chromsizes = util.read_chromsizes(chromsizes_file, all_names=True)
        bins = util.binnify(chromsizes, binsize)

    # Provided bins
    elif op.exists(arg):
        try:
            bins = pd.read_csv(
                arg,
                sep='\t',
                names=['chrom', 'start', 'end'],
                usecols=[0, 1, 2],
                dtype={'chrom': str})
        except pd.parser.CParserError as e:
            raise ValueError(
                'Failed to parse bins file "{}": {}'.format(arg, str(e)))

        chromtable = (
            bins.drop_duplicates(['chrom'], keep='last')[['chrom', 'end']]
                .reset_index(drop=True)
                .rename(columns={'chrom': 'name', 'end': 'length'})
        )
        chroms, lengths = list(chromtable['name']), list(chromtable['length'])
        chromsizes = pd.Series(index=chroms, data=lengths)
        
    else:
        raise ValueError(
            'Expected BINS to be either <Path to bins file> or '
            '<Path to chromsizes file>:<binsize in bp>.')

    return chromsizes, bins


# def check_ncpus(arg_value):
#     arg_value = int(arg_value)

#     if arg_value <= 0:
#         raise argparse.ArgumentTypeError("n_cpus must be >= 1")
#     else:
#         return min(arg_value, mp.cpu_count())
