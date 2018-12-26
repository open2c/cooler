from contextlib import contextmanager
from functools import wraps
import multiprocess as mp
import os.path as op
import errno
import sys
import os

import pandas as pd
import numpy as np
import click
import six

from .. import util


class DelimitedTuple(click.types.ParamType):
    def __init__(self, sep=',', type=str):
        self.sep = sep
        self.type = click.types.convert_type(type)

    @property
    def name(self):
        return "separated[%s]" % self.sep

    def convert(self, value, param, ctx):
        # needs to pass through value = None unchanged
        # needs to be idempotent
        # needs to be able to deal with param and context being None
        if value is None:
            return value
        elif isinstance(value, six.string_types):
            parts = value.split(',')
        else:
            parts = value
        return tuple(self.type(x, param, ctx) for x in parts)


def parse_kv_list_param(arg, item_sep=',', kv_sep='='):
    import yaml
    from io import StringIO

    if item_sep != ',':
        arg = arg.replace(item_sep, ',')
    arg = '{' + arg.replace(kv_sep, ': ') + '}'
    try:
        result = yaml.load(StringIO(arg))
    except yaml.YAMLError as e:
        raise click.BadParameter(
            "Error parsing key-value pairs: {}".format(arg))
    return result


def parse_field_param(arg, includes_colnum=True, includes_agg=True):
    parts = arg.split(':')
    prefix = parts[0]
    if len(parts) == 1:
        props = None
    elif len(parts) == 2:
        props = parts[1]
    else:
        raise click.BadParameter(arg)

    if includes_colnum:
        parts = prefix.split('=')
        name = parts[0]
        if len(parts) == 1:
            colnum = None
        elif len(parts) == 2:
            try:
                colnum = int(parts[1]) - 1
            except ValueError:
                raise click.BadParameter(
                    "Not a number: '{}'".format(parts[1]), param_hint=arg)
            if colnum < 0:
                raise click.BadParameter(
                    "Field numbers start at 1.", param_hint=arg)
        else:
            raise click.BadParameter(arg)
    else:
        name = parts[0]
        colnum = None

    dtype = None
    agg = None
    if props is not None:
        for item in props.split(','):
            try:
                prop, value = item.split('=')
            except ValueError:
                raise click.BadParameter(arg)
            if prop == 'dtype':
                dtype = np.dtype(value)
            elif prop == 'agg' and includes_agg:
                agg = value
            else:
                raise click.BadParameter(
                    "Invalid property: '{}'.".format(prop),
                    param_hint=arg)
    return name, colnum, dtype, agg


def parse_bins(arg):
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


def check_ncpus(arg_value):
    arg_value = int(arg_value)

    if arg_value <= 0:
        raise click.BadParameter("n_cpus must be >= 1")
    else:
        return min(arg_value, mp.cpu_count())


@contextmanager
def on_broken_pipe(handler):
    try:
        yield
    except IOError as e:
        if e.errno == errno.EPIPE:
            handler(e)
        else:
            # Not a broken pipe error. Bubble up.
            raise


def exit_on_broken_pipe(exit_code):
    """
    Decorator to catch a broken pipe (EPIPE) error and exit cleanly.

    Use this decorator to prevent the "[Errno 32] Broken pipe" output message.

    Notes
    -----
    A SIGPIPE signal is sent to a process writing to a pipe while the other
    end has been closed. For example, this happens when piping output to
    programs like head(1). Python traps this signal and translates it into an
    exception. It is presented as an IOError in PY2, and a subclass of OSError
    in PY3 (aliased by IOError), both using POSIX error number 32 (EPIPE).

    Some programs exit with 128 + signal.SIGPIPE == 141. However, according to
    the example in the docs, Python exits with the generic error code of 1
    on EPIPE.

    The equivalent system error when trying to write on a socket which has
    been shutdown for writing is ESHUTDOWN (108). It also raises
    BrokenPipeError on PY3.

    [1] https://docs.python.org/3.7/library/signal.html#note-on-sigpipe
    [2] https://www.quora.com/How-can-you-avoid-a-broken-pipe-error-on-Python

    """
    def decorator(func):
        @wraps(func)
        def decorated(*args, **kwargs):
            try:
                func(*args, **kwargs)
            except IOError as e:
                if e.errno == errno.EPIPE:
                    # We caught a broken pipe error.
                    # Python flushes standard streams on exit; redirect remaining
                    # output to devnull to avoid another BrokenPipeError at shutdown.
                    devnull = os.open(os.devnull, os.O_WRONLY)
                    os.dup2(devnull, sys.stdout.fileno())
                    sys.exit(exit_code)
                else:
                    # Not a broken pipe error. Bubble up.
                    raise
        return decorated
    return decorator
