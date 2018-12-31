# -*- coding: utf-8 -*-
from __future__ import division, print_function
import logging
import sys
import os
from .._version import __version__
from .._logging import get_logger
import click

# Monkey patch
click.core._verify_python3_env = lambda: None


CONTEXT_SETTINGS = {
    'help_option_names': ['-h', '--help'],
}

class UnsortedGroup(click.Group):
    def list_commands(self, ctx):
        return list(self.commands)


@click.version_option(__version__, '-V', '--version')
@click.group(context_settings=CONTEXT_SETTINGS, cls=UnsortedGroup)
@click.option(
    '-v', '--verbose',
    help="Verbose logging.",
    count=True)
@click.option(
    '-d', '--debug',
    help="On error, drop into the post-mortem debugger shell.",
    is_flag=True,
    default=False)
def cli(verbose, debug):
    """
    Type -h or --help after any subcommand for more information.

    """
    # Initialize logging to stderr
    logging.basicConfig(stream=sys.stderr)
    logging.captureWarnings(True)
    root_logger = get_logger()

    # Set verbosity level
    if verbose > 0:
        root_logger.setLevel(logging.DEBUG)
        if verbose > 1:
            try:
                import psutil
                import atexit

                @atexit.register
                def process_dump_at_exit():
                    process_attrs = [
                        'cmdline',
                        # 'connections',
                        'cpu_affinity',
                        'cpu_num',
                        'cpu_percent',
                        'cpu_times',
                        'create_time',
                        'cwd',
                        # 'environ',
                        'exe',
                        # 'gids',
                        'io_counters',
                        'ionice',
                        'memory_full_info',
                        # 'memory_info',
                        # 'memory_maps',
                        'memory_percent',
                        'name',
                        'nice',
                        'num_ctx_switches',
                        'num_fds',
                        'num_threads',
                        'open_files',
                        'pid',
                        'ppid',
                        'status',
                        'terminal',
                        'threads',
                        # 'uids',
                        'username',
                    ]
                    p = psutil.Process()
                    info = p.as_dict(process_attrs, ad_value='')
                    for key in process_attrs:
                        root_logger.debug(
                            "PSINFO:'{}': {}".format(key, info[key]))

            except ImportError:
                root_logger.warning(
                    "Install psutil to see process information.")

    else:
        root_logger.setLevel(logging.INFO)

    # Set hook for postmortem debugging
    if debug:
        import traceback
        try:
            import ipdb as pdb
        except ImportError:
            import pdb
        def _excepthook(exc_type, value, tb):
            traceback.print_exception(exc_type, value, tb)
            print()
            pdb.pm()
        sys.excepthook = _excepthook


from . import (
    cload,
    load,
    merge,
    coarsen,
    zoomify,
    balance,
    info,
    dump,
    show,
    makebins,
    digest,
    csort,
    fileops
)
