# -*- coding: utf-8 -*-
from __future__ import division, print_function
import logging
import sys
from .. import __version__, get_logger
import click


# Monkey patch
click.core._verify_python3_env = lambda: None


class UnsortedGroup(click.Group):
    def list_commands(self, ctx):
        return list(self.commands)

    # def format_commands(self, ctx, formatter):
    #     """Extra format methods for multi methods that adds all the commands
    #     after the options.
    #     """
    #     commands = []
    #     for subcommand in self.list_commands(ctx):
    #         cmd = self.get_command(ctx, subcommand)
    #         # What is this, the tool lied about a command.  Ignore it
    #         if cmd is None:
    #             continue
    #         if cmd.hidden:
    #             continue

    #         commands.append((subcommand, cmd))

    #     # allow for 3 times the default spacing
    #     if len(commands):
    #         limit = formatter.width - 6 - max(len(cmd[0]) for cmd in commands)

    #         rows = []
    #         for subcommand, cmd in commands:
    #             help = cmd.get_short_help_str(limit)
    #             rows.append((subcommand, help))

    #         if rows:
    #             with formatter.section('Commands'):
    #                 formatter.write_dl(rows)



CONTEXT_SETTINGS = {
    'help_option_names': ['-h', '--help'],
}


@click.version_option(version=__version__)
@click.group(context_settings=CONTEXT_SETTINGS, cls=UnsortedGroup)
@click.option(
    '--debug/--no-debug',
    help="Verbose logging",
    default=False)
@click.option(
    '-pm', '--post-mortem',
    help="Post mortem debugging",
    is_flag=True,
    default=False)
def cli(debug, post_mortem):
    """
    Type -h or --help after any subcommand for more information.

    """
    logging.basicConfig(stream=sys.stderr)
    logging.captureWarnings(True)
    root_logger = get_logger()
    if debug:
        root_logger.setLevel(logging.DEBUG)
    else:
        root_logger.setLevel(logging.INFO)

    if post_mortem:
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
    info,
    dump,
    show,
    balance,
    merge,
    coarsen,
    zoomify,
    makebins,
    digest,
    csort,
    fileops
)
