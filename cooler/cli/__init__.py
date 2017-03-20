# -*- coding: utf-8 -*-
from __future__ import division, print_function
import logging
import sys

import click
from .. import __version__, get_logger

logging.basicConfig(stream=sys.stderr)
logger = get_logger()
logger.setLevel(logging.INFO)


# Monkey patch
click.core._verify_python3_env = lambda: None


CONTEXT_SETTINGS = {
    'help_option_names': ['-h', '--help'],
}


@click.version_option(version=__version__)
@click.group(context_settings=CONTEXT_SETTINGS)
@click.option(
    '--debug/--no-debug', 
    help="Verbose logging", 
    default=False)
def cli(debug):
    if debug:
        logger.setLevel(logging.DEBUG)


from . import (
    makebins,
    digest,
    csort,
    cload,
    load,
    balance,
    dump,
    show,
    info,
    coarsegrain,
)
