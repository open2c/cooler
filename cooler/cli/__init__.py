# -*- coding: utf-8 -*-
from __future__ import division, print_function
import click
from .. import __version__


# Monkey patch
click.core._verify_python3_env = lambda: None


CONTEXT_SETTINGS = {
    'help_option_names': ['-h', '--help'],
}


@click.version_option(version=__version__)
@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass


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
