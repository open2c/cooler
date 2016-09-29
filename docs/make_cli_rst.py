# -*- coding: utf-8 -*-
"""
Generate reST markup for CLI documentation by gathering command help texts.

"""
from cooler.cli import cli
import click


COMMANDS = [
    'makebins',
    'digest',
    'csort',
    'cload',
    'balance',
    'info',
    'dump',
    'show'
]


TEMPLATE = """\
.. _cli-reference:

CLI Reference
=============

.. toctree::
   :maxdepth: 1


cooler
------

::

{cooler}


"""

for cmd in COMMANDS:
    TEMPLATE += """\
cooler {0}
----------------

::

{{{0}}}


""".format(cmd)


def indent(s, width):
    return '\n'.join((' ' * width) + line for line in s.split('\n'))


def helptext(subcommand=None):
    if subcommand is None:
        ctx = click.Context(
            cli,
            info_name='cooler')
    else:
        ctx = click.Context(
            cli.commands.get(subcommand),
            info_name='cooler ' + subcommand)

    return indent(ctx.get_help(), 4)


text = TEMPLATE.format(
    cooler=helptext(),
    **{cmd: helptext(cmd) for cmd in COMMANDS}
)


with open('cli.rst', 'w') as f:
    f.write(text)
