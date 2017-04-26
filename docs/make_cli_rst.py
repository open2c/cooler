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
    'load',
    'list',
    'info',
    'copy',
    'dump',
    'show',
    'balance',
    'merge',
    'coarsen',
    'zoomify',
]

SUBCOMMANDS = {
    'cload': ['hiclib', 'pairix', 'tabix']
}


TEMPLATE = """\
.. _cli-reference:

CLI Reference
=============

.. toctree::
   :maxdepth: 1


Quick reference
---------------

+------------------------+
| Helper commands        |
+========================+
| `cooler makebins`_     |
+------------------------+
| `cooler digest`_       |
+------------------------+
| `cooler csort`_        |
+------------------------+

+------------------------+
| Ingesting data         |
+========================+
| `cooler cload`_        |
+------------------------+
| `cooler load`_         |
+------------------------+

+------------------------+
| File manipulation/info |
+========================+
| `cooler list`_         |
+------------------------+
| `cooler info`_         |
+------------------------+
| `cooler copy`_         |
+------------------------+

+------------------------+
| Export/visualization   |
+========================+
| `cooler dump`_         |
+------------------------+
| `cooler show`_         |
+------------------------+

+------------------------+
| Operations             |
+========================+
| `cooler balance`_      |
+------------------------+
| `cooler merge`_        |
+------------------------+
| `cooler coarsen`_      |
+------------------------+
| `cooler zoomify`_      |
+------------------------+


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


def helptext(command=None):
    if command is None:
        ctx = click.Context(
            cli,
            info_name='cooler')
        text = indent(ctx.get_help(), 4)
    else:
        ctx = click.Context(
            cli.commands.get(command),
            info_name='cooler ' + command)
        text = indent(ctx.get_help(), 4)

        if command in SUBCOMMANDS:
            for subcommand in SUBCOMMANDS[command]:
                info_name = 'cooler {} {}'.format(
                        command, subcommand)
                ctx = click.Context(
                    cli.commands
                       .get(command)
                       .commands
                       .get(subcommand),
                    info_name=info_name)
                subtext = '\n\n' + info_name + '\n'
                subtext += '~' * len(info_name) + '\n'
                subtext += ctx.get_help() + '\n'
                text += indent(subtext, 8) 

    return text


text = TEMPLATE.format(
    cooler=helptext(),
    **{cmd: helptext(cmd) for cmd in COMMANDS}
)


with open('cli.rst', 'w') as f:
    f.write(text)
