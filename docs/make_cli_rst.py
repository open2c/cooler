"""
Code adapted from the sphinx-click project.

<https://github.com/click-contrib/sphinx-click>

"""
import click
from docutils import statemachine

import cooler.cli


def _indent(text, level=1):
    prefix = " " * (4 * level)

    def prefixed_lines():
        for line in text.splitlines(True):
            yield (prefix + line if line.strip() else line)

    return "".join(prefixed_lines())


def _get_usage(ctx):
    """Alternative, non-prefixed version of 'get_usage'."""
    formatter = ctx.make_formatter()
    pieces = ctx.command.collect_usage_pieces(ctx)
    formatter.write_usage(ctx.command_path, " ".join(pieces), prefix="")
    return formatter.getvalue().rstrip("\n")


def _get_help_record(opt):
    """Re-implementation of click.Opt.get_help_record.
    The variant of 'get_help_record' found in Click makes uses of slashes to
    separate multiple opts, and formats option arguments using upper case. This
    is not compatible with Sphinx's 'option' directive, which expects
    comma-separated opts and option arguments surrounded by angle brackets [1].
    [1] http://www.sphinx-doc.org/en/stable/domains.html#directive-option
    """

    def _write_opts(opts):
        rv, _ = click.formatting.join_options(opts)
        if not opt.is_flag and not opt.count:
            rv += f" <{opt.name}>"
        return rv

    rv = [_write_opts(opt.opts)]
    if opt.secondary_opts:
        rv.append(_write_opts(opt.secondary_opts))

    help = opt.help or ""
    extra = []
    if opt.default is not None and opt.show_default:
        extra.append(
            "default: {}".format(
                ", ".join("%s" % d for d in opt.default)
                if isinstance(opt.default, (list, tuple))
                else opt.default
            )
        )
    if opt.required:
        extra.append("required")
    if extra:
        help = "{}[{}]".format(help and help + "  " or "", "; ".join(extra))

    return ", ".join(rv), help


def _format_description(ctx):
    """Format the description for a given `click.Command`.
    We parse this as reStructuredText, allowing users to embed rich
    information in their help messages if they so choose.
    """
    help_string = ctx.command.help or ctx.command.short_help
    if not help_string:
        return

    bar_enabled = False
    for line in statemachine.string2lines(
        help_string, tab_width=4, convert_whitespace=True
    ):
        if line == "\b":
            bar_enabled = True
            continue
        if line == "":
            bar_enabled = False
        line = "| " + line if bar_enabled else line
        yield line
    yield ""


def _format_usage(ctx):
    """Format the usage for a `click.Command`."""
    yield ".. code-block:: shell"
    yield ""
    for line in _get_usage(ctx).splitlines():
        yield _indent(line)
    yield ""


def _format_option(opt):
    """Format the output for a `click.Option`."""
    opt = _get_help_record(opt)

    yield f".. option:: {opt[0]}"
    if opt[1]:
        yield ""
        for line in statemachine.string2lines(
            opt[1], tab_width=4, convert_whitespace=True
        ):
            yield _indent(line)


def _format_options(ctx):
    """Format all `click.Option` for a `click.Command`."""
    # the hidden attribute is part of click 7.x only hence use of getattr
    params = [
        x
        for x in ctx.command.params
        if isinstance(x, click.Option) and not getattr(x, "hidden", False)
    ]

    for param in params:
        yield from _format_option(param)
        yield ""


def _format_argument(arg):
    """Format the output of a `click.Argument`."""
    yield f".. option:: {arg.human_readable_name}"
    yield ""
    yield _indent(
        "{} argument{}".format(
            "Required" if arg.required else "Optional", "(s)" if arg.nargs != 1 else ""
        )
    )


def _format_arguments(ctx):
    """Format all `click.Argument` for a `click.Command`."""
    params = [x for x in ctx.command.params if isinstance(x, click.Argument)]

    for param in params:
        yield from _format_argument(param)
        yield ""


def _format_envvar(param):
    """Format the envvars of a `click.Option` or `click.Argument`."""
    yield f".. envvar:: {param.envvar}"
    yield "   :noindex:"
    yield ""
    if isinstance(param, click.Argument):
        param_ref = param.human_readable_name
    else:
        # if a user has defined an opt with multiple "aliases", always use the
        # first. For example, if '--foo' or '-f' are possible, use '--foo'.
        param_ref = param.opts[0]

    yield _indent(f"Provide a default for :option:`{param_ref}`")


def _format_envvars(ctx):
    """Format all envvars for a `click.Command`."""
    params = [x for x in ctx.command.params if getattr(x, "envvar")]

    for param in params:
        yield ".. _{command_name}-{param_name}-{envvar}:".format(
            command_name=ctx.command_path.replace(" ", "-"),
            param_name=param.name,
            envvar=param.envvar,
        )
        yield ""
        yield from _format_envvar(param)
        yield ""


def _format_subcommand(command):
    """Format a sub-command of a `click.Command` or `click.Group`."""
    yield f".. object:: {command.name}"

    if command.short_help:
        yield ""
        for line in statemachine.string2lines(
            command.short_help, tab_width=4, convert_whitespace=True
        ):
            yield _indent(line)


def _get_lazyload_commands(multicommand):
    commands = {}
    for command in multicommand.list_commands(multicommand):
        commands[command] = multicommand.get_command(multicommand, command)

    return commands


def _filter_commands(ctx, commands=None):
    """Return list of used commands."""
    lookup = getattr(ctx.command, "commands", {})
    if not lookup and isinstance(ctx.command, click.MultiCommand):
        lookup = _get_lazyload_commands(ctx.command)

    if commands is None:
        return sorted(lookup.values(), key=lambda item: item.name)

    names = [name.strip() for name in commands.split(",")]
    return [lookup[name] for name in names if name in lookup]


def format_command(ctx, show_nested, commands=None):
    """Format the output of `click.Command`."""
    # description

    yield from _format_description(ctx)

    yield f".. program:: {ctx.command_path}"

    # usage

    yield from _format_usage(ctx)

    # arguments

    lines = list(_format_arguments(ctx))
    if lines:
        yield ".. rubric:: Arguments"
        yield ""

    yield from lines

    # options

    lines = list(_format_options(ctx))
    if lines:
        # we use rubric to provide some separation without exploding the table
        # of contents
        yield ".. rubric:: Options"
        yield ""

    yield from lines

    # environment variables

    lines = list(_format_envvars(ctx))
    if lines:
        yield ".. rubric:: Environment variables"
        yield ""

    yield from lines

    # if we're nesting commands, we need to do this slightly differently
    if show_nested:
        return

    commands = _filter_commands(ctx, commands)

    if commands:
        yield ".. rubric:: Commands"
        yield ""
        yield ".. hlist::"
        yield f"  :columns: {len(commands)}"
        yield ""
        for command in commands:
            # for line in _format_subcommand(command):
            #     yield line
            yield f"  * .. object:: {command.name}"
        yield ""


def get_command_docs(name):
    if name in ["tree", "attrs", "cp", "mv", "ls", "ln"]:
        command = getattr(cooler.cli.fileops, name)
    elif name in ["cload pairs", "cload pairix", "cload tabix", "cload hiclib"]:
        command = getattr(cooler.cli.cload, name.split(" ")[1])
    else:
        command = getattr(getattr(cooler.cli, name), name)
    ctx = click.Context(command, info_name="cooler " + name)
    return "\n".join(format_command(ctx, show_nested=False))


COMMANDS = [
    "cload",
    "cload pairs",
    "cload pairix",
    "cload tabix",
    "cload hiclib",
    "load",
    "merge",
    "coarsen",
    "zoomify",
    "balance",
    "info",
    "dump",
    "show",
    "tree",
    "attrs",
    "ls",
    "cp",
    "mv",
    "ln",
    "makebins",
    "digest",
    "csort",
]

SUBCOMMANDS = {"cload": ["pairs", "pairix", "tabix", "hiclib"]}


TEMPLATE = """\
.. _cli-reference:

CLI Reference
=============

.. toctree::
   :maxdepth: 1


Quick reference
---------------

.. program:: cooler
.. code-block:: shell

    cooler [OPTIONS] COMMAND [ARGS]...


.. list-table::
    :widths: 25 100
    :align: left
    :header-rows: 1

    * - Data ingest
      -
    * - `cooler cload`_
      - Create a cooler from genomic point pairs and bins.
    * - `cooler load`_
      - Create a cooler from a pre-binned matrix.

.. list-table::
    :widths: 25 100
    :align: left
    :header-rows: 1

    * - Reduction
      -
    * - `cooler merge`_
      - Merge multiple coolers with identical axes.
    * - `cooler coarsen`_
      - Coarsen a cooler to a lower resolution.
    * - `cooler zoomify`_
      - Generate a multi-resolution cooler file by coarsening.

.. list-table::
    :widths: 25 100
    :align: left
    :header-rows: 1

    * -  Normalization
      -
    * - `cooler balance`_
      - Out-of-core matrix balancing.

.. list-table::
    :widths: 25 100
    :align: left
    :header-rows: 1

    * -  Export/visualization
      -
    * - `cooler info`_
      - Display a cooler`s info and metadata.
    * - `cooler dump`_
      - Dump a cooler`s data to a text stream.
    * - `cooler show`_
      - Display and browse a cooler with matplotlib.

.. list-table::
    :widths: 25 100
    :align: left
    :header-rows: 1

    * -  File manipulation/info
      -
    * - `cooler tree`_
      - Display a file`s data hierarchy.
    * - `cooler attrs`_
      - Display a file`s attribute hierarchy.
    * - `cooler ls`_
      - List all coolers inside a file.
    * - `cooler cp`_
      - Copy a cooler from one file to another or within the same file.
    * - `cooler mv`_
      - Rename a cooler within the same file.
    * - `cooler ln`_
      - Create a hard, soft or external link to a cooler.

.. list-table::
    :widths: 25 100
    :align: left
    :header-rows: 1

    * -  Helper commands
      -
    * - `cooler makebins`_
      - Generate fixed-width genomic bins.
    * - `cooler digest`_
      - Generate fragment-delimited genomic bins.
    * - `cooler csort`_
      - Sort and index a contact list.


.. rubric:: Options

.. option:: -v, --verbose

    Verbose logging.

.. option:: -d, --debug

    On error, drop into the post-mortem debugger shell.

.. option:: -V, --version

    Show the version and exit.

See the cooler_cli.ipynb Jupyter Notebook for specific examples on usage: (https://github.com/open2c/cooler-binder).

----

"""

for cmd in COMMANDS:
    TEMPLATE += f"""\
cooler {cmd}
{"-" * (7 + len(cmd))}

{{{cmd}}}

----

"""


TEMPLATE = TEMPLATE[:-1]  # Drop the very last newline
text = TEMPLATE.format(**{cmd: get_command_docs(cmd) for cmd in COMMANDS})


with open("cli.rst", "w") as f:
    f.write(text)
