import logging
import sys

import click

from .._logging import get_logger, setup_cli_logging
from .._version import __version__

CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


class UnsortedGroup(click.Group):
    """A click Group that lists commands in the order they were added."""

    def list_commands(self, ctx):
        return list(self.commands)


def setup_psutil_at_exit(ctx, logger):
    try:
        import psutil
    except ImportError:
        logger.warning("Install psutil to see process information.")
        return

    # A small set of attributes available on all platforms and stable across
    # psutil releases. Avoid `environ`, which can leak secrets into logs.
    attrs = [
        "pid",
        "ppid",
        "name",
        "cmdline",
        "status",
        "create_time",
        "cpu_times",
        "memory_info",
        "num_threads",
    ]

    # Runs when the command finishes, before its logging is torn down
    @ctx.call_on_close
    def process_dump_at_exit():
        try:
            process = psutil.Process()
            process_info = process.as_dict(attrs, ad_value="")
            for attr in attrs:
                logger.debug(f"PSINFO:'{attr}': {process_info[attr]}")
        except psutil.NoSuchProcess:
            logger.error("PSINFO: Error - Process no longer exists.")


def setup_postmortem_debugger():
    import pdb  # noqa: T100
    import traceback

    def _excepthook(exc_type, value, tb):
        traceback.print_exception(exc_type, value, tb)
        print()
        pdb.pm()

    sys.excepthook = _excepthook


@click.version_option(__version__, "-V", "--version")
@click.group(context_settings=CONTEXT_SETTINGS, cls=UnsortedGroup)
@click.option("-v", "--verbose", help="Verbose logging.", count=True)
@click.option(
    "-d",
    "--debug",
    help="On error, drop into the post-mortem debugger shell.",
    is_flag=True,
    default=False,
)
def cli(verbose, debug):
    """
    Type -h or --help after any subcommand for more information.

    """
    ctx = click.get_current_context()
    # INFO by default, DEBUG with -v
    level = logging.DEBUG if verbose else logging.INFO
    restore_logging = setup_cli_logging(level)
    ctx.call_on_close(restore_logging)
    logger = get_logger()

    if verbose >= 2:
        # Dump process info when the command finishes
        setup_psutil_at_exit(ctx, logger)

    if debug:
        # Set hook for postmortem debugging
        setup_postmortem_debugger()


# Load and register cli subcommands
from . import (  # noqa: E402,F401
    balance,
    cload,
    coarsen,
    digest,
    dump,
    fileops,
    info,
    load,
    makebins,
    merge,
    show,
    zoomify,
)
