import atexit
import sys

import click

from .._logging import get_logger, set_logging_context, set_verbosity_level
from .._version import __version__

CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


class UnsortedGroup(click.Group):
    """A click Group that lists commands in the order they were added."""

    def list_commands(self, ctx):
        return list(self.commands)


def setup_psutil_at_exit(logger):
    try:
        import psutil
    except ImportError:
        logger.warning("Install psutil to see process information.")
        return

    attrs_available = {
        x
        for x in dir(psutil.Process)
        if not x.startswith("_")
        and x
        not in {
            "send_signal",
            "suspend",
            "resume",
            "terminate",
            "kill",
            "wait",
            "is_running",
            "as_dict",
            "parent",
            "parents",
            "children",
            "rlimit",
            "memory_info_ex",
            "oneshot",
        }
    }

    attrs = [
        "cmdline",
        "connections",
        "cpu_affinity",
        "cpu_num",
        "cpu_percent",
        "cpu_times",
        "create_time",
        "cwd",
        "environ",
        "exe",
        "gids",
        "io_counters",
        "ionice",
        "memory_full_info",
        "memory_info",
        "memory_maps",
        "memory_percent",
        "name",
        "nice",
        "num_ctx_switches",
        "num_fds",
        "num_threads",
        "open_files",
        "pid",
        "ppid",
        "status",
        "terminal",
        # "threads",  # RuntimeError on MacOS Big Sur
        "uids",
        "username",
    ]

    attrs = [attr for attr in attrs if attr in attrs_available]

    @atexit.register
    def process_dump_at_exit():
        try:
            process = psutil.Process()
            process_info = process.as_dict(attrs, ad_value="")
            for attr in attrs:
                logger.debug(f"PSINFO:'{attr}': {process_info[attr]}")
        except psutil.NoSuchProcess:
            logger.error("PSINFO: Error - Process no longer exists.")


def setup_postmortem_debugger():
    import pdb
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
    set_logging_context("cli")
    set_verbosity_level(min(verbose + 1, 2))
    logger = get_logger()

    if verbose >= 2:  # pragma: no cover
        # Dump process info at exit
        setup_psutil_at_exit(logger)

    if debug:  # pragma: no cover
        # Set hook for postmortem debugging
        setup_postmortem_debugger()


# Load and register cli subcommands
from . import (  # noqa: E402,F401
    balance,
    cload,
    coarsen,
    csort,
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
