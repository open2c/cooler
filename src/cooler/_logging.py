"""
Logging for cooler
==================

Library usage
-------------
As a library, cooler doesn't configure where its log records go: the "cooler"
logger only has a ``NullHandler``, and records propagate to Python's root
logger (``logging.getLogger()``) as usual.

The records are shown by any registered handler along the way: on a "cooler"
logger or on the root logger (e.g. added by ``logging.basicConfig``). If there
aren't any other handlers, the log records are dropped: the ``NullHandler``
stops Python's last-resort handler from printing them. This silencing is
considered good practice for libraries.

If log records are handled by anything else, which ones get logged depends on
the "cooler" logger's level. It isn't set by default, so it inherits the root
logger's value, which is WARNING unless configured. Problems the user can act
on are issued with ``warnings.warn`` rather than a logger, so they are shown
without any logging configuration.

In a script or notebook, you can use :func:`log_to_stderr` to add a stderr
handler to the "cooler" logger to print log messages to stderr. It also sets
the "cooler" logger's level to INFO by default.

CLI usage
---------
The CLI sets up a stderr handler for the duration of each command.
"""
from __future__ import annotations

import logging
import sys
import warnings
from typing import Callable

_LOGGER_NAME = "cooler"
_FORMAT = "{asctime} {levelname:<8} {name}: {message}"
_DATEFMT = "%Y-%m-%d %H:%M:%S"

# Only for the deprecated verbosity functions
_verbosity_to_loglevel = {
    -3: logging.NOTSET,
    -2: logging.CRITICAL,
    -1: logging.ERROR,
    0: logging.WARNING,
    1: logging.INFO,
    2: logging.DEBUG,
}

_loglevel_to_verbosity = {
    logging.NOTSET: -3,
    logging.CRITICAL: -2,
    logging.ERROR: -1,
    logging.WARNING: 0,
    logging.INFO: 1,
    logging.DEBUG: 2,
}

logging.getLogger(_LOGGER_NAME).addHandler(logging.NullHandler())


class _StderrHandler(logging.StreamHandler):
    """
    A handler that writes to whatever ``sys.stderr`` is when a record is
    emitted, rather than the stream at the time it was created. This keeps it
    working when stderr is replaced, e.g. by click's test runner.
    """

    def __init__(self, level: int = logging.NOTSET):
        logging.Handler.__init__(self, level)

    @property
    def stream(self):
        return sys.stderr

    @stream.setter
    def stream(self, value):
        pass


def _make_stderr_handler() -> logging.Handler:
    handler = _StderrHandler()
    handler.setFormatter(logging.Formatter(_FORMAT, _DATEFMT, style="{"))
    return handler


def get_logger(name: str = _LOGGER_NAME) -> logging.Logger:
    """Return a logger in the cooler hierarchy. Same as ``logging.getLogger``."""
    return logging.getLogger(name)


def set_verbosity_level(level: int) -> None:
    """
    Set the level of the "cooler" logger from a verbosity level: -2
    (critical) to 2 (debug), where 0 is warnings.

    .. deprecated:: 0.11.1
        Use :func:`log_to_stderr` or ``logging.getLogger("cooler").setLevel``.
        Will be removed in cooler 0.12.
    """
    warnings.warn(
        "set_verbosity_level is deprecated and will be removed in cooler 0.12. "
        "Use cooler.log_to_stderr(level) or "
        "logging.getLogger('cooler').setLevel(level) instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    logger = logging.getLogger(_LOGGER_NAME)
    try:
        loglevel = _verbosity_to_loglevel[level]
    except KeyError:
        raise ValueError(
            f"Verbosity level must be one of: -2, -1, 0, 1, 2; got '{level}'."
        ) from None
    logger.setLevel(loglevel)


def get_verbosity_level() -> int:
    """
    Get the level of the "cooler" logger as a verbosity level.

    .. deprecated:: 0.11.1
        Use ``logging.getLogger("cooler").level``. Will be removed in cooler
        0.12.
    """
    warnings.warn(
        "get_verbosity_level is deprecated and will be removed in cooler 0.12. "
        "Use logging.getLogger('cooler').level instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    logger = logging.getLogger(_LOGGER_NAME)
    return _loglevel_to_verbosity[logger.level]


def log_to_stderr(level: int | str | None = logging.INFO) -> logging.Handler:
    """
    Print cooler's log messages to stderr.

    Adds a stderr handler to the "cooler" logger, replacing one added by an
    earlier call, and sets the logger's level. Useful in scripts and notebooks
    that don't otherwise configure logging. If Python's root logger also has a
    handler, messages will appear twice.

    Parameters
    ----------
    level : int or str, optional
        Log level, e.g. ``logging.DEBUG`` or "DEBUG". Default is INFO. None
        leaves the level unchanged.

    Returns
    -------
    logging.Handler
        The handler, which can be removed with
        ``logging.getLogger("cooler").removeHandler(handler)``.

    """
    logger = logging.getLogger(_LOGGER_NAME)
    for handler in logger.handlers[:]:
        if getattr(handler, "_cooler_log_to_stderr", False):
            logger.removeHandler(handler)
    handler = _make_stderr_handler()
    handler._cooler_log_to_stderr = True
    logger.addHandler(handler)
    if level is not None:
        logger.setLevel(level)
    return handler


def _log_warning(message, category, filename, lineno, file=None, line=None):
    # One line per warning, unlike logging.captureWarnings, which logs the
    # source file, line number and line of code as well
    logging.getLogger(_LOGGER_NAME).warning(f"{category.__name__}: {message}")


def setup_cli_logging(level: int) -> Callable[[], None]:
    """
    Send log records and warnings to stderr for the duration of a CLI command.

    Returns a function that undoes the setup.
    """
    logger = logging.getLogger(_LOGGER_NAME)
    handler = _make_stderr_handler()
    saved = (logger.level, logger.propagate, warnings.showwarning)

    logger.addHandler(handler)
    logger.setLevel(level)
    # Don't also pass records to Python's root logger, which may have handlers
    logger.propagate = False
    warnings.showwarning = _log_warning

    def restore() -> None:
        logger.removeHandler(handler)
        handler.close()
        logger.setLevel(saved[0])
        logger.propagate = saved[1]
        warnings.showwarning = saved[2]

    return restore
