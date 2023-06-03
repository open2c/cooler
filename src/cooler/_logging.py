import logging
import sys

_logging_context = None
_loggers = {}

verbosity_to_loglevel = {
    -3: logging.NOTSET,
    -2: logging.CRITICAL,
    -1: logging.ERROR,
    0: logging.WARNING,
    1: logging.INFO,
    2: logging.DEBUG,
}

loglevel_to_verbosity = {
    logging.NOTSET: -3,
    logging.CRITICAL: -2,
    logging.ERROR: -1,
    logging.WARNING: 0,
    logging.INFO: 1,
    logging.DEBUG: 2,
}


def configure(
    logger,
    stream=None,
    filename=None,
    open_kws=None,
    handlers=None,
    level=logging.WARNING,
    format="{levelname}:{name}:{message}",
    datefmt="%Y-%m-%d %I:%M:%S %p",
    style="{",
    propagate=False,
):
    """
    Configure a logger for a stream or file or a custom set of handlers.

    Based on `logging.basicConfig` but works on any logger, not just the root one.

    Parameters
    ----------
    logger : :class:`logging.Logger`
        A logger.
    stream : file-like, optional
        A stream, like ``sys.stdout``.
    filename : str, optional
        Path to a file, instead of a stream.
    open_kws : dict, optionsl
        Keyword args to ``open`` if using a filename instead of a stream.
        Defaults to using mode='a' and for text streams: encoding='utf-8' and
        errors='backslashreplace'.
    handlers : sequence of logging Handlers
        Arbitrary logging handlers to register. Cannot be used with ``filename``
        or ``stream``.
    level : int, optional
        The log level.
    format : str, optional
        A format string for log records.
    datefmt : str, optional
        A format string for the ``asctime`` variable when used in log records.
    style : {"{", "$", "%"}, optional
        The templating style of the log record format string.
    propagate : bool, optional
        Whether the logger should propagate log records up to its parent.

    Notes
    -----
    Any of the logger's existing handlers will be closed and destroyed.

    For logging level values, see
    https://docs.python.org/3/howto/logging.html#logging-levels

    For a list of variables that can go into log records, see
    https://docs.python.org/3/library/logging.html#logrecord-attributes

    """
    if handlers is None:
        if stream is not None and filename is not None:
            raise ValueError(
                "'stream' and 'filename' should not be specified together"
            )
    else:
        if stream is not None or filename is not None:
            raise ValueError(
                "'stream' or 'filename' should not be specified together with "
                "'handlers'"
            )

    # Set up the new handlers
    if filename is not None:
        open_kws = {} if open_kws is None else open_kws
        open_kws.setdefault("mode", "a")
        if "b" in open_kws["mode"]:
            open_kws["encoding"] = None
            open_kws["errors"] = None
        else:
            open_kws.setdefault("encoding", "utf-8")
            open_kws.setdefault("errors", "backslashreplace")
        handlers = [logging.FileHandler(filename, **open_kws)]
    elif stream is not None:
        handlers = [logging.StreamHandler(stream)]

    # Wipe out the old handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
        handler.close()

    # Set up the formatter and register it on all new handlers
    formatter = logging.Formatter(format, datefmt, style)
    for handler in handlers:
        if handler.formatter is None:
            handler.setFormatter(formatter)
        logger.addHandler(handler)

    if level is not None:
        logger.setLevel(level)

    logger.propagate = propagate


def set_logging_context(ctx):
    global _logging_context

    logger = logging.getLogger("cooler")

    if _logging_context != ctx:
        if ctx == "lib":
            configure(
                logger, stream=sys.stdout, level=logging.WARNING, format="{message}"
            )
            logging.captureWarnings(False)
        elif ctx == "cli":
            configure(logger, stream=sys.stderr, level=logging.INFO)
            logging.captureWarnings(True)
        elif ctx == "none":
            for handler in logger.handlers[:]:
                logger.removeHandler(handler)
                handler.close()
            logging.captureWarnings(False)
        else:
            raise ValueError(f"Unknown logging context: '{ctx}'")
        _logging_context = ctx


def get_logging_context():
    return _logging_context


def set_verbosity_level(level):
    logger = logging.getLogger("cooler")
    try:
        loglevel = verbosity_to_loglevel[level]
    except KeyError:
        raise ValueError(
            f"Verbosity level must be one of: -2, -1, 0, 1, 2; got '{level}'."
        ) from None
    logger.setLevel(loglevel)


def get_verbosity_level():
    logger = logging.getLogger("cooler")
    return loglevel_to_verbosity[logger.level]


def get_logger(name="cooler"):
    # Based on ipython traitlets
    global _loggers, _logging_context

    if _logging_context is None:
        set_logging_context("lib")

    if name not in _loggers:
        _loggers[name] = logging.getLogger(name)
        # Add a NullHandler to silence warnings about not being
        # initialized, per best practice for libraries.
        _loggers[name].addHandler(logging.NullHandler())

    return _loggers[name]
