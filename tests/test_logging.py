"""
Tests for logging in the library and the CLI.
"""
import logging
import os.path as op
import warnings

from click.testing import CliRunner

import cooler
from cooler._logging import get_logger
from cooler.cli import cli

testdir = op.dirname(op.realpath(__file__))
COOL = op.join(testdir, "data", "hg19.GM12878-MboI.matrix.2000kb.cool")


def _logging_state():
    logger = logging.getLogger("cooler")
    return (
        list(logger.handlers),
        logger.level,
        logger.propagate,
        warnings.showwarning,
    )


def test_library_does_not_configure_output(caplog):
    logger = logging.getLogger("cooler")
    assert [type(h) for h in logger.handlers] == [logging.NullHandler]
    assert logger.propagate
    # Records reach the root logger, so the user's configuration applies
    with caplog.at_level(logging.INFO, logger="cooler"):
        get_logger("cooler.test").info("hello")
    assert [r.getMessage() for r in caplog.records] == ["hello"]


def test_log_to_stderr(capsys):
    logger = logging.getLogger("cooler")
    level = logger.level
    try:
        cooler.log_to_stderr()
        handler = cooler.log_to_stderr("DEBUG")  # replaces the first one
        get_logger("cooler.test").debug("hello")
        captured = capsys.readouterr()
        assert captured.err.count("hello") == 1
        assert captured.out == ""
    finally:
        logger.removeHandler(handler)
        logger.setLevel(level)


def test_verbosity_level():
    level = logging.getLogger("cooler").level
    try:
        cooler.set_verbosity_level(2)
        assert cooler.get_verbosity_level() == 2
        assert logging.getLogger("cooler").level == logging.DEBUG
    finally:
        logging.getLogger("cooler").setLevel(level)


def test_cli_logs_to_stderr_and_cleans_up():
    before = _logging_state()
    try:
        runner = CliRunner(mix_stderr=False)  # click < 8.2 mixes them by default
    except TypeError:
        runner = CliRunner()  # click >= 8.2 always captures stderr separately
    args = ["balance", COOL, "--stdout", "-p", "1", "--max-iters", "1"]
    # Twice in one process: the second run must not log to a stale stream
    for _ in range(2):
        result = runner.invoke(cli, args)
        assert result.exit_code == 0, result.output
        assert "INFO     cooler.cli.balance: Balancing" in result.stderr
        # Library warnings are logged too, through the same handler
        warning = "WARNING  cooler: ConvergenceWarning: Iteration limit reached"
        assert warning in result.stderr
        assert "INFO" not in result.stdout
        assert "Iteration limit" not in result.stdout
        assert _logging_state() == before
