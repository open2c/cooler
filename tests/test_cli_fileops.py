import os.path as op

from click.testing import CliRunner

from cooler.cli import cli

testdir = op.realpath(op.dirname(__file__))
datadir = op.join(testdir, "data")


def test_cp():
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(
            cli,
            [
                "cp",
                op.join(datadir, "toy.symm.upper.2.cool"),
                "test.cool",
            ],
        )
        assert result.exit_code == 0

        result = runner.invoke(
            cli,
            [
                "cp",
                "test.cool",
                "test2.cool::some/old_path",
            ],
        )
        assert result.exit_code == 0

        result = runner.invoke(
            cli,
            [
                "mv",
                "test2.cool::some/old_path",
                "test2.cool::some/path",
            ],
        )
        assert result.exit_code == 0

        # mv only renames within a file
        result = runner.invoke(cli, ["mv", "test.cool", "test4.cool::x"])
        assert result.exit_code != 0
        assert not op.exists("test4.cool")

        result = runner.invoke(
            cli,
            [
                "ln",
                "test2.cool::some/path",
                "test2.cool::hard/link",
            ],
        )
        assert result.exit_code == 0

        result = runner.invoke(
            cli,
            [
                "ln",
                "-s",
                "test2.cool::some/path",
                "test2.cool::soft/link",
            ],
        )
        assert result.exit_code == 0

        result = runner.invoke(
            cli,
            [
                "ln",
                "-s",
                "test2.cool::some/path",
                "test3.cool::ext/link",
            ],
        )
        assert result.exit_code == 0


def test_list_coolers():
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "ls",
            op.join(datadir, "toy.symm.upper.2.cool"),
        ],
    )
    assert result.exit_code == 0

    result = runner.invoke(
        cli,
        [
            "ls",
            "-l",
            op.join(datadir, "toy.symm.upper.2.cool"),
        ],
    )
    assert result.exit_code == 0


def test_tree():
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "tree",
            op.join(datadir, "toy.symm.upper.2.cool"),
        ],
    )
    assert result.exit_code == 0


def test_attrs():
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "attrs",
            op.join(datadir, "toy.symm.upper.2.cool"),
        ],
    )
    assert result.exit_code == 0
