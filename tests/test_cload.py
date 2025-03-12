from click.testing import CliRunner

from cooler.cli.cload import pairs

runner = CliRunner()


def test_valid_indices(tmpdir):
    """Test that valid column indices process successfully."""
    bins_path = tmpdir.join("bins.bed")
    pairs_path = tmpdir.join("pairs.txt")
    out_path = tmpdir.join("out.cool")

    bins_path.write("chr1\t0\t1000000\nchr1\t1000000\t2000000\nchr2\t0\t1000000\n")
    pairs_path.write("chr1\t1000\tchr1\t2000\nchr1\t1500\tchr2\t3000\n")

    result = runner.invoke(pairs, [
        "-c1", "1",
        "-p1", "2",
        "-c2", "3",
        "-p2", "4",
        str(bins_path),
        str(pairs_path),
        str(out_path),
    ])

    assert result.exit_code == 0, f"Command failed: {result.output}"


def test_invalid_index(tmpdir):
    """Test that an invalid column index raises an error."""
    bins_path = tmpdir.join("bins.bed")
    pairs_path = tmpdir.join("pairs.txt")
    out_path = tmpdir.join("out.cool")

    bins_path.write("chr1\t0\t1000000\n")
    pairs_path.write("chr1\t1000\tchr1\t2000\n")

    result = runner.invoke(pairs, [
        "-c1", "1",
        "-p1", "2",
        "-c2", "3",
        "-p2", "5",
        str(bins_path),
        str(pairs_path),
        str(out_path),
    ])

    assert result.exit_code != 0, "Command should have failed"
    assert isinstance(result.exception, ValueError), (
        f"Expected ValueError, got {result.exception}"
    )
    assert (
        "Column index 5 (pos2) exceeds number of columns (4)" in str(result.exception)
    ), (
        f"Unexpected error: {result.exception}"
    )


def test_stdin_valid(tmpdir):
    """Test valid indices with stdin input."""
    bins_path = tmpdir.join("bins.bed")
    out_path = tmpdir.join("out.cool")

    bins_path.write("chr1\t0\t1000000\n")
    pairs_data = "chr1\t1000\tchr1\t2000\n"

    result = runner.invoke(
        pairs,
        [
            "-c1", "1",
            "-p1", "2",
            "-c2", "3",
            "-p2", "4",
            str(bins_path),
            "-",
            str(out_path),
        ],
        input=pairs_data,
    )

    assert result.exit_code == 0, f"Command failed: {result.output} {result.exception}"


def test_stdin_invalid(tmpdir):
    """Test invalid index with stdin raises an error."""
    bins_path = tmpdir.join("bins.bed")
    out_path = tmpdir.join("out.cool")

    bins_path.write("chr1\t0\t1000000\n")
    pairs_data = "chr1\t1000\tchr1\t2000\n"

    result = runner.invoke(
        pairs,
        [
            "-c1", "1",
            "-p1", "2",
            "-c2", "3",
            "-p2", "5",
            str(bins_path),
            "-",
            str(out_path),
        ],
        input=pairs_data,
    )

    assert result.exit_code != 0, "Command should have failed"
    assert isinstance(result.exception, ValueError), (
        f"Expected ValueError, got {result.exception}"
    )
    assert (
        "Column index 5 (pos2) exceeds number of columns (4)" in str(result.exception)
    ), (
        f"Unexpected error: {result.exception}"
    )


def test_header_skipping(tmpdir):
    """Test that headers are skipped and validation works."""
    bins_path = tmpdir.join("bins.bed")
    pairs_path = tmpdir.join("pairs.txt")
    out_path = tmpdir.join("out.cool")

    bins_path.write("chr1\t0\t1000000\n")
    pairs_path.write("#comment1\n#comment2\nchr1\t1000\tchr1\t2000\n")

    result = runner.invoke(pairs, [
        "-c1", "1",
        "-p1", "2",
        "-c2", "3",
        "-p2", "4",
        str(bins_path),
        str(pairs_path),
        str(out_path),
    ])

    assert result.exit_code == 0, f"Command failed: {result.output}"


def test_empty_file(tmpdir):
    """Test that an empty file or header-only file raises an error."""
    bins_path = tmpdir.join("bins.bed")
    pairs_path = tmpdir.join("pairs.txt")
    out_path = tmpdir.join("out.cool")

    bins_path.write("chr1\t0\t1000000\n")
    pairs_path.write("#comment1\n#comment2\n")

    result = runner.invoke(pairs, [
        "-c1", "1",
        "-p1", "2",
        "-c2", "3",
        "-p2", "4",
        str(bins_path),
        str(pairs_path),
        str(out_path),
    ])

    assert result.exit_code != 0, "Command should have failed"
    assert isinstance(result.exception, ValueError), (
        f"Expected ValueError, got {result.exception}"
    )
    assert "is empty or contains only header lines" in str(result.exception), (
        f"Unexpected error: {result.exception}"
    )


def test_extra_field(tmpdir):
    """Test validation with an additional --field column."""
    bins_path = tmpdir.join("bins.bed")
    pairs_path = tmpdir.join("pairs.txt")
    out_path = tmpdir.join("out.cool")

    bins_path.write("chr1\t0\t1000000\nchr2\t0\t1000000\n")
    pairs_path.write("chr1\t1000\tchr1\t2000\t1.5\nchr1\t1500\tchr2\t3000\t2.0\n")

    result = runner.invoke(pairs, [
        "-c1", "1",
        "-p1", "2",
        "-c2", "3",
        "-p2", "4",
        "--field", "value=5",
        str(bins_path),
        str(pairs_path),
        str(out_path),
    ])

    assert result.exit_code == 0, f"Command failed: {result.output}"

    # Test invalid extra field
    result = runner.invoke(pairs, [
        "-c1", "1",
        "-p1", "2",
        "-c2", "3",
        "-p2", "4",
        "--field", "value=6",
        str(bins_path),
        str(pairs_path),
        str(out_path),
    ])

    assert result.exit_code != 0, "Command should have failed"
    assert isinstance(result.exception, ValueError), (
        f"Expected ValueError, got {result.exception}"
    )
    assert (
        "Column index 6 (value) exceeds number of columns (5)" in str(result.exception)
    ), (
        f"Unexpected error: {result.exception}"
    )
