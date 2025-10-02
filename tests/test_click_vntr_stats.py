"""Tests for muconeup analyze vntr-stats CLI command.

Tests cover:
- Command help and options
- Integration with Click CLI
- File I/O and output formatting
- Error handling and exit codes
- Real-world usage scenarios
"""

import json

import pytest
from click.testing import CliRunner

from muc_one_up.cli.click_main import cli


@pytest.fixture
def runner():
    """Click test runner."""
    return CliRunner()


@pytest.fixture
def config_with_repeats(tmp_path):
    """Create config file with known repeats."""
    config = {
        "repeats": {
            "1": "ACGT",
            "2": "GCTA",
            "3": "TGCA",
            "X": "ATCG",
            "A": "CGAT",
            "B": "TAGC",
            "6p": "CTAG",
            "7": "GATC",
            "8": "TCAG",
            "9": "AGCT",
        },
        "constants": {"left": "AAA", "right": "TTT"},
    }
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps(config, indent=2))
    return str(config_file)


@pytest.fixture
def vntr_data_with_header(tmp_path):
    """Create test VNTR data file with header."""
    data = """publication\tid\tallele\tvntr
PMID:123\tS1\t1\t1-2-3-X-A-B
PMID:123\tS2\t2\tX-A-B-6p-7
PMID:456\tS3\t1\t1-2-3-X
"""
    data_file = tmp_path / "vntr_data.tsv"
    data_file.write_text(data)
    return str(data_file)


@pytest.fixture
def vntr_data_no_header(tmp_path):
    """Create test VNTR data file without header."""
    data = """1-2-3-X
X-A-B
1-2-3
"""
    data_file = tmp_path / "vntr_noheader.tsv"
    data_file.write_text(data)
    return str(data_file)


@pytest.mark.unit
@pytest.mark.cli
class TestVNTRStatsCommand:
    """Tests for vntr-stats subcommand."""

    def test_command_help(self, runner, config_with_repeats):
        """Test vntr-stats help displays correctly."""
        result = runner.invoke(
            cli, ["--config", config_with_repeats, "analyze", "vntr-stats", "--help"]
        )
        assert result.exit_code == 0
        assert "transition probabilities" in result.output.lower()
        assert "--structure-column" in result.output
        assert "--header" in result.output
        assert "--delimiter" in result.output
        assert "--output" in result.output

    def test_analyze_with_header(self, runner, config_with_repeats, vntr_data_with_header):
        """Test analysis with header row."""
        result = runner.invoke(
            cli,
            [
                "--config",
                config_with_repeats,
                "analyze",
                "vntr-stats",
                vntr_data_with_header,
                "--header",
                "--structure-column",
                "vntr",
            ],
        )

        assert result.exit_code == 0
        # Extract JSON from output (skip log lines before and after)
        lines = result.output.strip().split("\n")
        json_start = next(i for i, line in enumerate(lines) if line.strip().startswith("{"))
        # Find the closing brace (line that is exactly "}")
        json_end = json_start
        for i in range(len(lines) - 1, json_start - 1, -1):
            if lines[i] == "}":
                json_end = i
                break
        json_output = "\n".join(lines[json_start : json_end + 1])
        output = json.loads(json_output)
        assert "min_repeats" in output
        assert "max_repeats" in output
        assert "mean_repeats" in output
        assert "median_repeats" in output
        assert "probabilities" in output
        assert "repeats" in output

    def test_analyze_without_header(self, runner, config_with_repeats, vntr_data_no_header):
        """Test analysis without header (column index)."""
        result = runner.invoke(
            cli,
            [
                "--config",
                config_with_repeats,
                "analyze",
                "vntr-stats",
                vntr_data_no_header,
                "--structure-column",
                "0",
            ],
        )

        assert result.exit_code == 0
        # Extract JSON from output (skip log lines before and after)
        lines = result.output.strip().split("\n")
        json_start = next(i for i, line in enumerate(lines) if line.strip().startswith("{"))
        # Find the closing brace (line that is exactly "}")
        json_end = json_start
        for i in range(len(lines) - 1, json_start - 1, -1):
            if lines[i] == "}":
                json_end = i
                break
        json_output = "\n".join(lines[json_start : json_end + 1])
        output = json.loads(json_output)
        assert output["min_repeats"] == 3
        assert output["max_repeats"] == 4

    def test_output_to_file(self, runner, config_with_repeats, vntr_data_with_header, tmp_path):
        """Test writing output to file."""
        output_file = tmp_path / "stats.json"

        result = runner.invoke(
            cli,
            [
                "--config",
                config_with_repeats,
                "analyze",
                "vntr-stats",
                vntr_data_with_header,
                "--header",
                "--output",
                str(output_file),
            ],
        )

        assert result.exit_code == 0
        assert output_file.exists()

        with output_file.open() as f:
            data = json.load(f)
        assert "probabilities" in data
        assert "repeats" in data

    def test_output_short_flag(self, runner, config_with_repeats, vntr_data_with_header, tmp_path):
        """Test -o short flag for output."""
        output_file = tmp_path / "stats.json"

        result = runner.invoke(
            cli,
            [
                "--config",
                config_with_repeats,
                "analyze",
                "vntr-stats",
                vntr_data_with_header,
                "--header",
                "-o",
                str(output_file),
            ],
        )

        assert result.exit_code == 0
        assert output_file.exists()

    def test_custom_delimiter(self, runner, config_with_repeats, tmp_path):
        """Test analysis with custom delimiter."""
        data = """publication,id,vntr
PMID:123,S1,1-2-3
"""
        data_file = tmp_path / "comma.csv"
        data_file.write_text(data)

        result = runner.invoke(
            cli,
            [
                "--config",
                config_with_repeats,
                "analyze",
                "vntr-stats",
                str(data_file),
                "--header",
                "--delimiter",
                ",",
                "--structure-column",
                "vntr",
            ],
        )

        assert result.exit_code == 0
        # Extract JSON from output (skip log lines before and after)
        lines = result.output.strip().split("\n")
        json_start = next(i for i, line in enumerate(lines) if line.strip().startswith("{"))
        # Find the closing brace (line that is exactly "}")
        json_end = json_start
        for i in range(len(lines) - 1, json_start - 1, -1):
            if lines[i] == "}":
                json_end = i
                break
        json_output = "\n".join(lines[json_start : json_end + 1])
        output = json.loads(json_output)
        assert "min_repeats" in output

    def test_nonexistent_file_error(self, runner, config_with_repeats):
        """Test error handling for nonexistent file."""
        result = runner.invoke(
            cli,
            [
                "--config",
                config_with_repeats,
                "analyze",
                "vntr-stats",
                "/nonexistent/file.tsv",
            ],
        )

        assert result.exit_code != 0
        assert "does not exist" in result.output.lower() or "error" in result.output.lower()

    def test_empty_file_error(self, runner, config_with_repeats, tmp_path):
        """Test error handling for empty file."""
        empty_file = tmp_path / "empty.tsv"
        empty_file.write_text("vntr\n")

        result = runner.invoke(
            cli,
            [
                "--config",
                config_with_repeats,
                "analyze",
                "vntr-stats",
                str(empty_file),
                "--header",
            ],
        )

        assert result.exit_code == 1
        assert "failed" in result.output.lower() or "no valid" in result.output.lower()

    def test_config_without_repeats_error(self, runner, tmp_path, vntr_data_with_header):
        """Test error when config lacks repeats key."""
        bad_config = tmp_path / "bad_config.json"
        bad_config.write_text(json.dumps({"constants": {"left": "A"}}))

        result = runner.invoke(
            cli,
            [
                "--config",
                str(bad_config),
                "analyze",
                "vntr-stats",
                vntr_data_with_header,
                "--header",
            ],
        )

        assert result.exit_code == 1
        assert "repeats" in result.output.lower()

    def test_default_structure_column(self, runner, config_with_repeats, tmp_path):
        """Test that default structure column is 'vntr'."""
        data = """vntr
1-2-3
X-A-B
"""
        data_file = tmp_path / "default_col.tsv"
        data_file.write_text(data)

        result = runner.invoke(
            cli,
            [
                "--config",
                config_with_repeats,
                "analyze",
                "vntr-stats",
                str(data_file),
                "--header",
                # Not specifying --structure-column, should default to "vntr"
            ],
        )

        assert result.exit_code == 0
        # Extract JSON from output (skip log lines before and after)
        lines = result.output.strip().split("\n")
        json_start = next(i for i, line in enumerate(lines) if line.strip().startswith("{"))
        # Find the closing brace (line that is exactly "}")
        json_end = json_start
        for i in range(len(lines) - 1, json_start - 1, -1):
            if lines[i] == "}":
                json_end = i
                break
        json_output = "\n".join(lines[json_start : json_end + 1])
        output = json.loads(json_output)
        assert output["min_repeats"] == 3


@pytest.mark.unit
@pytest.mark.cli
class TestVNTRStatsIntegration:
    """Integration tests for vntr-stats command."""

    def test_real_world_vntr_database(self, runner):
        """Test with realistic VNTR database (uses actual example file if available)."""
        # This test uses the actual example file
        from pathlib import Path

        example_file = Path("data/examples/vntr_database.tsv")
        config_file = Path("config.json")

        if not example_file.exists() or not config_file.exists():
            pytest.skip("Example data or config not available")

        result = runner.invoke(
            cli,
            [
                "--config",
                str(config_file),
                "analyze",
                "vntr-stats",
                str(example_file),
                "--header",
                "--structure-column",
                "vntr",
            ],
        )

        assert result.exit_code == 0
        # Extract JSON from output (skip log lines before and after)
        lines = result.output.strip().split("\n")
        json_start = next(i for i, line in enumerate(lines) if line.strip().startswith("{"))
        # Find the closing brace (line that is exactly "}")
        json_end = json_start
        for i in range(len(lines) - 1, json_start - 1, -1):
            if lines[i] == "}":
                json_end = i
                break
        json_output = "\n".join(lines[json_start : json_end + 1])
        output = json.loads(json_output)

        # Verify expected structure
        assert "min_repeats" in output
        assert "max_repeats" in output
        assert "mean_repeats" in output
        assert "probabilities" in output
        assert "repeats" in output

        # Verify reasonable values (based on known data)
        assert output["min_repeats"] >= 10  # Real VNTR data has long sequences
        assert output["max_repeats"] <= 200
        assert "END" in str(output["probabilities"])  # Should have END state

    def test_transition_matrix_structure(self, runner, config_with_repeats, tmp_path):
        """Test that transition matrix has correct structure."""
        data = """vntr
1-2-3
X-A-B
"""
        data_file = tmp_path / "matrix_test.tsv"
        data_file.write_text(data)

        result = runner.invoke(
            cli,
            [
                "--config",
                config_with_repeats,
                "analyze",
                "vntr-stats",
                str(data_file),
                "--header",
            ],
        )

        assert result.exit_code == 0
        # Extract JSON from output (skip log lines before and after)
        lines = result.output.strip().split("\n")
        json_start = next(i for i, line in enumerate(lines) if line.strip().startswith("{"))
        # Find the closing brace (line that is exactly "}")
        json_end = json_start
        for i in range(len(lines) - 1, json_start - 1, -1):
            if lines[i] == "}":
                json_end = i
                break
        json_output = "\n".join(lines[json_start : json_end + 1])
        output = json.loads(json_output)

        probs = output["probabilities"]

        # Check transitions exist
        assert "1" in probs
        assert "2" in probs["1"]  # 1 -> 2 transition
        assert "END" in probs["3"]  # 3 -> END transition
        assert "END" in probs["B"]  # B -> END transition

        # Check probabilities sum to 1 for each source
        for _source, targets in probs.items():
            total = sum(targets.values())
            assert abs(total - 1.0) < 0.001  # Within floating point precision
