"""Tests for muconeup analyze snapshot-validate CLI command.

Tests cover:
- Command help and options
- Integration with Click CLI
- Input validation and error handling
- JSON output formatting
- Exit codes (0 on successful execution)
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
def config_with_snapshot(tmp_path):
    """Create minimal config with snapshot_validation section."""
    config = {
        "repeats": {"X": "GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA"},
        "constants": {"left": "AAA", "right": "TTT"},
        "probabilities": {"X": {"X": 1.0}},
        "length_model": {
            "distribution": "normal",
            "min_repeats": 10,
            "max_repeats": 30,
            "mean_repeats": 20,
            "median_repeats": 20,
        },
        "mutations": {},
        "tools": {"samtools": "samtools"},
        "read_simulation": {"human_reference": "/dev/null", "threads": 1},
        "snapshot_validation": {
            "dupC": {
                "description": "Test 8C mutation",
                "pcr": {
                    "forward_primer": "GGCCGGCCCCGGGCTCCACC",
                    "reverse_primer": "TGTCACCTCGGCCCCGGA",
                    "reverse_needs_rc": False,
                    "max_products": 100,
                    "size_range": {"min": 50, "max": 65},
                },
                "digest": {
                    "enzyme": "MwoI",
                    "recognition_site": "GCNNNNNNNGC",
                    "expected_survivors": "mutant_only",
                },
                "snapshot": {
                    "primers": {
                        "primer_7c": "CGGGCTCCACCGCCCCCCC",
                        "primer_repeat_r": "TGTCACCTCGGCCCCGGA",
                    },
                    "fluorophore_map": {
                        "A": {"color": "Green", "dye": "dR6G"},
                        "C": {"color": "Black", "dye": "dTAMRA"},
                        "G": {"color": "Blue", "dye": "dR110"},
                        "T": {"color": "Red", "dye": "dROX"},
                    },
                },
                "validation": {
                    "mutant_pattern": "GCCCCCCCCAGC",
                    "normal_pattern": "GCCCCCCCAGC",
                    "expected_mutant_fluorescence": "Black",
                },
            }
        },
    }
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps(config, indent=2))
    return str(config_file)


@pytest.fixture
def sample_fasta(tmp_path):
    """Create simple test FASTA file."""
    fasta_content = """>haplotype_1
AAAGCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCATTT
>haplotype_2
AAAGCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCATTT
"""
    fasta_file = tmp_path / "sample.fa"
    fasta_file.write_text(fasta_content)
    return str(fasta_file)


@pytest.mark.unit
@pytest.mark.cli
class TestSnapshotValidateCommand:
    """Tests for snapshot-validate subcommand."""

    def test_command_help(self, runner, config_with_snapshot):
        """Test snapshot-validate help displays correctly."""
        result = runner.invoke(
            cli, ["--config", config_with_snapshot, "analyze", "snapshot-validate", "--help"]
        )
        assert result.exit_code == 0
        assert "Validate SNaPshot assay" in result.output
        assert "--mutation" in result.output
        assert "--output" in result.output

    def test_missing_input_file(self, runner, config_with_snapshot):
        """Test error when input file is missing."""
        result = runner.invoke(
            cli,
            [
                "--config",
                config_with_snapshot,
                "analyze",
                "snapshot-validate",
                "nonexistent.fa",
                "--mutation",
                "dupC",
            ],
        )
        assert result.exit_code != 0
        assert "does not exist" in result.output or "Error" in result.output

    def test_missing_mutation_option(self, runner, config_with_snapshot, sample_fasta):
        """Test error when --mutation option is missing."""
        result = runner.invoke(
            cli, ["--config", config_with_snapshot, "analyze", "snapshot-validate", sample_fasta]
        )
        assert result.exit_code != 0
        assert "Error" in result.output or "required" in result.output.lower()

    def test_unknown_mutation(self, runner, config_with_snapshot, sample_fasta):
        """Test error when mutation not in config."""
        result = runner.invoke(
            cli,
            [
                "--config",
                config_with_snapshot,
                "analyze",
                "snapshot-validate",
                sample_fasta,
                "--mutation",
                "unknownMutation",
            ],
        )
        assert result.exit_code == 1
        # Command should handle this gracefully

    def test_valid_run_json_output(self, runner, config_with_snapshot, sample_fasta):
        """Test successful run with JSON output to stdout."""
        result = runner.invoke(
            cli,
            [
                "--config",
                config_with_snapshot,
                "analyze",
                "snapshot-validate",
                sample_fasta,
                "--mutation",
                "dupC",
            ],
        )
        # Should succeed
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
        output_data = json.loads(json_output)
        assert "mutation" in output_data
        assert "haplotypes" in output_data
        assert "overall_detection" in output_data
        assert output_data["mutation"] == "dupC"

    def test_valid_run_json_to_file(self, runner, config_with_snapshot, sample_fasta, tmp_path):
        """Test successful run with JSON output to file."""
        output_file = tmp_path / "results.json"

        result = runner.invoke(
            cli,
            [
                "--config",
                config_with_snapshot,
                "analyze",
                "snapshot-validate",
                sample_fasta,
                "--mutation",
                "dupC",
                "--output",
                str(output_file),
            ],
        )

        # Should succeed
        assert result.exit_code == 0

        # Output file should exist and contain valid JSON
        assert output_file.exists()
        with open(output_file) as f:
            output_data = json.load(f)
            assert "mutation" in output_data
            assert "haplotypes" in output_data
            assert output_data["mutation"] == "dupC"

    def test_json_structure(self, runner, config_with_snapshot, sample_fasta):
        """Test that JSON output has expected structure."""
        result = runner.invoke(
            cli,
            [
                "--config",
                config_with_snapshot,
                "analyze",
                "snapshot-validate",
                sample_fasta,
                "--mutation",
                "dupC",
            ],
        )

        assert result.exit_code == 0
        # Extract JSON from output (skip log lines before and after)
        lines = result.output.strip().split("\n")
        json_start = next(i for i, line in enumerate(lines) if line.strip().startswith("{"))
        json_end = json_start
        for i in range(len(lines) - 1, json_start - 1, -1):
            if lines[i] == "}":
                json_end = i
                break
        json_output = "\n".join(lines[json_start : json_end + 1])
        output_data = json.loads(json_output)

        # Top-level keys
        assert "mutation" in output_data
        assert "input_file" in output_data
        assert "haplotypes" in output_data
        assert "overall_detection" in output_data

        # Haplotype structure
        for _hap_name, hap_data in output_data["haplotypes"].items():
            assert "mutation_detected" in hap_data
            assert "pcr_products" in hap_data
            assert "digest_survivors" in hap_data
            assert "summary" in hap_data

    def test_exit_code_logic(self, runner, config_with_snapshot, sample_fasta, tmp_path):
        """Test that exit code is 0 on successful execution for negative and positive controls."""
        # 1. Negative control (wild-type sample, mutation not detected)
        result_neg = runner.invoke(
            cli,
            [
                "--config",
                config_with_snapshot,
                "analyze",
                "snapshot-validate",
                sample_fasta,
                "--mutation",
                "dupC",
            ],
        )
        assert result_neg.exit_code == 0, "Exit code should be 0 on negative control"
        lines_neg = result_neg.output.strip().split("\n")
        json_start_neg = next(i for i, line in enumerate(lines_neg) if line.strip().startswith("{"))
        json_end_neg = max(i for i, line in enumerate(lines_neg) if line.strip() == "}")
        data_neg = json.loads("\n".join(lines_neg[json_start_neg : json_end_neg + 1]))
        assert data_neg["overall_detection"] is False

        # 2. Positive control (mutant sample with 8C insertion, mutation detected)
        mutant_fasta = tmp_path / "mutant.fa"
        mutant_seq = "GGCCGGCCCCGGGCTCCACCGCCCCCCCCAGCCCACGGTTGTCACCTCGGCCCCGGA"
        mutant_fasta.write_text(f">haplotype_1\n{mutant_seq}\n")
        result_pos = runner.invoke(
            cli,
            [
                "--config",
                config_with_snapshot,
                "analyze",
                "snapshot-validate",
                str(mutant_fasta),
                "--mutation",
                "dupC",
            ],
        )
        assert result_pos.exit_code == 0, "Exit code should be 0 when mutation detected"
        lines_pos = result_pos.output.strip().split("\n")
        json_start_pos = next(i for i, line in enumerate(lines_pos) if line.strip().startswith("{"))
        json_end_pos = max(i for i, line in enumerate(lines_pos) if line.strip() == "}")
        data_pos = json.loads("\n".join(lines_pos[json_start_pos : json_end_pos + 1]))
        assert data_pos["overall_detection"] is True
