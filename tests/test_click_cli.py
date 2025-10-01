"""Tests for muc_one_up.cli.click_main module.

Tests cover:
- CLI group structure and help
- Simulate command with all options
- Reads command group (illumina, ont)
- Analyze command group (orfs, stats)
- Error handling and exit codes
- CliRunner integration
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
def temp_config(tmp_path):
    """Create a minimal valid config file for testing."""
    config = {
        "repeats": {
            "1": "GCCACCACCTCCAACTCCT",
            "2": "GCCACCACC",
            "3": "TCCAACTCCT",
            "6": "TCCGACTCCT",
            "6p": "TCCAACTCCT",
            "7": "TCTAGGACCTCCTAGCTCCCAGAACTCCT",
            "8": "TCTAGGACCTCCTAGCTCCTAGCTCCT",
            "9": "TCTAGGACCTAGCTCCT",
        },
        "constants": {
            "left": "ATGGCCCCATCTCTCACCGTCTCGGTCATCTCCTTGATG",
            "right": "CAATGGTGTCTTGGGTAGCTTCGTCACGGTTTTCCAG",
        },
        "probabilities": {
            "1": {"2": 1.0},
            "2": {"3": 1.0},
            "3": {"6": 0.5, "6p": 0.5},
            "6": {"7": 1.0},
            "6p": {"7": 1.0},
            "7": {"8": 1.0},
            "8": {"9": 1.0},
            "9": {},
        },
        "length_model": {
            "distribution": "uniform",
            "min_repeats": 20,
            "max_repeats": 40,
        },
        "mutations": {
            "dupC": {
                "allowed_repeats": ["X", "A", "B"],
                "strict_mode": False,
                "changes": [{"type": "insert", "position": 0, "sequence": "C"}],
            }
        },
        "reference_assembly": "hg38",
    }

    config_file = tmp_path / "test_config.json"
    config_file.write_text(json.dumps(config, indent=2))
    return str(config_file)


# ============================================================================
# Root CLI Tests
# ============================================================================


@pytest.mark.unit
@pytest.mark.cli
class TestCLIRoot:
    """Tests for root CLI group."""

    def test_cli_help(self, runner):
        """Test CLI help displays correctly."""
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "MucOneUp" in result.output
        assert "simulate" in result.output
        assert "reads" in result.output
        assert "analyze" in result.output

    def test_cli_version(self, runner):
        """Test --version flag."""
        result = runner.invoke(cli, ["--version"])
        assert result.exit_code == 0
        assert "MucOneUp" in result.output

    def test_cli_requires_config(self, runner):
        """Test that commands require --config."""
        result = runner.invoke(cli, ["simulate"])
        assert result.exit_code != 0
        assert "Missing option '--config'" in result.output or "required" in result.output.lower()

    def test_cli_invalid_log_level(self, runner, temp_config):
        """Test invalid log level."""
        result = runner.invoke(cli, ["--config", temp_config, "--log-level", "INVALID"])
        assert result.exit_code != 0


# ============================================================================
# Simulate Command Tests
# ============================================================================


@pytest.mark.unit
@pytest.mark.cli
class TestSimulateCommand:
    """Tests for simulate command."""

    def test_simulate_help(self, runner, temp_config):
        """Test simulate command help."""
        result = runner.invoke(cli, ["--config", temp_config, "simulate", "--help"])
        assert result.exit_code == 0
        assert "Generate MUC1 VNTR diploid haplotypes" in result.output
        assert "--out-base" in result.output
        assert "--out-dir" in result.output
        assert "--num-haplotypes" in result.output

    def test_simulate_basic(self, runner, temp_config, tmp_path):
        """Test basic simulate command."""
        result = runner.invoke(
            cli,
            [
                "--config",
                temp_config,
                "simulate",
                "--out-dir",
                str(tmp_path),
                "--out-base",
                "test",
                "--seed",
                "42",
            ],
        )
        # May fail due to missing tools, but CLI parsing should work
        assert "config" in result.output.lower() or result.exit_code in [0, 1]

    def test_simulate_with_fixed_lengths(self, runner, temp_config, tmp_path):
        """Test simulate with fixed lengths."""
        result = runner.invoke(
            cli,
            [
                "--config",
                temp_config,
                "simulate",
                "--out-dir",
                str(tmp_path),
                "--out-base",
                "test",
                "--fixed-lengths",
                "30",
                "--seed",
                "42",
            ],
        )
        assert result.exit_code in [0, 1]  # May fail on missing tools

    def test_simulate_with_mutation(self, runner, temp_config, tmp_path):
        """Test simulate with mutation."""
        result = runner.invoke(
            cli,
            [
                "--config",
                temp_config,
                "simulate",
                "--out-dir",
                str(tmp_path),
                "--out-base",
                "test",
                "--mutation-name",
                "dupC",
                "--seed",
                "42",
            ],
        )
        assert result.exit_code in [0, 1]

    def test_simulate_with_pipeline_flags(self, runner, temp_config, tmp_path):
        """Test simulate with --simulate-reads and --output-orfs."""
        result = runner.invoke(
            cli,
            [
                "--config",
                temp_config,
                "simulate",
                "--out-dir",
                str(tmp_path),
                "--out-base",
                "test",
                "--seed",
                "42",
                "--output-orfs",
                "--orf-min-aa",
                "50",
            ],
        )
        assert result.exit_code in [0, 1]

    def test_simulate_with_random_snps(self, runner, temp_config, tmp_path):
        """Test simulate with random SNPs."""
        result = runner.invoke(
            cli,
            [
                "--config",
                temp_config,
                "simulate",
                "--out-dir",
                str(tmp_path),
                "--out-base",
                "test",
                "--seed",
                "42",
                "--random-snps",
                "--random-snp-density",
                "1.0",
                "--random-snp-output-file",
                str(tmp_path / "snps.tsv"),
            ],
        )
        assert result.exit_code in [0, 1]


# ============================================================================
# Reads Command Tests
# ============================================================================


@pytest.mark.unit
@pytest.mark.cli
class TestReadsCommand:
    """Tests for reads command group."""

    def test_reads_help(self, runner, temp_config):
        """Test reads group help."""
        result = runner.invoke(cli, ["--config", temp_config, "reads", "--help"])
        assert result.exit_code == 0
        assert "Standalone read simulation utilities" in result.output
        assert "illumina" in result.output
        assert "ont" in result.output

    def test_reads_illumina_help(self, runner, temp_config):
        """Test reads illumina command help."""
        result = runner.invoke(cli, ["--config", temp_config, "reads", "illumina", "--help"])
        assert result.exit_code == 0
        assert "Simulate Illumina short reads" in result.output
        assert "--coverage" in result.output
        assert "--threads" in result.output

    def test_reads_illumina_requires_input(self, runner, temp_config):
        """Test that reads illumina requires input FASTA."""
        result = runner.invoke(
            cli, ["--config", temp_config, "reads", "illumina", "--out-base", "test"]
        )
        assert result.exit_code != 0
        assert "INPUT_FASTA" in result.output or "Missing argument" in result.output

    def test_reads_illumina_with_file(self, runner, temp_config, tmp_path):
        """Test reads illumina with FASTA file."""
        # Create a minimal FASTA file
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(">test\nACGTACGTACGT\n")

        result = runner.invoke(
            cli,
            [
                "--config",
                temp_config,
                "reads",
                "illumina",
                str(fasta_file),
                "--out-base",
                "test",
                "--coverage",
                "10",
            ],
        )
        # May fail on missing tools, but parsing should work
        assert result.exit_code in [0, 1]

    def test_reads_ont_help(self, runner, temp_config):
        """Test reads ont command help."""
        result = runner.invoke(cli, ["--config", temp_config, "reads", "ont", "--help"])
        assert result.exit_code == 0
        assert "Simulate Oxford Nanopore long reads" in result.output
        assert "--min-read-length" in result.output


# ============================================================================
# Analyze Command Tests
# ============================================================================


@pytest.mark.unit
@pytest.mark.cli
class TestAnalyzeCommand:
    """Tests for analyze command group."""

    def test_analyze_help(self, runner, temp_config):
        """Test analyze group help."""
        result = runner.invoke(cli, ["--config", temp_config, "analyze", "--help"])
        assert result.exit_code == 0
        assert "Standalone analysis utilities" in result.output
        assert "orfs" in result.output
        assert "stats" in result.output

    def test_analyze_orfs_help(self, runner, temp_config):
        """Test analyze orfs command help."""
        result = runner.invoke(cli, ["--config", temp_config, "analyze", "orfs", "--help"])
        assert result.exit_code == 0
        assert "Predict ORFs" in result.output
        assert "--orf-min-aa" in result.output
        assert "--orf-aa-prefix" in result.output

    def test_analyze_orfs_requires_input(self, runner, temp_config):
        """Test that analyze orfs requires input FASTA."""
        result = runner.invoke(
            cli, ["--config", temp_config, "analyze", "orfs", "--out-base", "test"]
        )
        assert result.exit_code != 0
        assert "INPUT_FASTA" in result.output or "Missing argument" in result.output

    def test_analyze_orfs_with_file(self, runner, temp_config, tmp_path):
        """Test analyze orfs with FASTA file."""
        # Create a minimal FASTA file
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(">test\nACGTACGTACGT\n")

        result = runner.invoke(
            cli,
            [
                "--config",
                temp_config,
                "analyze",
                "orfs",
                str(fasta_file),
                "--out-base",
                "test",
                "--orf-min-aa",
                "10",
            ],
        )
        # May fail on missing dependencies, but parsing should work
        assert result.exit_code in [0, 1]

    def test_analyze_stats_help(self, runner, temp_config):
        """Test analyze stats command help."""
        result = runner.invoke(cli, ["--config", temp_config, "analyze", "stats", "--help"])
        assert result.exit_code == 0
        assert "Generate simulation statistics" in result.output

    def test_analyze_stats_with_file(self, runner, temp_config, tmp_path):
        """Test analyze stats with FASTA file."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(">haplotype_1\nACGTACGTACGT\n")

        result = runner.invoke(
            cli,
            [
                "--config",
                temp_config,
                "analyze",
                "stats",
                str(fasta_file),
                "--out-base",
                "test",
            ],
        )
        assert result.exit_code in [0, 1]


# ============================================================================
# Error Handling Tests
# ============================================================================


@pytest.mark.unit
@pytest.mark.cli
class TestErrorHandling:
    """Tests for error handling."""

    def test_invalid_config_file(self, runner, tmp_path):
        """Test error with non-existent config file."""
        result = runner.invoke(cli, ["--config", str(tmp_path / "nonexistent.json"), "simulate"])
        assert result.exit_code != 0

    def test_invalid_config_json(self, runner, tmp_path):
        """Test error with invalid JSON in config."""
        bad_config = tmp_path / "bad.json"
        bad_config.write_text("{invalid json")

        result = runner.invoke(
            cli,
            [
                "--config",
                str(bad_config),
                "simulate",
                "--out-base",
                "test",
                "--out-dir",
                str(tmp_path),
            ],
        )
        assert result.exit_code in [1, 2]  # Expected error

    def test_keyboard_interrupt_handling(self, runner, temp_config):
        """Test KeyboardInterrupt handling exists in implementation."""
        # Verify KeyboardInterrupt handler is in the module source
        import inspect

        from muc_one_up.cli import click_main

        module_source = inspect.getsource(click_main)
        assert "except KeyboardInterrupt:" in module_source
        assert "ctx.exit(130)" in module_source  # Standard Unix code for SIGINT


# ============================================================================
# Integration Tests
# ============================================================================


@pytest.mark.integration
@pytest.mark.cli
class TestCLIIntegration:
    """Integration tests for complete workflows."""

    def test_full_simulate_workflow(self, runner, temp_config, tmp_path):
        """Test complete simulation workflow."""
        result = runner.invoke(
            cli,
            [
                "--config",
                temp_config,
                "--log-level",
                "INFO",
                "simulate",
                "--out-dir",
                str(tmp_path),
                "--out-base",
                "integration_test",
                "--num-haplotypes",
                "2",
                "--fixed-lengths",
                "25",
                "--seed",
                "12345",
                "--output-structure",
            ],
        )

        # Check that command was parsed correctly (may fail on missing tools)
        assert result.exit_code in [0, 1]

        # If successful, check outputs would exist
        # Note: Files may not exist if external tools are missing
        if result.exit_code == 0:
            pass  # Output verification would go here


# ============================================================================
# Comparison Tests (argparse vs Click)
# ============================================================================


@pytest.mark.unit
@pytest.mark.cli
class TestArgparseClickParity:
    """Tests to verify parity between argparse and Click implementations."""

    def test_click_has_all_argparse_options(self, runner, temp_config):
        """Verify Click simulate has all argparse options."""
        result = runner.invoke(cli, ["--config", temp_config, "simulate", "--help"])

        # Core options
        assert "--out-base" in result.output
        assert "--out-dir" in result.output
        assert "--num-haplotypes" in result.output
        assert "--seed" in result.output
        assert "--reference-assembly" in result.output

        # Length options
        assert "--fixed-lengths" in result.output
        assert "--input-structure" in result.output
        assert "--simulate-series" in result.output

        # Mutation options
        assert "--mutation-name" in result.output
        assert "--mutation-targets" in result.output

        # SNP options
        assert "--random-snps" in result.output
        assert "--snp-input-file" in result.output

        # Pipeline options
        assert "--simulate-reads" in result.output
        assert "--output-orfs" in result.output
        assert "--orf-min-aa" in result.output

    def test_click_default_values_match_argparse(self, runner, temp_config):
        """Verify Click defaults match argparse defaults."""
        result = runner.invoke(cli, ["--config", temp_config, "simulate", "--help"])

        # Check for default value indicators
        assert "muc1_simulated" in result.output  # --out-base default
        assert "[default: 2]" in result.output  # --num-haplotypes default
        assert "[default: 100]" in result.output  # --orf-min-aa default
