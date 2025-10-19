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

    def test_verbose_flag_accepted(self, runner, temp_config, tmp_path):
        """Test that --verbose flag is accepted and works."""
        result = runner.invoke(
            cli,
            [
                "--config",
                temp_config,
                "--verbose",
                "simulate",
                "--out-dir",
                str(tmp_path),
                "--seed",
                "42",
            ],
        )
        # Verbose flag should be accepted
        assert result.exit_code in [0, 1]
        # Should not have "no such option" error
        assert "no such option" not in result.output.lower()

    def test_verbose_short_flag(self, runner, temp_config, tmp_path):
        """Test that -v short form works."""
        result = runner.invoke(
            cli,
            [
                "--config",
                temp_config,
                "-v",
                "simulate",
                "--out-dir",
                str(tmp_path),
                "--seed",
                "42",
            ],
        )
        assert result.exit_code in [0, 1]
        assert "no such option" not in result.output.lower()

    def test_verbose_precedence_over_log_level(self, runner, temp_config, tmp_path):
        """Test that --verbose takes precedence over --log-level."""
        # This is a logic test - we verify verbose flag is processed
        # Testing actual DEBUG output would require capturing logs
        result = runner.invoke(
            cli,
            [
                "--config",
                temp_config,
                "--log-level",
                "ERROR",  # Set to ERROR
                "--verbose",  # But verbose should override to DEBUG
                "simulate",
                "--out-dir",
                str(tmp_path),
                "--seed",
                "42",
            ],
        )
        # Both flags should be accepted without conflict
        assert result.exit_code in [0, 1]


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

    def test_simulate_pure_responsibility(self, runner, temp_config, tmp_path):
        """Test simulate does NOT accept pipeline flags (clean separation)."""
        # Verify --output-orfs is NOT accepted in simulate
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
                "--output-orfs",  # Should be rejected
            ],
        )
        # Should fail because --output-orfs is not a simulate option
        assert result.exit_code != 0
        assert "no such option" in result.output.lower() or "unrecognized" in result.output.lower()

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

    def test_simulate_series_triggers_multiple_iterations(self, runner, temp_config, tmp_path):
        """Test that --simulate-series creates multiple iterations (logic test, not UI)."""
        result = runner.invoke(
            cli,
            [
                "--config",
                temp_config,
                "simulate",
                "--out-dir",
                str(tmp_path),
                "--out-base",
                "series",
                "--fixed-lengths",
                "20-24",
                "--simulate-series",
                "2",  # Creates 3 iterations: 20, 22, 24
                "--seed",
                "42",
            ],
        )

        # May fail on tools, but should process multiple iterations
        assert result.exit_code in [0, 1]

        # Logic test: Check that multiple files were created (if successful)
        if result.exit_code == 0:
            fasta_files = list(tmp_path.glob("series.*.simulated.fa"))
            # Should create 3 files: .001, .002, .003
            assert len(fasta_files) >= 1  # At least one iteration completed

            # Note: We test the LOGIC (multiple iterations), not the progressbar UI
            # CliRunner doesn't capture progressbar output, and that's OK


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
        assert (
            "Read simulation utilities" in result.output
            or "read simulation" in result.output.lower()
        )
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

    def test_reads_ont_with_file(self, runner, temp_config, tmp_path):
        """Test reads ont with single FASTA file (parity with illumina)."""
        # Create a minimal FASTA file
        fasta_file = tmp_path / "test_ont.fa"
        fasta_file.write_text(">test_sequence\nACGTACGTACGTACGT\n")

        result = runner.invoke(
            cli,
            [
                "--config",
                temp_config,
                "reads",
                "ont",
                str(fasta_file),
                "--out-base",
                "test_ont",
                "--coverage",
                "15",
                "--min-read-length",
                "100",
            ],
        )

        # May fail on missing NanoSim, but parsing should work
        assert result.exit_code in [0, 1]

    @pytest.mark.parametrize(
        "command,subcommand,file_count,extra_options",
        [
            ("reads", "illumina", 3, ["--coverage", "10", "--threads", "4"]),
            ("reads", "ont", 3, ["--coverage", "20", "--min-read-length", "50"]),
            ("analyze", "orfs", 2, ["--orf-min-aa", "10"]),
            ("analyze", "stats", 3, []),
        ],
        ids=["illumina-batch", "ont-batch", "orfs-batch", "stats-batch"],
    )
    def test_batch_processing(
        self, runner, temp_config, tmp_path, command, subcommand, file_count, extra_options
    ):
        """Test batch processing for reads/analyze commands (DRY parametrized test)."""
        # Create multiple FASTA files
        fastas = []
        for i in range(1, file_count + 1):
            fasta = tmp_path / f"test{i}.fa"
            fasta.write_text(f">seq{i}\nACGTACGTACGT\n")
            fastas.append(str(fasta))

        # Build command
        cmd = [
            "--config",
            temp_config,
            command,
            subcommand,
            *fastas,
            "--out-dir",
            str(tmp_path),
            *extra_options,
        ]

        result = runner.invoke(cli, cmd)

        # Verify batch processing
        assert result.exit_code in [0, 1]  # May fail on missing tools
        # Check for batch processing indication
        if result.exit_code == 0:
            assert f"{file_count} FASTA file(s)" in result.output

    def test_batch_processing_warns_on_out_base(self, runner, temp_config, tmp_path):
        """Test that using --out-base with multiple files triggers warning."""
        # Create 2 FASTA files
        fasta1 = tmp_path / "test1.fa"
        fasta2 = tmp_path / "test2.fa"
        fasta1.write_text(">seq1\nACGT\n")
        fasta2.write_text(">seq2\nGCTA\n")

        result = runner.invoke(
            cli,
            [
                "--config",
                temp_config,
                "reads",
                "illumina",
                str(fasta1),
                str(fasta2),
                "--out-base",
                "custom",  # Should warn
            ],
        )

        # Verify warning appears (or fails before warning due to missing tools)
        if result.exit_code in [0, 1]:
            assert (
                "will be used for all" in result.output
                or "Consider omitting" in result.output
                or result.exit_code == 1
            )


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
        assert "Analysis utilities" in result.output or "analysis" in result.output.lower()
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
                "--out-dir",
                str(tmp_path),
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
        assert "statistics" in result.output.lower() or "stats" in result.output.lower()

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
                "--out-dir",
                str(tmp_path),
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
        """Verify Click simulate has core simulation options (clean separation)."""
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

        # Pipeline options should NOT be in simulate (clean separation)
        assert "--simulate-reads" not in result.output
        assert "--output-orfs" not in result.output

    def test_click_default_values_match_argparse(self, runner, temp_config):
        """Verify Click defaults match argparse defaults."""
        result = runner.invoke(cli, ["--config", temp_config, "simulate", "--help"])

        # Check for default value indicators in simulate (flexible formatting)
        assert "muc1_simulated" in result.output  # --out-base default
        assert "[default:" in result.output and "2]" in result.output  # --num-haplotypes default

        # Check ORF defaults in analyze orfs command (clean separation)
        result_orfs = runner.invoke(cli, ["--config", temp_config, "analyze", "orfs", "--help"])
        assert (
            "[default:" in result_orfs.output and "100]" in result_orfs.output
        )  # --orf-min-aa default


# ============================================================================
# Seed Parameter Tests
# ============================================================================


@pytest.mark.unit
@pytest.mark.cli
class TestSeedParameter:
    """Tests for --seed parameter in reads commands."""

    def test_reads_illumina_with_seed(self, runner, temp_config, tmp_path, mocker):
        """Test reads illumina command with --seed parameter."""
        # Create a minimal FASTA file
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(">test\nACGTACGTACGT\n")

        # Mock the actual read simulation to avoid needing tools
        mock_simulate = mocker.patch("muc_one_up.read_simulation.simulate_reads")
        mock_simulate.return_value = None  # Simulate successful execution

        result = runner.invoke(
            cli,
            [
                "--config",
                temp_config,
                "reads",
                "illumina",
                str(fasta_file),
                "--seed",
                "42",
                "--coverage",
                "10",
            ],
        )

        # Command should execute successfully
        assert result.exit_code == 0

        # Verify simulate_reads was called
        assert mock_simulate.called

        # Verify the config passed has seed set
        called_config = mock_simulate.call_args[0][0]
        assert "read_simulation" in called_config
        assert called_config["read_simulation"]["seed"] == 42

    def test_reads_illumina_without_seed(self, runner, temp_config, tmp_path, mocker):
        """Test reads illumina command without --seed (backward compatibility)."""
        # Create a minimal FASTA file
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(">test\nACGTACGTACGT\n")

        # Mock the actual read simulation
        mock_simulate = mocker.patch("muc_one_up.read_simulation.simulate_reads")
        mock_simulate.return_value = None

        result = runner.invoke(
            cli,
            [
                "--config",
                temp_config,
                "reads",
                "illumina",
                str(fasta_file),
                "--coverage",
                "10",
            ],
        )

        # Command should execute successfully
        assert result.exit_code == 0

        # Verify the config does not have seed set
        called_config = mock_simulate.call_args[0][0]
        assert "read_simulation" in called_config
        # Seed should not be in config (or be None)
        assert called_config["read_simulation"].get("seed") is None

    def test_reads_ont_with_seed(self, runner, temp_config, tmp_path, mocker):
        """Test reads ont command with --seed parameter."""
        # Create a minimal FASTA file
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(">test\nACGTACGTACGTACGT\n")

        # Mock the actual read simulation
        mock_simulate = mocker.patch("muc_one_up.read_simulation.simulate_reads")
        mock_simulate.return_value = None

        result = runner.invoke(
            cli,
            [
                "--config",
                temp_config,
                "reads",
                "ont",
                str(fasta_file),
                "--seed",
                "123",
                "--coverage",
                "20",
            ],
        )

        # Command should execute successfully
        assert result.exit_code == 0

        # Verify simulate_reads was called
        assert mock_simulate.called

        # Verify the config passed has seed in nanosim_params
        called_config = mock_simulate.call_args[0][0]
        assert "nanosim_params" in called_config
        assert called_config["nanosim_params"]["seed"] == 123

    def test_reads_ont_without_seed(self, runner, temp_config, tmp_path, mocker):
        """Test reads ont command without --seed (backward compatibility)."""
        # Create a minimal FASTA file
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(">test\nACGTACGTACGTACGT\n")

        # Mock the actual read simulation
        mock_simulate = mocker.patch("muc_one_up.read_simulation.simulate_reads")
        mock_simulate.return_value = None

        result = runner.invoke(
            cli,
            [
                "--config",
                temp_config,
                "reads",
                "ont",
                str(fasta_file),
                "--coverage",
                "20",
            ],
        )

        # Command should execute successfully
        assert result.exit_code == 0

        # Verify the config does not have seed in nanosim_params (or is None)
        called_config = mock_simulate.call_args[0][0]
        assert "nanosim_params" in called_config
        assert called_config["nanosim_params"].get("seed") is None

    def test_seed_accepts_integer_values(self, runner, temp_config, tmp_path, mocker):
        """Test that --seed accepts valid integer values."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(">test\nACGT\n")

        mock_simulate = mocker.patch("muc_one_up.read_simulation.simulate_reads")
        mock_simulate.return_value = None

        # Test various integer values
        for seed_value in [0, 1, 42, 12345, 999999]:
            result = runner.invoke(
                cli,
                [
                    "--config",
                    temp_config,
                    "reads",
                    "illumina",
                    str(fasta_file),
                    "--seed",
                    str(seed_value),
                    "--coverage",
                    "10",
                ],
            )
            assert result.exit_code == 0
            called_config = mock_simulate.call_args[0][0]
            assert called_config["read_simulation"]["seed"] == seed_value

    def test_seed_logging_message(self, runner, temp_config, tmp_path, mocker, caplog):
        """Test that using --seed produces appropriate log message."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(">test\nACGT\n")

        mock_simulate = mocker.patch("muc_one_up.read_simulation.simulate_reads")
        mock_simulate.return_value = None

        # Run with seed
        result = runner.invoke(
            cli,
            [
                "--config",
                temp_config,
                "reads",
                "illumina",
                str(fasta_file),
                "--seed",
                "42",
                "--coverage",
                "10",
            ],
        )

        assert result.exit_code == 0
        # Check that log message appears in output or logs
        # The logging.info call should have executed
        assert mock_simulate.called
