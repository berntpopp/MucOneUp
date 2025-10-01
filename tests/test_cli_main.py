"""Tests for muc_one_up.cli.main module.

Tests cover:
- Argument parser creation and validation
- Logging configuration
- Main entry point orchestration
- Error handling and exit codes
"""

import logging

import pytest

from muc_one_up.cli.main import build_parser, configure_logging


@pytest.mark.unit
class TestBuildParser:
    """Tests for build_parser function."""

    def test_build_parser_creates_parser(self):
        """Test that build_parser returns an ArgumentParser."""
        parser = build_parser()
        assert parser is not None
        assert hasattr(parser, "parse_args")

    def test_parser_has_required_config_argument(self):
        """Test that config argument is required."""
        parser = build_parser()

        # Should fail without --config
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_parser_config_argument(self):
        """Test parsing --config argument."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json"])
        assert args.config == "config.json"

    def test_parser_out_base_default(self):
        """Test that out-base has correct default."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json"])
        assert args.out_base == "muc1_simulated"

    def test_parser_out_base_custom(self):
        """Test custom out-base."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--out-base", "custom"])
        assert args.out_base == "custom"

    def test_parser_out_dir_default(self):
        """Test that out-dir has correct default."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json"])
        assert args.out_dir == "."

    def test_parser_out_dir_custom(self):
        """Test custom out-dir."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--out-dir", "/tmp/output"])
        assert args.out_dir == "/tmp/output"

    def test_parser_num_haplotypes_default(self):
        """Test that num-haplotypes defaults to 2."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json"])
        assert args.num_haplotypes == 2

    def test_parser_num_haplotypes_custom(self):
        """Test custom num-haplotypes."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--num-haplotypes", "4"])
        assert args.num_haplotypes == 4

    def test_parser_fixed_lengths_single(self):
        """Test --fixed-lengths with single value."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--fixed-lengths", "60"])
        assert args.fixed_lengths == ["60"]

    def test_parser_fixed_lengths_multiple(self):
        """Test --fixed-lengths with multiple values."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--fixed-lengths", "20-40", "30-50"])
        assert args.fixed_lengths == ["20-40", "30-50"]

    def test_parser_input_structure(self):
        """Test --input-structure argument."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--input-structure", "structure.txt"])
        assert args.input_structure == "structure.txt"

    def test_parser_seed(self):
        """Test --seed argument."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--seed", "42"])
        assert args.seed == 42

    def test_parser_simulate_series_default(self):
        """Test --simulate-series with default step."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--simulate-series"])
        assert args.simulate_series == 1

    def test_parser_simulate_series_custom_step(self):
        """Test --simulate-series with custom step."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--simulate-series", "5"])
        assert args.simulate_series == 5

    def test_parser_mutation_name(self):
        """Test --mutation-name argument."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--mutation-name", "dupC"])
        assert args.mutation_name == "dupC"

    def test_parser_mutation_name_dual(self):
        """Test --mutation-name with dual simulation."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--mutation-name", "normal,dupC"])
        assert args.mutation_name == "normal,dupC"

    def test_parser_mutation_targets_single(self):
        """Test --mutation-targets with single target."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--mutation-targets", "1,25"])
        assert args.mutation_targets == ["1,25"]

    def test_parser_mutation_targets_multiple(self):
        """Test --mutation-targets with multiple targets."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--mutation-targets", "1,25", "2,30"])
        assert args.mutation_targets == ["1,25", "2,30"]

    def test_parser_output_structure_flag(self):
        """Test --output-structure flag."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--output-structure"])
        assert args.output_structure is True

    def test_parser_output_structure_default(self):
        """Test --output-structure defaults to False."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json"])
        assert args.output_structure is False

    def test_parser_output_orfs_flag(self):
        """Test --output-orfs flag."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--output-orfs"])
        assert args.output_orfs is True

    def test_parser_orf_min_aa_default(self):
        """Test --orf-min-aa default value."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json"])
        assert args.orf_min_aa == 100

    def test_parser_orf_min_aa_custom(self):
        """Test custom --orf-min-aa."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--orf-min-aa", "50"])
        assert args.orf_min_aa == 50

    def test_parser_orf_aa_prefix_default(self):
        """Test --orf-aa-prefix with default value."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--orf-aa-prefix"])
        assert args.orf_aa_prefix == "MTSSV"

    def test_parser_orf_aa_prefix_custom(self):
        """Test custom --orf-aa-prefix."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--orf-aa-prefix", "MTSSS"])
        assert args.orf_aa_prefix == "MTSSS"

    def test_parser_simulate_reads_default(self):
        """Test --simulate-reads with default value."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--simulate-reads"])
        assert args.simulate_reads == "illumina"

    def test_parser_simulate_reads_illumina(self):
        """Test --simulate-reads illumina."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--simulate-reads", "illumina"])
        assert args.simulate_reads == "illumina"

    def test_parser_simulate_reads_ont(self):
        """Test --simulate-reads ont."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--simulate-reads", "ont"])
        assert args.simulate_reads == "ont"

    def test_parser_log_level_default(self):
        """Test --log-level default."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json"])
        assert args.log_level == "INFO"

    def test_parser_log_level_debug(self):
        """Test --log-level DEBUG."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--log-level", "DEBUG"])
        assert args.log_level == "DEBUG"

    def test_parser_log_level_none(self):
        """Test --log-level NONE."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--log-level", "NONE"])
        assert args.log_level == "NONE"

    def test_parser_reference_assembly_hg38(self):
        """Test --reference-assembly hg38."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--reference-assembly", "hg38"])
        assert args.reference_assembly == "hg38"

    def test_parser_reference_assembly_hg19(self):
        """Test --reference-assembly hg19."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--reference-assembly", "hg19"])
        assert args.reference_assembly == "hg19"

    def test_parser_snp_input_file(self):
        """Test --snp-input-file argument."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--snp-input-file", "snps.tsv"])
        assert args.snp_input_file == "snps.tsv"

    def test_parser_random_snps_flag(self):
        """Test --random-snps flag."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--random-snps"])
        assert args.random_snps is True

    def test_parser_random_snp_density(self):
        """Test --random-snp-density argument."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--random-snp-density", "1.5"])
        assert args.random_snp_density == 1.5

    def test_parser_random_snp_output_file(self):
        """Test --random-snp-output-file argument."""
        parser = build_parser()
        args = parser.parse_args(
            ["--config", "config.json", "--random-snp-output-file", "snps_out.tsv"]
        )
        assert args.random_snp_output_file == "snps_out.tsv"

    def test_parser_random_snp_region_default(self):
        """Test --random-snp-region default value."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json"])
        assert args.random_snp_region == "constants_only"

    def test_parser_random_snp_region_all(self):
        """Test --random-snp-region all."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--random-snp-region", "all"])
        assert args.random_snp_region == "all"

    def test_parser_random_snp_haplotypes_default(self):
        """Test --random-snp-haplotypes default value."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json"])
        assert args.random_snp_haplotypes == "all"

    def test_parser_random_snp_haplotypes_1(self):
        """Test --random-snp-haplotypes 1."""
        parser = build_parser()
        args = parser.parse_args(["--config", "config.json", "--random-snp-haplotypes", "1"])
        assert args.random_snp_haplotypes == "1"

    def test_parser_version(self, capsys):
        """Test --version flag."""
        parser = build_parser()

        with pytest.raises(SystemExit) as exc_info:
            parser.parse_args(["--version"])

        assert exc_info.value.code == 0


@pytest.mark.unit
class TestConfigureLogging:
    """Tests for configure_logging function."""

    def test_configure_logging_info(self):
        """Test configuring logging at INFO level."""
        configure_logging("INFO")

        root_logger = logging.getLogger()
        assert root_logger.level == logging.INFO

    def test_configure_logging_debug(self):
        """Test configuring logging at DEBUG level."""
        configure_logging("DEBUG")

        root_logger = logging.getLogger()
        assert root_logger.level == logging.DEBUG

    def test_configure_logging_warning(self):
        """Test configuring logging at WARNING level."""
        configure_logging("WARNING")

        root_logger = logging.getLogger()
        assert root_logger.level == logging.WARNING

    def test_configure_logging_error(self):
        """Test configuring logging at ERROR level."""
        configure_logging("ERROR")

        root_logger = logging.getLogger()
        assert root_logger.level == logging.ERROR

    def test_configure_logging_critical(self):
        """Test configuring logging at CRITICAL level."""
        configure_logging("CRITICAL")

        root_logger = logging.getLogger()
        assert root_logger.level == logging.CRITICAL

    def test_configure_logging_none_disables_logging(self):
        """Test that NONE disables logging."""
        configure_logging("NONE")

        # Logging should be disabled
        assert logging.root.manager.disable == logging.CRITICAL

    def test_configure_logging_case_insensitive(self):
        """Test that logging configuration is case-insensitive."""
        configure_logging("info")

        root_logger = logging.getLogger()
        assert root_logger.level == logging.INFO


# Note: Integration tests for main() are complex due to heavy dependencies
# The unit tests above provide good coverage of argument parsing and logging
# Full integration testing is better suited for end-to-end tests
