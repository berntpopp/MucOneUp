"""Tests for help and version access without --config requirement.

This module tests the improvements to CLI UX where help and version
flags are accessible without providing the --config option, while
still enforcing config requirement for actual command execution.
"""

import pytest
from click.testing import CliRunner

from muc_one_up.cli.click_main import cli
from muc_one_up.version import __version__


class TestRootHelpAccess:
    """Test root-level help access patterns."""

    def test_root_help_long_flag(self):
        """muconeup --help should work without config."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "MucOneUp - MUC1 VNTR diploid reference simulator" in result.output
        assert "Philosophy: Each command does ONE thing" in result.output

    def test_root_help_short_flag(self):
        """muconeup -h should work without config."""
        runner = CliRunner()
        result = runner.invoke(cli, ["-h"])
        assert result.exit_code == 0
        assert "MucOneUp - MUC1 VNTR diploid reference simulator" in result.output
        assert "Philosophy: Each command does ONE thing" in result.output


class TestVersionAccess:
    """Test version flag access patterns."""

    def test_version_long_flag(self):
        """muconeup --version should work without config."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--version"])
        assert result.exit_code == 0
        assert "MucOneUp, version" in result.output
        assert __version__ in result.output

    def test_version_short_flag(self):
        """muconeup -V should work without config."""
        runner = CliRunner()
        result = runner.invoke(cli, ["-V"])
        assert result.exit_code == 0
        assert "MucOneUp, version" in result.output
        assert __version__ in result.output


class TestVerboseFlag:
    """Test verbose flag still works correctly."""

    def test_verbose_long_flag_with_help(self):
        """muconeup --verbose --help should work."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--verbose", "--help"])
        assert result.exit_code == 0
        assert "MucOneUp - MUC1 VNTR diploid reference simulator" in result.output

    def test_verbose_short_flag_with_help(self):
        """muconeup -v --help should work (-v is verbose, not version)."""
        runner = CliRunner()
        result = runner.invoke(cli, ["-v", "--help"])
        assert result.exit_code == 0
        assert "MucOneUp - MUC1 VNTR diploid reference simulator" in result.output


class TestSubcommandGroupHelp:
    """Test subcommand group help access without config."""

    def test_simulate_help_without_config(self):
        """muconeup simulate --help should work without config."""
        runner = CliRunner()
        result = runner.invoke(cli, ["simulate", "--help"])
        assert result.exit_code == 0
        assert "Generate MUC1 VNTR diploid haplotypes" in result.output
        assert "--out-base" in result.output

    def test_reads_group_help_without_config(self):
        """muconeup reads --help should work without config."""
        runner = CliRunner()
        result = runner.invoke(cli, ["reads", "--help"])
        assert result.exit_code == 0
        assert "Read simulation utilities" in result.output

    def test_reads_illumina_help_without_config(self):
        """muconeup reads illumina --help should work without config."""
        runner = CliRunner()
        result = runner.invoke(cli, ["reads", "illumina", "--help"])
        assert result.exit_code == 0
        assert "Simulate Illumina short reads" in result.output

    def test_reads_ont_help_without_config(self):
        """muconeup reads ont --help should work without config."""
        runner = CliRunner()
        result = runner.invoke(cli, ["reads", "ont", "--help"])
        assert result.exit_code == 0
        assert "Simulate Oxford Nanopore long reads" in result.output

    def test_reads_pacbio_help_without_config(self):
        """muconeup reads pacbio --help should work without config."""
        runner = CliRunner()
        result = runner.invoke(cli, ["reads", "pacbio", "--help"])
        assert result.exit_code == 0
        assert "Simulate PacBio HiFi reads" in result.output

    def test_analyze_group_help_without_config(self):
        """muconeup analyze --help should work without config."""
        runner = CliRunner()
        result = runner.invoke(cli, ["analyze", "--help"])
        assert result.exit_code == 0
        assert "Analysis utilities" in result.output

    def test_analyze_orfs_help_without_config(self):
        """muconeup analyze orfs --help should work without config."""
        runner = CliRunner()
        result = runner.invoke(cli, ["analyze", "orfs", "--help"])
        assert result.exit_code == 0
        assert "Predict ORFs" in result.output

    def test_analyze_stats_help_without_config(self):
        """muconeup analyze stats --help should work without config."""
        runner = CliRunner()
        result = runner.invoke(cli, ["analyze", "stats", "--help"])
        assert result.exit_code == 0
        assert "Generate basic sequence statistics" in result.output

    def test_analyze_vntr_stats_help_without_config(self):
        """muconeup analyze vntr-stats --help should work without config."""
        runner = CliRunner()
        result = runner.invoke(cli, ["analyze", "vntr-stats", "--help"])
        assert result.exit_code == 0
        assert "Analyze VNTR structures" in result.output

    def test_analyze_snapshot_validate_help_without_config(self):
        """muconeup analyze snapshot-validate --help should work without config."""
        runner = CliRunner()
        result = runner.invoke(cli, ["analyze", "snapshot-validate", "--help"])
        assert result.exit_code == 0
        assert "Validate SNaPshot assay" in result.output


class TestConfigRequirementEnforcement:
    """Test that commands still require config when executing."""

    def test_simulate_requires_config(self):
        """muconeup simulate should fail without config."""
        runner = CliRunner()
        result = runner.invoke(cli, ["simulate", "--out-base", "test"])
        assert result.exit_code != 0
        assert "Missing required option '--config'" in result.output
        assert "muconeup --config FILE COMMAND" in result.output

    @pytest.mark.parametrize(
        "command_args,file_ext,file_content,test_description",
        [
            (
                ["reads", "illumina"],
                "fa",
                ">test\nATGC\n",
                "reads illumina",
            ),
            (
                ["reads", "ont"],
                "fa",
                ">test\nATGC\n",
                "reads ont",
            ),
            (
                ["reads", "pacbio"],
                "fa",
                ">test\nATGC\n",
                "reads pacbio",
            ),
            (
                ["analyze", "orfs"],
                "fa",
                ">test\nATGC\n",
                "analyze orfs",
            ),
            (
                ["analyze", "vntr-stats"],
                "tsv",
                "vntr\n123\n",
                "analyze vntr-stats",
            ),
            (
                ["analyze", "snapshot-validate", "--mutation", "dupC"],
                "fa",
                ">test\nATGC\n",
                "analyze snapshot-validate",
            ),
        ],
        ids=lambda x: x if isinstance(x, str) and " " in x else None,
    )
    def test_commands_requiring_input_files_need_config(
        self, tmp_path, command_args, file_ext, file_content, test_description
    ):
        """Test that various commands requiring input files also require --config.

        This parametrized test reduces code duplication while maintaining clarity
        about which commands are being tested. Each test case is clearly identified
        in test output via the test_description parameter.
        """
        # Create dummy input file with appropriate content and extension
        test_file = tmp_path / f"test.{file_ext}"
        test_file.write_text(file_content)

        # Build command args with file path in correct position
        # For snapshot-validate: command file --mutation dupC
        # For others: command file
        if "snapshot-validate" in command_args:
            # File path goes after 'analyze snapshot-validate' but before '--mutation'
            mutation_idx = command_args.index("--mutation")
            args = [*command_args[:mutation_idx], str(test_file), *command_args[mutation_idx:]]
        else:
            args = [*command_args, str(test_file)]

        runner = CliRunner()
        result = runner.invoke(cli, args)

        assert result.exit_code != 0, f"{test_description} should fail without --config"
        assert "Missing required option '--config'" in result.output


class TestErrorMessages:
    """Test that error messages are clear and actionable."""

    def test_config_error_message_format(self):
        """Config error should provide clear usage information."""
        runner = CliRunner()
        result = runner.invoke(cli, ["simulate", "--out-base", "test"])
        assert result.exit_code != 0
        assert "Missing required option '--config'" in result.output
        assert "muconeup --config FILE COMMAND" in result.output
        assert "muconeup --help" in result.output

    def test_help_hint_in_error(self):
        """Error message should suggest using --help."""
        runner = CliRunner()
        result = runner.invoke(cli, ["simulate"])
        assert "Try 'muconeup --help'" in result.output or "muconeup --help" in result.output


class TestBackwardCompatibility:
    """Test that existing command patterns still work."""

    def test_config_before_command_pattern(self, tmp_path):
        """Standard pattern: muconeup --config FILE COMMAND should work."""
        # Create a minimal config file
        config_file = tmp_path / "config.json"
        config_file.write_text('{"repeats": {}, "constants": {}}')

        runner = CliRunner()
        # This should NOT fail with config error (may fail for other reasons)
        result = runner.invoke(cli, ["--config", str(config_file), "simulate", "--help"])
        assert result.exit_code == 0
        assert "Generate MUC1 VNTR diploid haplotypes" in result.output

    def test_verbose_with_config_pattern(self, tmp_path):
        """Pattern: muconeup -v --config FILE should work."""
        config_file = tmp_path / "config.json"
        config_file.write_text('{"repeats": {}, "constants": {}}')

        runner = CliRunner()
        result = runner.invoke(cli, ["-v", "--config", str(config_file), "--help"])
        assert result.exit_code == 0
        assert "MucOneUp - MUC1 VNTR diploid reference simulator" in result.output


class TestShortFlagUniqueness:
    """Test that short flags don't conflict."""

    def test_h_is_help_not_other(self):
        """-h should be help, not any other option."""
        runner = CliRunner()
        result = runner.invoke(cli, ["-h"])
        assert result.exit_code == 0
        assert "MucOneUp - MUC1 VNTR diploid reference simulator" in result.output

    def test_v_is_verbose_not_version(self):
        """-v should be verbose, not version."""
        runner = CliRunner()
        result = runner.invoke(cli, ["-v", "--version"])
        # -v sets verbose, --version shows version
        assert result.exit_code == 0
        assert "version" in result.output.lower()

    def test_V_is_version(self):
        """-V (capital) should be version."""
        runner = CliRunner()
        result = runner.invoke(cli, ["-V"])
        assert result.exit_code == 0
        assert "version" in result.output.lower()
        assert __version__ in result.output
