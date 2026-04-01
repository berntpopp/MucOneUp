"""Smoke tests for CLI import and basic operability.

These tests verify that the CLI can be imported and basic commands
work without a full environment (no config file, no external tools).
"""

from click.testing import CliRunner


def test_cli_module_imports():
    """The CLI module should import without errors."""
    from muc_one_up.cli import main  # noqa: F401


def test_cli_help_exits_zero():
    """muconeup --help should exit 0 without requiring config."""
    from muc_one_up.cli.click_main import cli

    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "MucOneUp" in result.output


def test_cli_version_exits_zero():
    """muconeup --version should exit 0 without requiring config."""
    from muc_one_up.cli.click_main import cli
    from muc_one_up.version import __version__

    runner = CliRunner()
    result = runner.invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.output


def test_simulate_help_exits_zero():
    """muconeup simulate --help should exit 0 without requiring config."""
    from muc_one_up.cli.click_main import cli

    runner = CliRunner()
    result = runner.invoke(cli, ["simulate", "--help"])
    assert result.exit_code == 0


def test_reads_help_exits_zero():
    """muconeup reads --help should exit 0 without requiring config."""
    from muc_one_up.cli.click_main import cli

    runner = CliRunner()
    result = runner.invoke(cli, ["reads", "--help"])
    assert result.exit_code == 0


def test_analyze_help_exits_zero():
    """muconeup analyze --help should exit 0 without requiring config."""
    from muc_one_up.cli.click_main import cli

    runner = CliRunner()
    result = runner.invoke(cli, ["analyze", "--help"])
    assert result.exit_code == 0
