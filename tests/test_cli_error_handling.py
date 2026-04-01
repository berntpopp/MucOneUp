"""Tests for CLI error handling decorator."""

import click
from click.testing import CliRunner

from muc_one_up.cli.error_handling import cli_error_handler
from muc_one_up.exceptions import (
    ConfigurationError,
    MucOneUpError,
    SimulationError,
)


def _make_test_command(exception_to_raise):
    @click.command()
    @cli_error_handler
    def cmd():
        raise exception_to_raise

    return cmd


class TestCLIErrorHandler:
    def test_keyboard_interrupt_exits_130(self):
        runner = CliRunner()
        result = runner.invoke(_make_test_command(KeyboardInterrupt()))
        assert result.exit_code == 130

    def test_muc_one_up_error_exits_1(self):
        runner = CliRunner()
        result = runner.invoke(_make_test_command(MucOneUpError("test")))
        assert result.exit_code == 1

    def test_simulation_error_exits_1(self):
        runner = CliRunner()
        result = runner.invoke(_make_test_command(SimulationError("fail")))
        assert result.exit_code == 1

    def test_configuration_error_exits_1(self):
        runner = CliRunner()
        result = runner.invoke(_make_test_command(ConfigurationError("bad")))
        assert result.exit_code == 1

    def test_generic_exception_exits_2(self):
        runner = CliRunner()
        result = runner.invoke(_make_test_command(RuntimeError("unexpected")))
        assert result.exit_code == 2

    def test_click_exception_passes_through(self):
        runner = CliRunner()
        result = runner.invoke(_make_test_command(click.ClickException("bad param")))
        assert result.exit_code == 1
        assert "bad param" in result.output

    def test_usage_error_passes_through(self):
        runner = CliRunner()
        result = runner.invoke(_make_test_command(click.UsageError("missing option")))
        assert result.exit_code == 2
        assert "missing option" in result.output

    def test_no_exception_exits_0(self):
        @click.command()
        @cli_error_handler
        def cmd():
            click.echo("ok")

        runner = CliRunner()
        result = runner.invoke(cmd)
        assert result.exit_code == 0
        assert "ok" in result.output
