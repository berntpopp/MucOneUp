"""Structural tests for CLI module organization."""

import inspect


def test_click_main_under_100_lines():
    """click_main.py should be registration-only, under 100 lines."""
    import muc_one_up.cli.click_main as mod

    source = inspect.getsource(mod)
    line_count = len(source.splitlines())
    assert line_count < 100, (
        f"click_main.py is {line_count} lines. Should be under 100 (registration only)."
    )


def test_no_simple_namespace_in_click_main():
    """No SimpleNamespace usage in click_main.py."""
    import muc_one_up.cli.click_main as mod

    source = inspect.getsource(mod)
    assert "SimpleNamespace" not in source


def test_commands_registered():
    """All command groups should be registered on the root CLI."""
    from muc_one_up.cli.click_main import cli

    command_names = list(cli.commands.keys())
    assert "simulate" in command_names
    assert "reads" in command_names
    assert "analyze" in command_names
