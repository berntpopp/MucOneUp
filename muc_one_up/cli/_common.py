"""Shared CLI utilities for MucOneUp.

Contains constants, logging configuration, and Click callbacks
shared across CLI modules.
"""

import logging

import click

from ..version import __version__

# ============================================================================
# CLI Context Settings
# ============================================================================

# Context settings for custom help flags (-h in addition to --help)
CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


# ============================================================================
# Logging Configuration
# ============================================================================


def configure_logging(level_str: str) -> None:
    """Configure root logging based on the provided level string.

    If level_str is 'NONE', disable logging.
    """
    if level_str.upper() == "NONE":
        logging.disable(logging.CRITICAL)
    else:
        level = getattr(logging, level_str.upper(), logging.INFO)
        root_logger = logging.getLogger()
        if root_logger.handlers:
            for handler in root_logger.handlers[:]:
                root_logger.removeHandler(handler)
        logging.basicConfig(level=level, format="%(asctime)s - %(levelname)s - %(message)s")
        logging.info("Logging configured at level: %s", level_str.upper())


def require_config(ctx: click.Context) -> None:
    """Ensure config is provided before executing a command.

    Call this at the start of each command that needs config.
    Provides clear error message if config is missing.

    Args:
        ctx: Click context

    Raises:
        click.UsageError: If config is not provided

    Note:
        Includes defensive check for ctx.obj existence, even though
        cli() calls ensure_object(dict). This follows defensive
        programming best practices.
    """
    if not ctx.obj or not ctx.obj.get("config_path"):
        raise click.UsageError(
            "Missing required option '--config'.\n"
            "Usage: muconeup --config FILE COMMAND [OPTIONS]\n"
            "Try 'muconeup --help' for more information."
        )


def print_version_callback(ctx: click.Context, param: click.Parameter, value: bool) -> None:
    """Print version and exit when --version or -V is used.

    This callback is marked as eager to ensure it runs before
    required options are validated, allowing --version to work
    without providing --config.

    Args:
        ctx: Click context
        param: Click parameter
        value: Whether the flag was set
    """
    if not value or ctx.resilient_parsing:
        return
    click.echo(f"MucOneUp, version {__version__}")
    ctx.exit()
