"""Centralized CLI error handling.

Maps domain exceptions to consistent exit codes:
- KeyboardInterrupt -> 130
- MucOneUpError (and subclasses) -> 1
- Unexpected exceptions -> 2
"""

from __future__ import annotations

import functools
import logging

import click
from click.exceptions import Exit as ClickExit

from ..exceptions import MucOneUpError


def cli_error_handler(func):
    """Decorator that catches exceptions and maps them to CLI exit codes."""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        ctx = click.get_current_context()
        try:
            return func(*args, **kwargs)
        except ClickExit:
            raise
        except click.ClickException:
            raise
        except KeyboardInterrupt:
            logging.info("Operation cancelled by user")
            ctx.exit(130)
        except MucOneUpError as e:
            logging.error("%s: %s", type(e).__name__, e)
            ctx.exit(1)
        except Exception as e:
            logging.error("Unexpected error: %s", e, exc_info=True)
            ctx.exit(2)

    return wrapper
