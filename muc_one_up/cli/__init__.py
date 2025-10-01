"""
CLI package for MucOneUp.

Click-based CLI following Unix philosophy and SOLID principles.

Public API:
    - main(): Entry point for the Click CLI application
"""

from .click_main import main

__all__ = ["main"]
