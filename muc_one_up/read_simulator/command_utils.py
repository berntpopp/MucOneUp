#!/usr/bin/env python3
"""
Command building utilities for external tool wrappers.

SECURITY: This module provides safe command construction for subprocess execution
without shell=True. All command building MUST use these utilities to prevent
command injection vulnerabilities (CWE-78, Bandit B602).

This module centralizes the logic for handling both simple tool commands and
complex multi-word commands (e.g., conda/mamba environment activation).
"""

import shlex
from typing import Any


def build_tool_command(tool_cmd: str, *args: Any) -> list[str]:
    """
    Build a command list for subprocess execution from a tool command and arguments.

    This function safely handles multi-word tool commands (conda/mamba environments)
    by parsing them with shlex.split() when necessary. It prevents command injection
    vulnerabilities by constructing proper command lists for subprocess with shell=False.

    SECURITY RATIONALE:
    - subprocess.Popen/run with shell=False requires commands as list[str]
    - Multi-word commands from config (e.g., "mamba run -n env tool") must be split
    - shlex.split() safely parses quoted strings, escapes, and whitespace
    - Using this function prevents B602 HIGH severity command injection vulnerabilities
    - Centralized implementation reduces security regression risk

    Examples:
        >>> # Simple command
        >>> build_tool_command("bwa", "mem", "-t", "4", "ref.fa", "reads.fq")
        ["bwa", "mem", "-t", "4", "ref.fa", "reads.fq"]

        >>> # Conda/mamba command (multi-word)
        >>> build_tool_command("mamba run -n env_wessim reseq", "replaceN", "-r", "in.fa")
        ["mamba", "run", "-n", "env_wessim", "reseq", "replaceN", "-r", "in.fa"]

        >>> # Handles quoted strings correctly
        >>> build_tool_command('tool --arg "value with spaces"', "input.txt")
        ["tool", "--arg", "value with spaces", "input.txt"]

        >>> # Numeric arguments auto-converted
        >>> build_tool_command("samtools", "view", "-@", 4, "-b", "input.bam")
        ["samtools", "view", "-@", "4", "-b", "input.bam"]

    Args:
        tool_cmd: Tool command string from config (may contain spaces for conda/mamba).
                  Examples: "bwa", "mamba run --no-capture-output -n env_wessim reseq"
        *args: Additional command arguments (will be converted to strings).
               Can be strings, integers, floats, Path objects, etc.

    Returns:
        List of command arguments suitable for subprocess.Popen/run with shell=False.

    Raises:
        ValueError: If tool_cmd is empty or None.

    References:
        - CWE-78: Improper Neutralization of Special Elements (Command Injection)
        - Bandit B602: subprocess with shell=True (HIGH severity)
        - Python subprocess security: https://docs.python.org/3/library/subprocess.html#security
        - shlex documentation: https://docs.python.org/3/library/shlex.html

    See Also:
        get_tool_executable: Extract executable name from tool command
    """
    if not tool_cmd:
        raise ValueError("tool_cmd cannot be empty")

    # Check if command contains spaces (likely conda/mamba or complex command)
    # Use shlex.split() to safely parse complex commands with quoted arguments, escapes, and whitespace
    # e.g., "mamba run -n env tool" -> ["mamba", "run", "-n", "env", "tool"]
    # Simple commands (no spaces) are wrapped in a list: "bwa" -> ["bwa"]
    cmd_list = shlex.split(tool_cmd) if " " in tool_cmd else [tool_cmd]

    # Add additional arguments (convert all to strings)
    # This allows numeric arguments like threads=4 to be passed directly
    cmd_list.extend(str(arg) for arg in args)

    return cmd_list


def get_tool_executable(tool_cmd: str) -> str:
    """
    Extract the actual executable name from a tool command.

    Useful for logging, error messages, and debugging. This function identifies
    the actual tool being invoked, even when wrapped in conda/mamba commands.

    Examples:
        >>> get_tool_executable("bwa")
        "bwa"

        >>> get_tool_executable("mamba run -n env_wessim reseq")
        "reseq"

        >>> get_tool_executable("mamba run --no-capture-output -n env_wessim samtools")
        "samtools"

        >>> get_tool_executable("")
        "unknown"

    Args:
        tool_cmd: Tool command string (may contain conda/mamba prefix).

    Returns:
        The final executable name. For multi-word commands, returns the last
        non-flag element (typically the actual tool name). Returns "unknown"
        if the command is empty or cannot be determined.

    See Also:
        build_tool_command: Build full command list for subprocess execution
    """
    if not tool_cmd:
        return "unknown"

    if " " in tool_cmd:
        # Multi-word command - parse and find the tool name
        parts = shlex.split(tool_cmd)

        # The tool name is usually the last non-flag element
        # e.g., "mamba run -n env_wessim reseq" -> "reseq"
        # Walk backwards to find first element that doesn't start with "-"
        for part in reversed(parts):
            if not part.startswith("-"):
                return part

        # Fallback: if all elements start with "-", return the last one
        return parts[-1] if parts else "unknown"

    # Simple command - just return it
    return tool_cmd
