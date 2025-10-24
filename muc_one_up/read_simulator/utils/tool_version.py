#!/usr/bin/env python3
"""Tool version detection for read simulation pipelines.

This module provides functionality to detect and capture versions of external
bioinformatics tools used in the simulation pipelines. Following KISS, DRY, and
SOLID principles, it provides simple string-based version information with
graceful degradation on failure.

Design Principles:
    - KISS: Simple string return ("samtools 1.17" or "N/A")
    - DRY: Reuses build_tool_command() for conda/mamba support
    - SOLID: Tool-specific parsers via strategy pattern (tool_map)
    - Never crashes: Graceful degradation with logging

Based on VariantCentrifuge battle-tested pattern and verified through real-world
testing in conda environments (see TESTING_RESULTS.md).
"""

import logging
import re
import subprocess
from collections.abc import Callable

from ..command_utils import build_tool_command

logger = logging.getLogger(__name__)


def get_tool_version(tool_cmd: str, tool_name: str) -> str:
    """
    Get version string for a bioinformatics tool.

    Args:
        tool_cmd: Command from config (e.g., "samtools" or "mamba run -n env samtools")
        tool_name: Canonical name for lookup ("samtools", "bwa", etc.)

    Returns:
        Version string like "samtools 1.17" or "N/A" on failure

    Examples:
        >>> get_tool_version("samtools", "samtools")
        "samtools 1.17"

        >>> get_tool_version("mamba run -n env_wessim reseq", "reseq")
        "ReSeq version 1.1"

        >>> get_tool_version("unknown_tool", "unknown")
        "N/A"

    Design:
        - KISS: Simple string return, not complex dict
        - DRY: Reuses build_tool_command() for conda/mamba
        - SOLID: Tool-specific parsers via strategy pattern (tool_map)
        - Never crashes (graceful degradation to "N/A")

    Special Cases:
        - bwa: Returns exit code 1, parses stderr
        - pbsim3: NO version in output, returns "pbsim3 (installed)"
        - faToTwoBit: NO version info, returns "faToTwoBit (version unknown)"
        - pblat: Parses version from usage line ("v. 36x2")
        - NanoSim: Binary is "simulator.py", not "nanosim"
    """

    # Tool-specific parser functions
    def parse_samtools(stdout: str, stderr: str) -> str:
        """Parse: samtools 1.17\\nUsing htslib..."""
        for line in stdout.splitlines():
            if line.startswith("samtools"):
                return line.strip()
        return "N/A"

    def parse_minimap2(stdout: str, stderr: str) -> str:
        """Parse: 2.28-r1209 (prepend tool name)"""
        line = stdout.strip()
        if line:
            return f"minimap2 {line}"
        return "N/A"

    def parse_bwa(stdout: str, stderr: str) -> str:
        """Parse: Version: 0.7.18-r1243 from stderr (exit code 1)"""
        for line in stderr.splitlines():
            if "version:" in line.lower():
                match = re.search(r"version:\s*([\w\.-]+)", line, re.IGNORECASE)
                if match:
                    return f"bwa {match.group(1)}"
        return "N/A"

    def parse_nanosim(stdout: str, stderr: str) -> str:
        """Parse: NanoSim 3.2.2 from simulator.py --version

        Note: Tool binary is "simulator.py", not "nanosim"
        """
        for line in stdout.splitlines():
            line_clean = line.strip()
            if line_clean and "nanosim" in line_clean.lower():
                return line_clean  # Returns "NanoSim 3.2.2"
        return "N/A"

    def parse_pbsim3(stdout: str, stderr: str) -> str:
        """Parse pbsim3 version.

        CRITICAL: pbsim3 conda package provides "pbsim" binary with NO version in output!
        Fallback: Return generic string since version not retrievable from tool itself.
        """
        # Check if this is pbsim output (usage message)
        if "pbsim" in stdout.lower() or "USAGE: pbsim" in stdout:
            # pbsim3 doesn't output version - return generic identifier
            return "pbsim3 (installed)"
        return "N/A"

    def parse_ccs(stdout: str, stderr: str) -> str:
        """Parse: ccs 6.4.0 (commit v6.4.0)"""
        for line in stdout.splitlines():
            if line.strip().startswith("ccs"):
                return line.strip()  # Returns "ccs 6.4.0 (commit v6.4.0)"
        return "N/A"

    def parse_reseq(stdout: str, stderr: str) -> str:
        """Parse: ReSeq version 1.1"""
        for line in stdout.splitlines():
            if "reseq" in line.lower() and "version" in line.lower():
                return line.strip()  # Returns "ReSeq version 1.1"
        return "N/A"

    def parse_fatotwobit(stdout: str, stderr: str) -> str:
        """Parse faToTwoBit version.

        UCSC tool with NO version flag or version in usage output.
        Return generic string to indicate tool is present.
        """
        if "fatotwobit" in stdout.lower() or "fatotwobit" in stderr.lower():
            return "faToTwoBit (version unknown)"
        return "N/A"

    def parse_pblat(stdout: str, stderr: str) -> str:
        """Parse: pblat - BLAT with parallel supports v. 36x2"""
        for line in stdout.splitlines():
            if "pblat" in line.lower():
                # Extract version using regex: "v. 36x2"
                match = re.search(r"v\.\s*([\w]+)", line)
                if match:
                    return f"pblat v.{match.group(1)}"  # Returns "pblat v.36x2"
        return "N/A"

    # Tool command configuration (Strategy Pattern)
    # NOTE: Tool names must match config.json "tools" keys
    tool_map: dict[str, dict] = {
        "samtools": {"args": ["--version"], "parse": parse_samtools},
        "minimap2": {"args": ["--version"], "parse": parse_minimap2},
        "nanosim": {"args": ["--version"], "parse": parse_nanosim},  # Binary: simulator.py
        "pbsim3": {"args": [], "parse": parse_pbsim3},  # Binary: pbsim, NO version output!
        "ccs": {"args": ["--version"], "parse": parse_ccs},
        "reseq": {"args": ["--version"], "parse": parse_reseq},  # --version works, not just --help
        "bwa": {"args": [], "parse": parse_bwa},  # No flag, stderr output
        "faToTwoBit": {"args": [], "parse": parse_fatotwobit},  # No version info available
        "pblat": {"args": [], "parse": parse_pblat},  # No --version flag, parse usage
    }

    # Lookup tool configuration
    if tool_name not in tool_map:
        logger.warning(f"No version logic for {tool_name}. Returning 'N/A'.")
        return "N/A"

    config = tool_map[tool_name]
    parse_func: Callable[[str, str], str] = config["parse"]
    args: list[str] = config["args"]

    try:
        # Build command (handles conda/mamba from config)
        cmd = build_tool_command(tool_cmd, *args)
        logger.debug(f"Getting version: {' '.join(cmd)}")

        # Run with check=False (tools may return non-zero)
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,  # CRITICAL: Don't raise on non-zero exit (bwa, pbsim, faToTwoBit)
            timeout=5,
        )

        # Parse version from stdout/stderr
        version = parse_func(result.stdout, result.stderr)

        if version == "N/A":
            logger.warning(
                f"Could not parse version for {tool_name}. "
                f"stdout: {result.stdout[:100]}, stderr: {result.stderr[:100]}"
            )
        else:
            logger.debug(f"Got version for {tool_name}: {version}")

        return version

    except subprocess.TimeoutExpired:
        logger.warning(f"Timeout getting version for {tool_name} (5s)")
        return "N/A"
    except Exception as e:
        logger.warning(f"Failed to get version for {tool_name}: {e}")
        return "N/A"


def capture_tool_versions(tools_config: dict[str, str]) -> dict[str, str]:
    """
    Capture versions for all configured tools.

    Args:
        tools_config: Dictionary mapping tool names to commands
            Example: {"samtools": "samtools", "reseq": "mamba run -n env reseq"}

    Returns:
        Dictionary mapping tool names to version strings
            Example: {"samtools": "samtools 1.17", "reseq": "ReSeq version 1.1"}

    Examples:
        >>> tools = {"samtools": "samtools", "minimap2": "minimap2"}
        >>> versions = capture_tool_versions(tools)
        >>> versions["samtools"]
        "samtools 1.17"

    Note:
        This function captures versions for all tools in parallel-safe manner.
        Failed version captures return "N/A" without affecting other tools.
    """
    return {name: get_tool_version(cmd, name) for name, cmd in tools_config.items()}


def log_tool_versions(versions: dict[str, str]) -> None:
    """
    Log versions in formatted table.

    Args:
        versions: Dictionary mapping tool names to version strings

    Examples:
        >>> versions = {"samtools": "samtools 1.17", "bwa": "bwa 0.7.18-r1243"}
        >>> log_tool_versions(versions)
        # Logs:
        # Tool Versions:
        # ┌──────────────┬─────────────────┐
        # │ Tool         │ Version         │
        # ├──────────────┼─────────────────┤
        # │ bwa          │ bwa 0.7.18-r1243│
        # │ samtools     │ samtools 1.17   │
        # └──────────────┴─────────────────┘
    """
    if not versions:
        return

    # Calculate column widths
    max_tool = max(len(t) for t in versions)
    max_ver = max(len(v) for v in versions.values())
    tw, vw = max(max_tool, 12), max(max_ver, 15)

    # Print table
    logger.info("Tool Versions:")
    logger.info(f"┌{'─' * (tw + 2)}┬{'─' * (vw + 2)}┐")
    logger.info(f"│ {'Tool':<{tw}} │ {'Version':<{vw}} │")
    logger.info(f"├{'─' * (tw + 2)}┼{'─' * (vw + 2)}┤")
    for tool, ver in sorted(versions.items()):
        logger.info(f"│ {tool:<{tw}} │ {ver:<{vw}} │")
    logger.info(f"└{'─' * (tw + 2)}┴{'─' * (vw + 2)}┘")
