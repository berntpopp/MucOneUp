#!/usr/bin/env python3
"""
Utility functions for the read simulation pipeline.

This module provides utility functions that are used across the read simulation
process, including process execution, file management, and tool validation.
"""

import logging
import os
import shutil
import signal
import subprocess
import threading
from pathlib import Path

from ..exceptions import ExternalToolError, ValidationError


def run_command(
    cmd: list[str] | str,
    shell: bool = False,
    timeout: int | None = None,
    stderr_log_level: int = logging.ERROR,
    stderr_prefix: str = "",
) -> int:
    """
    Run a command in its own process group so that it can be killed on timeout.
    Capture both stdout and stderr live and log them line-by-line.

    Args:
        cmd: Command (list or string) to run.
        shell: Whether to run the command in the shell.
        timeout: Timeout in seconds (None means wait indefinitely).

    Returns:
        The process return code.

    Raises:
        SystemExit: If the process fails to start or returns non-zero.
    """
    if isinstance(cmd, list) and not shell and " " in cmd[0]:
        shell = True
        cmd = " ".join(cmd)
    cmd_str = cmd if isinstance(cmd, str) else " ".join(cmd)
    logging.info("Running command: %s (timeout=%s)", cmd_str, timeout)

    try:
        proc = subprocess.Popen(
            cmd,
            shell=shell,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            # Use binary mode and handle encoding manually to avoid decode errors
            universal_newlines=False,
            preexec_fn=os.setsid,
        )
    except Exception as e:
        raise ExternalToolError(
            tool="command",
            exit_code=1,
            stderr=str(e),
            cmd=cmd_str
        ) from e

    def log_stream(stream, log_func, level=None, prefix=""):
        """Read lines from a stream and log them using log_func."""
        # Handle binary stream and decode with error handling
        for line_bytes in iter(stream.readline, b""):
            try:
                line = line_bytes.decode("utf-8", errors="replace").rstrip()
                if level is not None:
                    # Use logger.log method with specified level
                    log_func(level, f"{prefix}{line}")
                else:
                    # Use direct logging function (info/error)
                    log_func(f"{prefix}{line}")
            except Exception as decode_error:
                error_msg = f"{prefix}[Error decoding output: {decode_error}]"
                if level is not None:
                    log_func(level, error_msg)
                else:
                    log_func(error_msg)
        stream.close()

    # Log stdout always as INFO
    stdout_thread = threading.Thread(target=log_stream, args=(proc.stdout, logging.info))
    # Log stderr with configurable level
    logger = logging.getLogger()
    stderr_thread = threading.Thread(
        target=log_stream,
        args=(proc.stderr, logger.log, stderr_log_level, stderr_prefix),
    )
    stdout_thread.start()
    stderr_thread.start()

    try:
        proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        logging.warning("Command timed out after %s seconds. Killing process group.", timeout)
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
        except Exception:
            logging.exception("Error killing process group")
        proc.wait()
    stdout_thread.join()
    stderr_thread.join()

    if proc.returncode != 0:
        # Check if this is a timeout-related exit code (-15 is typical for SIGTERM)
        if timeout is not None and proc.returncode == -15:
            logging.info("Command terminated due to timeout (%d seconds): %s", timeout, cmd_str)
        else:
            logging.error(
                "Command exited with non-zero exit code %d: %s",
                proc.returncode,
                cmd_str,
            )
        raise ExternalToolError(
            tool="command",
            exit_code=proc.returncode,
            stderr="Command failed",
            cmd=cmd_str
        )
    return proc.returncode


def fix_field(field: str, desired_length: int, pad_char: str) -> str:
    """
    Ensure a field is exactly desired_length characters long.

    Args:
        field: Input string.
        desired_length: Desired length.
        pad_char: Character to pad with.

    Returns:
        String of length desired_length.
    """
    if len(field) > desired_length:
        return field[:desired_length]
    else:
        return field + pad_char * (desired_length - len(field))


def check_external_tools(tools: dict[str, str]) -> None:
    """
    Check that all external tools required for read simulation are available.

    For each tool (e.g. reseq, faToTwoBit, samtools, pblat, bwa), this function:
      - Validates that the executable is in the system PATH or accessible via conda/mamba.
      - In some cases, makes a simple test call to ensure the tool works correctly.

    Args:
        tools: Dictionary mapping tool names to executable paths/commands.

    Raises:
        SystemExit: If any required tool is missing or misconfigured.
    """
    required_tools = [
        "reseq",
        "faToTwoBit",
        "samtools",
        "pblat",
        "bwa",
    ]

    # Check if tools dictionary has all required tools
    for tool in required_tools:
        if tool not in tools:
            raise ValidationError(f"Required tool '{tool}' not found in configuration")

    # Check if the provided commands involve conda/mamba
    using_conda = any(cmd.startswith(("conda", "mamba")) for cmd in tools.values())

    if using_conda:
        logging.info("Conda/mamba environment detected for tools")
        # For conda environments, we'll do minimal checking here
        # The actual validation will happen when the tools are first used
        return

    # For non-conda tools, verify they exist in PATH
    for tool, command in tools.items():
        if tool not in required_tools:
            continue

        cmd = command.split()[0]  # Get just the executable part
        if not (shutil.which(cmd) or Path(cmd).is_file()):
            raise ValidationError(f"Required tool '{tool}' (command: '{cmd}') not found in PATH")

    # Skip detailed validation if using conda (it will be checked when used)
    # Basic validation for specific tools
    try:
        # Only check samtools as it's usually the most reliable
        logging.info("Checking samtools availability...")
        run_command([tools["samtools"], "--version"], timeout=10)

        # Don't check BWA directly as it gives an error code when run without arguments
        # which is actually expected behavior

        logging.info("External tool validation completed successfully.")
    except Exception as e:
        raise ValidationError(
            f"External tool validation failed. If using conda/mamba environments, ensure they're properly set up: {e}"
        ) from e


def cleanup_files(file_list: list[str]) -> None:
    """
    Remove files in the provided list if they exist.

    Args:
        file_list: List of filenames.
    """
    for file in file_list:
        if file:
            file_path = Path(file)
            if file_path.exists():
                try:
                    file_path.unlink()
                    logging.info("Removed intermediate file: %s", file)
                except Exception as e:
                    logging.warning("Failed to remove file %s: %s", file, str(e))
