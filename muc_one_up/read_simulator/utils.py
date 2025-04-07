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
import sys
import threading
from typing import Dict, List, Optional, Union


def run_command(
    cmd: Union[List[str], str], shell: bool = False, timeout: Optional[int] = None
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
    except Exception:
        logging.exception("Failed to start command: %s", cmd_str)
        sys.exit(1)

    def log_stream(stream, log_func):
        """Read lines from a stream and log them using log_func."""
        # Handle binary stream and decode with error handling
        for line_bytes in iter(stream.readline, b""):
            try:
                line = line_bytes.decode("utf-8", errors="replace").rstrip()
                log_func(line)
            except Exception as decode_error:
                log_func(f"[Error decoding output: {decode_error}]")
        stream.close()

    stdout_thread = threading.Thread(
        target=log_stream, args=(proc.stdout, logging.info)
    )
    stderr_thread = threading.Thread(
        target=log_stream, args=(proc.stderr, logging.error)
    )
    stdout_thread.start()
    stderr_thread.start()

    try:
        proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        logging.warning(
            "Command timed out after %s seconds. Killing process group.", timeout
        )
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
        except Exception:
            logging.exception("Error killing process group")
        proc.wait()
    stdout_thread.join()
    stderr_thread.join()

    if proc.returncode != 0:
        logging.error(
            "Command exited with non-zero exit code %d: %s", proc.returncode, cmd_str
        )
        sys.exit(proc.returncode)
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


def check_external_tools(tools: Dict[str, str]) -> None:
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
            logging.error(f"Required tool '{tool}' not found in configuration.")
            sys.exit(1)

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
        if not (shutil.which(cmd) or os.path.isfile(cmd)):
            logging.error(
                f"Required tool '{tool}' (command: '{cmd}') not found in PATH."
            )
            sys.exit(1)

    # Skip detailed validation if using conda (it will be checked when used)
    # Basic validation for specific tools
    try:
        # Only check samtools as it's usually the most reliable
        logging.info("Checking samtools availability...")
        run_command([tools["samtools"], "--version"], timeout=10)

        # Don't check BWA directly as it gives an error code when run without arguments
        # which is actually expected behavior

        logging.info("External tool validation completed successfully.")
    except SystemExit:
        logging.error("External tool validation failed.")
        logging.info(
            "If using conda/mamba environments, ensure they're properly set up"
        )
        sys.exit(1)


def cleanup_files(file_list: List[str]) -> None:
    """
    Remove files in the provided list if they exist.

    Args:
        file_list: List of filenames.
    """
    for file in file_list:
        if file and os.path.exists(file):
            try:
                os.remove(file)
                logging.info("Removed intermediate file: %s", file)
            except Exception as e:
                logging.warning("Failed to remove file %s: %s", file, str(e))
