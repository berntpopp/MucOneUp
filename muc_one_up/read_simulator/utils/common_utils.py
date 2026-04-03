#!/usr/bin/env python3
"""
Utility functions for the read simulation pipeline.

This module provides utility functions that are used across the read simulation
process, including process execution, file management, and tool validation.
"""

import contextlib
import logging
import os
import shutil
import signal
import subprocess
import threading
from dataclasses import dataclass
from pathlib import Path

from ...exceptions import ExternalToolError, ValidationError


@dataclass(frozen=True, slots=True)
class RunResult:
    """Result of running an external command.

    Supports comparison with int for backward compatibility:
    callers that do ``if run_command(...) == 0:`` still work.

    When returned from capture mode, stdout and stderr are always str
    (never None). They are None only in streaming mode.
    """

    returncode: int
    stdout: str | None
    stderr: str | None
    command: str

    def __eq__(self, other: object) -> bool:
        if isinstance(other, int):
            return self.returncode == other
        if isinstance(other, RunResult):
            return (self.returncode, self.stdout, self.stderr, self.command) == (
                other.returncode,
                other.stdout,
                other.stderr,
                other.command,
            )
        return NotImplemented

    def __bool__(self) -> bool:
        return bool(self.returncode)

    def __hash__(self) -> int:
        return hash((self.returncode, self.stdout, self.stderr, self.command))


def run_command(
    cmd: list[str],
    timeout: int | None = None,
    stderr_log_level: int = logging.ERROR,
    stderr_prefix: str = "",
    capture: bool = False,
    cwd: Path | None = None,
    stdout_path: Path | None = None,
) -> RunResult:
    """
    Run a command in its own process group so that it can be killed on timeout.
    Capture both stdout and stderr live and log them line-by-line.

    When *capture=True*, uses ``subprocess.run`` to collect stdout/stderr as
    strings instead of streaming them to the logger.

    When *stdout_path* is provided, stdout is streamed directly to the file
    (no in-memory buffering) regardless of the capture flag. This avoids OOM
    for commands that produce large output (e.g., samtools fasta, minimap2).

    SECURITY: This function NEVER uses shell=True to prevent command injection.
    All commands must be provided as lists, not strings.

    Args:
        cmd: Command as a list of strings (e.g., ["bwa", "mem", "-t", "4"]).
        timeout: Timeout in seconds (None means wait indefinitely).
        stderr_log_level: Logging level for stderr output (default: ERROR).
        stderr_prefix: Prefix for stderr log messages.
        capture: If True, capture stdout/stderr instead of streaming.
        cwd: Working directory for the subprocess.
        stdout_path: If provided, stream stdout directly to this file path.

    Returns:
        A RunResult with returncode (and captured output when capture=True).

    Raises:
        ExternalToolError: If the process fails to start or returns non-zero.
    """
    if not isinstance(cmd, list):
        raise TypeError("cmd must be a list of strings, not a string (security requirement)")

    cmd_str = " ".join(cmd)
    logging.info("Running command: %s (timeout=%s)", cmd_str, timeout)

    # ---- stdout-to-file mode: stream directly to disk (no memory buffering) ----
    if stdout_path is not None:
        try:
            with Path(stdout_path).open("w") as stdout_fh:
                proc = subprocess.run(
                    cmd,
                    shell=False,
                    stdout=stdout_fh,
                    stderr=subprocess.PIPE,
                    text=True,
                    errors="replace",
                    timeout=timeout,
                    cwd=cwd,
                )
        except subprocess.TimeoutExpired as e:
            raise ExternalToolError(
                tool="command",
                exit_code=-1,
                stderr=f"Command timed out after {timeout}s",
                cmd=cmd_str,
            ) from e
        except FileNotFoundError:
            raise
        except Exception as e:
            raise ExternalToolError(tool="command", exit_code=1, stderr=str(e), cmd=cmd_str) from e

        if proc.returncode != 0:
            logging.error("Command exited with code %d: %s", proc.returncode, cmd_str)
            raise ExternalToolError(
                tool="command",
                exit_code=proc.returncode,
                stderr=proc.stderr or "Command failed",
                cmd=cmd_str,
            )
        return RunResult(
            returncode=proc.returncode,
            stdout=None,
            stderr=proc.stderr,
            command=cmd_str,
        )

    # ---- capture mode: collect output as strings ----
    if capture:
        try:
            proc = subprocess.run(
                cmd,
                shell=False,
                capture_output=True,
                text=True,
                errors="replace",
                timeout=timeout,
                cwd=cwd,
            )
        except subprocess.TimeoutExpired as e:
            raise ExternalToolError(
                tool="command",
                exit_code=-1,
                stderr=f"Command timed out after {timeout}s",
                cmd=cmd_str,
            ) from e
        except FileNotFoundError:
            raise  # Let callers handle "tool not found" directly
        except Exception as e:
            raise ExternalToolError(tool="command", exit_code=1, stderr=str(e), cmd=cmd_str) from e

        result = RunResult(
            returncode=proc.returncode,
            stdout=proc.stdout,
            stderr=proc.stderr,
            command=cmd_str,
        )
        if proc.returncode != 0:
            logging.error("Command exited with code %d: %s", proc.returncode, cmd_str)
            raise ExternalToolError(
                tool="command",
                exit_code=proc.returncode,
                stderr=proc.stderr or "Command failed",
                stdout=proc.stdout or "",
                cmd=cmd_str,
            )
        return result

    # ---- streaming mode: log stdout/stderr line-by-line via threads ----
    try:
        popen_proc = subprocess.Popen(
            cmd,
            shell=False,  # SECURITY: Never use shell=True
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            # Use binary mode and handle encoding manually to avoid decode errors
            universal_newlines=False,
            preexec_fn=os.setsid,  # Unix only
            cwd=cwd,
        )
    except Exception as e:
        raise ExternalToolError(tool="command", exit_code=1, stderr=str(e), cmd=cmd_str) from e

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
    stdout_thread = threading.Thread(target=log_stream, args=(popen_proc.stdout, logging.info))
    # Log stderr with configurable level
    logger = logging.getLogger()
    stderr_thread = threading.Thread(
        target=log_stream,
        args=(popen_proc.stderr, logger.log, stderr_log_level, stderr_prefix),
    )
    stdout_thread.start()
    stderr_thread.start()

    try:
        popen_proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        logging.warning("Command timed out after %s seconds. Killing process group.", timeout)
        try:
            os.killpg(os.getpgid(popen_proc.pid), signal.SIGTERM)  # Unix only
        except Exception:
            logging.exception("Error killing process group")
        popen_proc.wait()
    stdout_thread.join()
    stderr_thread.join()

    if popen_proc.returncode != 0:
        # Check if this is a timeout-related exit code (-15 is typical for SIGTERM)
        if timeout is not None and popen_proc.returncode == -15:
            logging.info("Command terminated due to timeout (%d seconds): %s", timeout, cmd_str)
        else:
            logging.error(
                "Command exited with non-zero exit code %d: %s",
                popen_proc.returncode,
                cmd_str,
            )
        raise ExternalToolError(
            tool="command",
            exit_code=popen_proc.returncode,
            stderr="Command failed",
            cmd=cmd_str,
        )
    return RunResult(
        returncode=popen_proc.returncode,
        stdout=None,
        stderr=None,
        command=cmd_str,
    )


def run_pipeline(
    cmds: list[list[str]],
    *,
    capture: bool = True,
    timeout: int | None = None,
    cwd: Path | None = None,
) -> RunResult:
    """Run a pipeline of commands connected by pipes.

    Connects stdout of each process to stdin of the next.
    Returns the RunResult from the final process.

    Args:
        cmds: List of commands, each as a list of strings.
        capture: If True, capture final stdout/stderr.
        timeout: Timeout in seconds for entire pipeline.
        cwd: Working directory for all processes.

    Returns:
        RunResult from the final process in the pipeline.

    Raises:
        ExternalToolError: If any process in the pipeline fails.
        ValueError: If cmds is empty.
        TypeError: If any command is not a list.
    """
    if not cmds:
        raise ValueError("Pipeline must have at least one command")

    for cmd in cmds:
        if not isinstance(cmd, list):
            raise TypeError("Each command must be a list of strings")

    pipeline_str = " | ".join(" ".join(cmd) for cmd in cmds)
    logging.info("Running pipeline: %s (timeout=%s)", pipeline_str, timeout)

    # Single command: delegate to run_command
    if len(cmds) == 1:
        return run_command(cmds[0], timeout=timeout, capture=capture, cwd=cwd)

    processes: list[subprocess.Popen] = []
    try:
        # Start all processes in the pipeline
        for i, cmd in enumerate(cmds):
            stdin_source = processes[-1].stdout if i > 0 else None
            # Only capture stderr from last process
            stderr_dest = subprocess.PIPE if (i == len(cmds) - 1) else subprocess.DEVNULL

            proc = subprocess.Popen(
                cmd,
                shell=False,
                stdin=stdin_source,
                stdout=subprocess.PIPE,
                stderr=stderr_dest,
                cwd=cwd,
                start_new_session=True,
            )
            processes.append(proc)

            # Close previous process's stdout so it gets SIGPIPE on reader close
            if i > 0 and processes[-2].stdout:
                processes[-2].stdout.close()

        # Wait for the final process
        last = processes[-1]
        try:
            stdout_bytes, stderr_bytes = last.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            for p in processes:
                with contextlib.suppress(Exception):
                    os.killpg(os.getpgid(p.pid), signal.SIGTERM)
            last.communicate()
            raise ExternalToolError(
                tool="pipeline",
                exit_code=-1,
                stderr=f"Pipeline timed out after {timeout}s",
                cmd=pipeline_str,
            ) from None

        # Wait for upstream processes to finish (with timeout to avoid hanging)
        for p in processes[:-1]:
            try:
                p.wait(timeout=10)
            except subprocess.TimeoutExpired:
                with contextlib.suppress(Exception):
                    os.killpg(os.getpgid(p.pid), signal.SIGTERM)
                p.wait(timeout=5)

        # Check for failures
        for i, p in enumerate(processes):
            if p.returncode != 0:
                raise ExternalToolError(
                    tool="pipeline",
                    exit_code=p.returncode,
                    stderr=(
                        f"Pipeline stage {i} ({' '.join(cmds[i])}) "
                        f"failed with exit code {p.returncode}"
                    ),
                    cmd=pipeline_str,
                )

        stdout_str = (
            stdout_bytes.decode("utf-8", errors="replace") if capture and stdout_bytes else None
        )
        stderr_str = (
            stderr_bytes.decode("utf-8", errors="replace") if capture and stderr_bytes else None
        )

        return RunResult(
            returncode=0,
            stdout=stdout_str,
            stderr=stderr_str,
            command=pipeline_str,
        )
    except ExternalToolError:
        raise
    except Exception as e:
        # Clean up any running processes
        for p in processes:
            with contextlib.suppress(Exception):
                os.killpg(os.getpgid(p.pid), signal.SIGTERM)
        raise ExternalToolError(
            tool="pipeline", exit_code=1, stderr=str(e), cmd=pipeline_str
        ) from e


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
