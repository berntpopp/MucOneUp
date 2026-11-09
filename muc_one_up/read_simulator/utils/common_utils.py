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
import tempfile
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ...exceptions import ExternalToolError, ValidationError
from ..command_utils import build_tool_command


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

    def __hash__(self) -> int:
        return hash((self.returncode, self.stdout, self.stderr, self.command))


def _kill_process_group(proc: subprocess.Popen) -> None:
    """SIGTERM the process group (start_new_session=True: pgid == pid), then SIGKILL it.

    The SIGKILL is sent even if the leader exited, so grandchildren die too.
    """
    with contextlib.suppress(Exception):
        os.killpg(proc.pid, signal.SIGTERM)
    with contextlib.suppress(subprocess.TimeoutExpired):
        proc.wait(timeout=5)
    with contextlib.suppress(Exception):
        os.killpg(proc.pid, signal.SIGKILL)
    with contextlib.suppress(Exception):
        proc.wait()


def _communicate(
    cmd: list[str], cmd_str: str, timeout: float | None, cwd: Path | None, stdout: Any
) -> tuple[int, str | None, str]:
    """Run *cmd* in its own process group, collecting stderr (and stdout if PIPE).

    On timeout the whole process group is killed, so grandchildren (e.g. tools
    started via ``mamba run``) do not outlive the call.
    """
    try:
        proc = subprocess.Popen(
            cmd,
            shell=False,
            stdout=stdout,
            stderr=subprocess.PIPE,
            text=True,
            errors="replace",
            cwd=cwd,
            start_new_session=True,
        )
    except FileNotFoundError:
        raise  # Let callers handle "tool not found" directly
    except Exception as e:
        raise ExternalToolError(tool="command", exit_code=1, stderr=str(e), cmd=cmd_str) from e
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired as e:
        _kill_process_group(proc)
        with contextlib.suppress(Exception):
            proc.communicate()
        raise ExternalToolError(
            tool="command",
            exit_code=-1,
            stderr=f"Command timed out after {timeout}s",
            cmd=cmd_str,
        ) from e
    return proc.returncode, out, err or ""


def run_command(
    cmd: list[str],
    timeout: float | None = None,
    stderr_log_level: int = logging.ERROR,
    stderr_prefix: str = "",
    capture: bool = False,
    cwd: Path | None = None,
    stdout_path: Path | None = None,
    failure_log_level: int = logging.ERROR,
) -> RunResult:
    """
    Run a command in its own process group so that it can be killed on timeout.
    Capture both stdout and stderr live and log them line-by-line.

    When *capture=True*, collects stdout/stderr as strings instead of
    streaming them to the logger.

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
        failure_log_level: Logging level for the non-zero exit message (default:
            ERROR). Use DEBUG for probes where a non-zero exit is expected.

    Returns:
        A RunResult with returncode (and captured output when capture=True).

    Raises:
        ExternalToolError: If the process fails to start or returns non-zero.
    """
    if not isinstance(cmd, list):
        raise TypeError("cmd must be a list of strings, not a string (security requirement)")

    cmd_str = " ".join(cmd)
    logging.info("Running command: %s (timeout=%s)", cmd_str, timeout)

    # ---- stdout-to-file mode (no memory buffering) or capture mode ----
    if stdout_path is not None or capture:
        if stdout_path is not None:
            with Path(stdout_path).open("w") as stdout_fh:
                returncode, stdout, stderr = _communicate(cmd, cmd_str, timeout, cwd, stdout_fh)
        else:
            returncode, stdout, stderr = _communicate(cmd, cmd_str, timeout, cwd, subprocess.PIPE)
        if returncode != 0:
            logging.log(failure_log_level, "Command exited with code %d: %s", returncode, cmd_str)
            raise ExternalToolError(
                tool="command",
                exit_code=returncode,
                stderr=stderr or "Command failed",
                stdout=stdout or "",
                cmd=cmd_str,
            )
        return RunResult(returncode=returncode, stdout=stdout, stderr=stderr, command=cmd_str)

    # ---- streaming mode: log stdout/stderr line-by-line via threads ----
    try:
        popen_proc = subprocess.Popen(
            cmd,
            shell=False,  # SECURITY: Never use shell=True
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            # Use binary mode and handle encoding manually to avoid decode errors
            universal_newlines=False,
            start_new_session=True,  # Unix only
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

    timed_out = False
    try:
        popen_proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        logging.warning("Command timed out after %s seconds. Killing process group.", timeout)
        _kill_process_group(popen_proc)
    stdout_thread.join()
    stderr_thread.join()

    if timed_out or popen_proc.returncode != 0:
        if timed_out:
            logging.info("Command terminated due to timeout (%d seconds): %s", timeout, cmd_str)
            stderr_msg = f"Command timed out after {timeout}s"
        else:
            logging.log(
                failure_log_level,
                "Command exited with non-zero exit code %d: %s",
                popen_proc.returncode,
                cmd_str,
            )
            stderr_msg = f"Command failed with exit code {popen_proc.returncode}"
        raise ExternalToolError(
            tool="command",
            exit_code=popen_proc.returncode,
            stderr=stderr_msg,
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
    timeout: float | None = None,
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
    # Upstream stderr goes to temp files (a PIPE nobody reads could deadlock)
    upstream_stderr = [tempfile.TemporaryFile() for _ in cmds[:-1]]
    try:
        # Start all processes in the pipeline
        for i, cmd in enumerate(cmds):
            stdin_source = processes[-1].stdout if i > 0 else None
            stderr_dest = subprocess.PIPE if (i == len(cmds) - 1) else upstream_stderr[i]

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
            with contextlib.suppress(Exception):
                if last.stdout:
                    last.stdout.close()
            with contextlib.suppress(Exception):
                if last.stderr:
                    last.stderr.close()
            pgids = {}
            for p in processes:
                with contextlib.suppress(Exception):
                    pgids[p] = os.getpgid(p.pid)
                with contextlib.suppress(Exception):
                    os.killpg(pgids.get(p, p.pid), signal.SIGTERM)
            for p in processes:
                with contextlib.suppress(subprocess.TimeoutExpired):
                    p.wait(timeout=5)
                with contextlib.suppress(Exception):
                    target_pgid = pgids.get(p)
                    if target_pgid is None:
                        target_pgid = os.getpgid(p.pid)
                    os.killpg(target_pgid, signal.SIGKILL)
                with contextlib.suppress(Exception):
                    p.wait()
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
                try:
                    p.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    with contextlib.suppress(Exception):
                        os.killpg(os.getpgid(p.pid), signal.SIGKILL)
                    with contextlib.suppress(Exception):
                        p.wait()

        # Check for failures
        for i, p in enumerate(processes):
            if p.returncode != 0:
                if i < len(upstream_stderr):
                    upstream_stderr[i].seek(0)
                    stderr_text = upstream_stderr[i].read().decode("utf-8", errors="replace")
                else:
                    stderr_text = (stderr_bytes or b"").decode("utf-8", errors="replace")
                stage_msg = (
                    f"Pipeline stage {i} ({' '.join(cmds[i])}) failed with exit code {p.returncode}"
                )
                stderr = f"{stage_msg}: {stderr_text.strip()}" if stderr_text.strip() else stage_msg
                raise ExternalToolError(
                    tool="pipeline",
                    exit_code=p.returncode,
                    stderr=stderr,
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
    finally:
        for fh in upstream_stderr:
            fh.close()


def check_external_tools(
    tools: dict[str, str],
    required_tools: list[str] | None = None,
) -> None:
    """
    Check that all external tools required for read simulation are available.

    For each tool (e.g. reseq, faToTwoBit, samtools, pblat, bwa), this function:
      - Validates that the executable is in the system PATH or accessible via conda/mamba.
      - In some cases, makes a simple test call to ensure the tool works correctly.

    Args:
        tools: Dictionary mapping tool names to executable paths/commands.
        required_tools: Optional subset of required tool keys to validate.

    Raises:
        ValidationError: If any required tool is missing or misconfigured.
    """
    if required_tools is None:
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

        cmd = build_tool_command(command)[0]  # Get just the executable part
        if not (shutil.which(cmd) or Path(cmd).is_file()):
            raise ValidationError(f"Required tool '{tool}' (command: '{cmd}') not found in PATH")

    # Skip detailed validation if using conda (it will be checked when used)
    # Basic validation for specific tools
    try:
        # Only check samtools as it's usually the most reliable
        logging.info("Checking samtools availability...")
        run_command(build_tool_command(tools["samtools"], "--version"), timeout=10)

        # Don't check BWA directly as it gives an error code when run without arguments
        # which is actually expected behavior

        logging.info("External tool validation completed successfully.")
    except Exception as e:
        raise ValidationError(
            f"External tool validation failed. If using conda/mamba environments, ensure they're properly set up: {e}"
        ) from e
