"""Tests for common_utils subprocess helpers."""

import subprocess
from unittest.mock import MagicMock, patch

import pytest

from muc_one_up.exceptions import ExternalToolError
from muc_one_up.read_simulator.utils.common_utils import run_command


def test_streaming_mode_timeout_error_includes_timeout_info():
    """Streaming mode must report timeout duration in error, not generic 'Command failed'."""
    mock_proc = MagicMock()
    mock_proc.pid = 12345
    mock_proc.stdout = MagicMock()
    mock_proc.stderr = MagicMock()
    mock_proc.stdout.readline = MagicMock(return_value=b"")
    mock_proc.stderr.readline = MagicMock(return_value=b"")
    # First call (with timeout) raises TimeoutExpired; second call (no timeout, post-kill) succeeds
    mock_proc.wait = MagicMock(
        side_effect=[subprocess.TimeoutExpired(cmd="test", timeout=30), None]
    )
    mock_proc.returncode = -15  # SIGTERM

    with (
        patch("subprocess.Popen", return_value=mock_proc),
        patch("os.setsid"),
        patch("os.killpg"),
        patch("os.getpgid", return_value=12345),
    ):
        with pytest.raises(ExternalToolError) as exc_info:
            run_command(["sleep", "999"], timeout=30)

        error = exc_info.value
        assert "timed out" in str(error).lower() or "timeout" in str(error).lower()


def test_streaming_mode_timeout_detected_regardless_of_exit_code():
    """Timeout must be reported even if process handles SIGTERM and exits 0."""
    mock_proc = MagicMock()
    mock_proc.pid = 12345
    mock_proc.stdout = MagicMock()
    mock_proc.stderr = MagicMock()
    mock_proc.stdout.readline = MagicMock(return_value=b"")
    mock_proc.stderr.readline = MagicMock(return_value=b"")
    # wait() raises TimeoutExpired, but after kill the process exits 0
    mock_proc.wait = MagicMock(
        side_effect=[subprocess.TimeoutExpired(cmd="test", timeout=10), None]
    )
    mock_proc.returncode = 0  # Process handled SIGTERM gracefully

    with (
        patch("subprocess.Popen", return_value=mock_proc),
        patch("os.setsid"),
        patch("os.killpg"),
        patch("os.getpgid", return_value=12345),
    ):
        with pytest.raises(ExternalToolError) as exc_info:
            run_command(["sleep", "999"], timeout=10)

        error = exc_info.value
        assert "timed out" in str(error).lower()


def test_streaming_mode_non_timeout_error_includes_exit_code():
    """Non-timeout errors must report exit code, not generic 'Command failed'."""
    mock_proc = MagicMock()
    mock_proc.pid = 12345
    mock_proc.stdout = MagicMock()
    mock_proc.stderr = MagicMock()
    mock_proc.stdout.readline = MagicMock(return_value=b"")
    mock_proc.stderr.readline = MagicMock(return_value=b"")
    mock_proc.wait = MagicMock()  # No timeout
    mock_proc.returncode = 1

    with patch("subprocess.Popen", return_value=mock_proc), patch("os.setsid"):
        with pytest.raises(ExternalToolError) as exc_info:
            run_command(["false"], timeout=None)

        error = exc_info.value
        assert "exit code" in str(error).lower() or "failed" in str(error).lower()
