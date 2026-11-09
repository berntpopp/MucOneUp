"""Tests for common_utils subprocess and pipeline execution."""

from __future__ import annotations

import signal
import subprocess
from unittest.mock import MagicMock, patch

import pytest

from muc_one_up.exceptions import ExternalToolError
from muc_one_up.read_simulator.utils.common_utils import run_command, run_pipeline


class TestCommonUtilsSubprocess:
    """Tests for run_command process management."""

    def test_run_command_uses_start_new_session(self):
        """run_command must start streaming subprocesses with start_new_session=True."""
        mock_proc = MagicMock()
        mock_proc.pid = 12345
        mock_proc.stdout = MagicMock()
        mock_proc.stderr = MagicMock()
        mock_proc.stdout.readline = MagicMock(return_value=b"")
        mock_proc.stderr.readline = MagicMock(return_value=b"")
        mock_proc.wait = MagicMock(return_value=0)
        mock_proc.returncode = 0

        with patch("subprocess.Popen", return_value=mock_proc) as mock_popen:
            run_command(["echo", "hi"])
            mock_popen.assert_called_once()
            assert mock_popen.call_args.kwargs.get("start_new_session") is True

    def test_run_command_escalates_to_sigkill_on_timeout(self):
        """When process does not exit after SIGTERM, run_command must escalate to SIGKILL."""
        mock_proc = MagicMock()
        mock_proc.pid = 12345
        mock_proc.stdout = MagicMock()
        mock_proc.stderr = MagicMock()
        mock_proc.stdout.readline = MagicMock(return_value=b"")
        mock_proc.stderr.readline = MagicMock(return_value=b"")
        # First wait: initial timeout; Second wait (post-SIGTERM): times out; Third wait: succeeds
        mock_proc.wait = MagicMock(
            side_effect=[
                subprocess.TimeoutExpired(cmd="test", timeout=10),
                subprocess.TimeoutExpired(cmd="test", timeout=5),
                0,
            ]
        )
        mock_proc.returncode = -9

        with (
            patch("subprocess.Popen", return_value=mock_proc),
            patch("os.killpg") as mock_killpg,
            patch("os.getpgid", return_value=12345),
            pytest.raises(ExternalToolError) as exc_info,
        ):
            run_command(["sleep", "100"], timeout=10)

        assert mock_killpg.call_count == 2
        # First call SIGTERM, second call SIGKILL
        assert mock_killpg.call_args_list[0].args == (12345, signal.SIGTERM)
        assert mock_killpg.call_args_list[1].args == (12345, signal.SIGKILL)
        assert "timed out" in str(exc_info.value).lower()


class TestRunPipeline:
    """Tests for run_pipeline execution."""

    def test_upstream_stderr_goes_to_temp_file(self):
        """Upstream stages write stderr to a temp file (not an unread PIPE); last stage is PIPE."""
        mock_proc1 = MagicMock()
        mock_proc1.stdout = MagicMock()
        mock_proc1.wait = MagicMock(return_value=0)
        mock_proc1.returncode = 0

        mock_proc2 = MagicMock()
        mock_proc2.communicate = MagicMock(return_value=(b"out", b"err"))
        mock_proc2.returncode = 0

        with patch("subprocess.Popen", side_effect=[mock_proc1, mock_proc2]) as mock_popen:
            run_pipeline([["cmd1"], ["cmd2"]])

            assert mock_popen.call_count == 2
            # First stage writes stderr to a file object, never an unread PIPE
            upstream_stderr = mock_popen.call_args_list[0].kwargs["stderr"]
            assert upstream_stderr not in (None, subprocess.PIPE)
            assert hasattr(upstream_stderr, "fileno")
            # Last stage has stderr=subprocess.PIPE
            assert mock_popen.call_args_list[1].kwargs["stderr"] == subprocess.PIPE

    def test_pipeline_timeout_escalates_to_sigkill(self):
        """When pipeline times out and processes do not terminate, escalate to SIGKILL."""
        mock_proc1 = MagicMock()
        mock_proc1.pid = 111
        mock_proc1.wait = MagicMock(
            side_effect=[subprocess.TimeoutExpired(cmd="cmd1", timeout=5), 0]
        )
        mock_proc1.returncode = -9

        mock_proc2 = MagicMock()
        mock_proc2.pid = 222
        mock_proc2.communicate = MagicMock(
            side_effect=subprocess.TimeoutExpired(cmd="cmd2", timeout=10)
        )
        mock_proc2.wait = MagicMock(
            side_effect=[subprocess.TimeoutExpired(cmd="cmd2", timeout=5), 0]
        )
        mock_proc2.returncode = -9

        with (
            patch("subprocess.Popen", side_effect=[mock_proc1, mock_proc2]),
            patch("os.killpg") as mock_killpg,
            patch("os.getpgid", side_effect=lambda pid: pid),
            pytest.raises(ExternalToolError) as exc_info,
        ):
            run_pipeline([["cmd1"], ["cmd2"]], timeout=10)

        # Check SIGKILL was sent to both processes
        calls = [(call.args[0], call.args[1]) for call in mock_killpg.call_args_list]
        assert (111, signal.SIGKILL) in calls
        assert (222, signal.SIGKILL) in calls
        assert "timed out" in str(exc_info.value).lower()

    def test_pipeline_error_preserves_final_stage_stderr(self):
        """ExternalToolError should contain the decoded stderr from final stage."""
        mock_proc1 = MagicMock()
        mock_proc1.wait = MagicMock(return_value=0)
        mock_proc1.returncode = 0

        mock_proc2 = MagicMock()
        mock_proc2.communicate = MagicMock(return_value=(b"", b"samtools: memory exhausted"))
        mock_proc2.returncode = 1

        with patch("subprocess.Popen", side_effect=[mock_proc1, mock_proc2]):
            with pytest.raises(ExternalToolError) as exc_info:
                run_pipeline([["bwa", "mem"], ["samtools", "sort"]])

            assert "samtools: memory exhausted" in exc_info.value.stderr

    def test_upstream_stage_failure_reports_its_own_stderr(self):
        """The failing upstream stage's stderr is in the error, not lost (audit #14)."""
        with pytest.raises(ExternalToolError) as exc_info:
            run_pipeline([["sh", "-c", "echo boom >&2; exit 3"], ["cat"]])

        assert "Pipeline stage 0 (sh -c echo boom >&2; exit 3) failed with exit code 3" in (
            exc_info.value.stderr
        )
        assert "boom" in exc_info.value.stderr

    def test_downstream_failure_reported_over_upstream_sigpipe(self):
        """An upstream SIGPIPE death must not mask the real downstream failure (#165)."""
        with pytest.raises(ExternalToolError) as exc_info:
            run_pipeline([["yes"], ["sh", "-c", "echo No space left on device >&2; exit 1"]])

        assert "Pipeline stage 1" in exc_info.value.stderr
        assert "No space left" in exc_info.value.stderr


@pytest.mark.parametrize("mode", ["capture", "streaming", "pipeline"])
def test_ctrl_c_kills_child_process_group(tmp_path, mode):
    """SIGINT during a tool run kills the tool instead of orphaning it (#162)."""
    import os
    import threading
    import time

    pid_file = tmp_path / "pid"
    cmd = ["sh", "-c", f"echo $$ > {pid_file}; exec sleep 30"]
    threading.Timer(1, os.kill, (os.getpid(), signal.SIGINT)).start()

    with pytest.raises(KeyboardInterrupt):
        if mode == "pipeline":
            run_pipeline([cmd, ["cat"]])
        else:
            run_command(cmd, capture=mode == "capture")

    child = int(pid_file.read_text())
    for _ in range(50):
        try:
            os.kill(child, 0)
        except ProcessLookupError:
            return
        time.sleep(0.1)
    pytest.fail(f"child {child} survived Ctrl-C")
