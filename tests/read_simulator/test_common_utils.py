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

    def test_upstream_stderr_is_none(self):
        """Upstream pipeline stages must inherit stderr (stderr=None), only last stage is PIPE."""
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
            # First stage has stderr=None (inherits stderr)
            assert mock_popen.call_args_list[0].kwargs["stderr"] is None
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

    def test_upstream_pipeline_stage_failure_preserves_stage_message(self):
        """ExternalToolError should identify upstream failed stage and preserve stderr."""
        mock_proc1 = MagicMock()
        mock_proc1.wait = MagicMock(return_value=1)
        mock_proc1.returncode = 1

        mock_proc2 = MagicMock()
        mock_proc2.communicate = MagicMock(return_value=(b"", b"samtools: view error"))
        mock_proc2.returncode = 0

        with patch("subprocess.Popen", side_effect=[mock_proc1, mock_proc2]):
            with pytest.raises(ExternalToolError) as exc_info:
                run_pipeline([["bwa", "mem"], ["samtools", "view"]])

            assert "Pipeline stage 0 (bwa mem) failed with exit code 1" in exc_info.value.stderr
            assert "samtools: view error" in exc_info.value.stderr
