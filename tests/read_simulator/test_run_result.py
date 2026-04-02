"""Tests for RunResult, run_command, and run_pipeline."""

from __future__ import annotations

import sys

import pytest

from muc_one_up.read_simulator.utils.common_utils import RunResult, run_command, run_pipeline


class TestRunResult:
    def test_create(self):
        r = RunResult(returncode=0, stdout="hello", stderr="", command="echo hello")
        assert r.returncode == 0
        assert r.stdout == "hello"

    def test_none_stdout_in_streaming_mode(self):
        r = RunResult(returncode=0, stdout=None, stderr=None, command="cmd")
        assert r.stdout is None

    def test_eq_int_for_backward_compat(self):
        """RunResult should compare equal to its returncode for backward compat."""
        r = RunResult(returncode=0, stdout=None, stderr=None, command="cmd")
        assert r == 0
        assert r != 1

    def test_bool_truthy_on_zero(self):
        """RunResult with returncode 0 should be falsy like int 0 for `if result:` patterns."""
        r = RunResult(returncode=0, stdout=None, stderr=None, command="cmd")
        # Note: this follows int semantics where 0 is falsy
        assert not r  # 0 is falsy


class TestRunCommandCapture:
    def test_capture_stdout(self):
        result = run_command(
            [sys.executable, "-c", "print('hello world')"],
            capture=True,
        )
        assert isinstance(result, RunResult)
        assert result.returncode == 0
        assert "hello world" in result.stdout

    def test_capture_stderr(self):
        result = run_command(
            [sys.executable, "-c", "import sys; sys.stderr.write('err msg\\n')"],
            capture=True,
        )
        assert result.returncode == 0
        assert "err msg" in result.stderr

    def test_capture_nonzero_raises(self):
        from muc_one_up.exceptions import ExternalToolError

        with pytest.raises(ExternalToolError):
            run_command(
                [sys.executable, "-c", "import sys; sys.exit(1)"],
                capture=True,
            )

    def test_streaming_mode_returns_run_result(self):
        result = run_command(
            [sys.executable, "-c", "print('streamed')"],
        )
        assert isinstance(result, RunResult)
        assert result.returncode == 0
        assert result.stdout is None

    def test_capture_timeout(self):
        from muc_one_up.exceptions import ExternalToolError

        with pytest.raises(ExternalToolError):
            run_command(
                [sys.executable, "-c", "import time; time.sleep(30)"],
                capture=True,
                timeout=1,
            )


class TestRunPipeline:
    def test_two_command_pipe(self):
        result = run_pipeline(
            [
                [sys.executable, "-c", "print('hello\\nworld\\nfoo')"],
                [
                    sys.executable,
                    "-c",
                    "import sys; [print(l, end='') for l in sys.stdin if 'world' in l]",
                ],
            ],
            capture=True,
        )
        assert result.returncode == 0
        assert "world" in result.stdout
        assert "hello" not in result.stdout

    def test_single_command_pipeline(self):
        result = run_pipeline(
            [[sys.executable, "-c", "print('solo')"]],
            capture=True,
        )
        assert result.returncode == 0
        assert "solo" in result.stdout

    def test_pipeline_failure_in_last_stage_raises(self):
        from muc_one_up.exceptions import ExternalToolError

        with pytest.raises(ExternalToolError):
            run_pipeline(
                [
                    [sys.executable, "-c", "print('data')"],
                    [sys.executable, "-c", "import sys; sys.exit(1)"],
                ],
                capture=True,
            )

    def test_empty_pipeline_raises(self):
        with pytest.raises(ValueError):
            run_pipeline([], capture=True)

    def test_pipeline_returns_run_result(self):
        result = run_pipeline(
            [[sys.executable, "-c", "print('ok')"]],
            capture=True,
        )
        assert isinstance(result, RunResult)
