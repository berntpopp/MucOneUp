"""Error-path tests for run_command and run_pipeline.

Covers non-zero exit codes, timeouts, missing binaries, and pipeline
stage failures across all three run_command modes (capture, streaming,
stdout-to-file) plus multi-stage pipelines.
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import patch

import pytest

from muc_one_up.exceptions import ExternalToolError
from muc_one_up.read_simulator.utils.common_utils import run_command, run_pipeline

# ---------------------------------------------------------------------------
# run_command — capture mode
# ---------------------------------------------------------------------------


class TestRunCommandCaptureErrors:
    """Error paths when capture=True."""

    def test_nonzero_exit_raises(self):
        with pytest.raises(ExternalToolError) as exc_info:
            run_command(
                [sys.executable, "-c", "import sys; sys.exit(42)"],
                capture=True,
            )
        assert exc_info.value.exit_code == 42

    def test_nonzero_exit_preserves_stderr(self):
        with pytest.raises(ExternalToolError) as exc_info:
            run_command(
                [sys.executable, "-c", "import sys; sys.stderr.write('boom\\n'); sys.exit(1)"],
                capture=True,
            )
        assert "boom" in exc_info.value.stderr

    def test_timeout_raises(self):
        with pytest.raises(ExternalToolError) as exc_info:
            run_command(
                [sys.executable, "-c", "import time; time.sleep(60)"],
                capture=True,
                timeout=1,
            )
        assert exc_info.value.exit_code == -1
        assert "timed out" in exc_info.value.stderr

    def test_missing_binary_raises_file_not_found(self):
        """FileNotFoundError is re-raised directly, not wrapped."""
        with pytest.raises(FileNotFoundError):
            run_command(
                ["this_binary_does_not_exist_1234567890"],
                capture=True,
            )

    def test_cmd_must_be_list(self):
        with pytest.raises(TypeError, match="must be a list"):
            run_command("echo hello", capture=True)  # type: ignore[arg-type]

    @pytest.mark.parametrize("exit_code", [1, 2, 127, 255])
    def test_various_exit_codes(self, exit_code: int):
        with pytest.raises(ExternalToolError) as exc_info:
            run_command(
                [sys.executable, "-c", f"import sys; sys.exit({exit_code})"],
                capture=True,
            )
        assert exc_info.value.exit_code == exit_code


# ---------------------------------------------------------------------------
# run_command — streaming mode (capture=False, no stdout_path)
# ---------------------------------------------------------------------------


class TestRunCommandStreamingErrors:
    """Error paths in the default streaming (Popen) mode."""

    def test_nonzero_exit_raises(self):
        with pytest.raises(ExternalToolError) as exc_info:
            run_command(
                [sys.executable, "-c", "import sys; sys.exit(7)"],
            )
        assert exc_info.value.exit_code == 7

    def test_timeout_kills_process(self):
        """Streaming timeout should kill the process and raise ExternalToolError."""
        with pytest.raises(ExternalToolError):
            run_command(
                [sys.executable, "-c", "import time; time.sleep(60)"],
                timeout=1,
            )

    def test_missing_binary_raises(self):
        """Popen wraps FileNotFoundError into ExternalToolError."""
        with pytest.raises(ExternalToolError):
            run_command(
                ["this_binary_does_not_exist_1234567890"],
            )


# ---------------------------------------------------------------------------
# run_command — stdout-to-file mode
# ---------------------------------------------------------------------------


class TestRunCommandStdoutPathErrors:
    """Error paths when stdout_path is provided."""

    def test_nonzero_exit_raises(self, tmp_path: Path):
        out = tmp_path / "out.txt"
        with pytest.raises(ExternalToolError) as exc_info:
            run_command(
                [sys.executable, "-c", "import sys; sys.exit(3)"],
                stdout_path=out,
            )
        assert exc_info.value.exit_code == 3

    def test_timeout_raises(self, tmp_path: Path):
        out = tmp_path / "out.txt"
        with pytest.raises(ExternalToolError) as exc_info:
            run_command(
                [sys.executable, "-c", "import time; time.sleep(60)"],
                stdout_path=out,
                timeout=1,
            )
        assert "timed out" in exc_info.value.stderr

    def test_missing_binary_raises_file_not_found(self, tmp_path: Path):
        out = tmp_path / "out.txt"
        with pytest.raises(FileNotFoundError):
            run_command(
                ["this_binary_does_not_exist_1234567890"],
                stdout_path=out,
            )


# ---------------------------------------------------------------------------
# run_pipeline errors
# ---------------------------------------------------------------------------


class TestRunPipelineErrors:
    """Error paths for run_pipeline."""

    def test_empty_pipeline_raises_value_error(self):
        with pytest.raises(ValueError, match="at least one command"):
            run_pipeline([], capture=True)

    def test_non_list_command_raises_type_error(self):
        with pytest.raises(TypeError, match="must be a list"):
            run_pipeline(["echo hello"], capture=True)  # type: ignore[list-item]

    def test_first_stage_failure(self):
        """First command in a two-stage pipeline fails."""
        with pytest.raises(ExternalToolError) as exc_info:
            run_pipeline(
                [
                    [sys.executable, "-c", "import sys; sys.exit(1)"],
                    [sys.executable, "-c", "import sys; sys.stdin.read()"],
                ],
                capture=True,
            )
        assert exc_info.value.exit_code == 1

    def test_second_stage_failure(self):
        """Second command fails while the first succeeds."""
        with pytest.raises(ExternalToolError) as exc_info:
            run_pipeline(
                [
                    [sys.executable, "-c", "print('data')"],
                    [sys.executable, "-c", "import sys; sys.exit(2)"],
                ],
                capture=True,
            )
        assert exc_info.value.exit_code == 2

    def test_pipeline_timeout(self):
        """Entire pipeline times out."""
        with pytest.raises(ExternalToolError) as exc_info:
            run_pipeline(
                [
                    [sys.executable, "-c", "import time; time.sleep(60)"],
                    [sys.executable, "-c", "import sys; sys.stdin.read()"],
                ],
                capture=True,
                timeout=1,
            )
        assert "timed out" in exc_info.value.stderr.lower()

    def test_single_command_pipeline_failure_delegates(self):
        """Single-command pipeline delegates to run_command, so errors propagate."""
        with pytest.raises(ExternalToolError) as exc_info:
            run_pipeline(
                [[sys.executable, "-c", "import sys; sys.exit(5)"]],
                capture=True,
            )
        assert exc_info.value.exit_code == 5

    def test_pipeline_missing_binary(self):
        """Missing binary in pipeline raises ExternalToolError."""
        with pytest.raises(ExternalToolError):
            run_pipeline(
                [
                    ["this_binary_does_not_exist_1234567890"],
                    [sys.executable, "-c", "import sys; sys.stdin.read()"],
                ],
                capture=True,
            )

    def test_pipeline_popen_generic_exception(self):
        """A generic exception during Popen startup is wrapped in ExternalToolError."""
        with patch(
            "muc_one_up.read_simulator.utils.common_utils.subprocess.Popen",
            side_effect=PermissionError("mocked permission denied"),
        ):
            with pytest.raises(ExternalToolError) as exc_info:
                run_pipeline(
                    [
                        [sys.executable, "-c", "print('a')"],
                        [sys.executable, "-c", "import sys; sys.stdin.read()"],
                    ],
                    capture=True,
                )
            assert "permission denied" in exc_info.value.stderr.lower()
