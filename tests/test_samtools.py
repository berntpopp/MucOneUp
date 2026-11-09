"""Comprehensive tests for samtools utility module.

Tests cover:
- Tool availability checking
- BAM validation
- Read counting
- Error handling
"""

import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

from muc_one_up.exceptions import ExternalToolError
from muc_one_up.read_simulator.utils.common_utils import RunResult
from muc_one_up.read_simulator.utils.samtools import (
    SamtoolsError,
    check_samtools_available,
    get_read_count,
    validate_bam,
)

# All tests mock run_command at the samtools module level since
# samtools.py now delegates to common_utils.run_command instead
# of calling subprocess.run directly.
MOCK_TARGET = "muc_one_up.read_simulator.utils.samtools.run_command"


def _ok(stdout: str = "", stderr: str = "") -> RunResult:
    """Helper to build a successful RunResult."""
    return RunResult(returncode=0, stdout=stdout, stderr=stderr, command="samtools ...")


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def mock_bam_file(temp_dir):
    """Create a mock BAM file path."""
    bam_path = temp_dir / "test.bam"
    bam_path.touch()
    return bam_path


@pytest.fixture
def mock_bed_file(temp_dir):
    """Create a mock BED file."""
    bed_path = temp_dir / "test.bed"
    bed_path.write_text("chr1\t1000\t2000\tregion1\n")
    return bed_path


class TestCheckSamtoolsAvailable:
    """Tests for samtools availability checking."""

    @patch(MOCK_TARGET)
    def test_samtools_available(self, mock_run):
        """Test detection of available samtools."""
        mock_run.return_value = _ok()
        assert check_samtools_available() is True
        mock_run.assert_called_once()

    @patch(MOCK_TARGET)
    def test_samtools_not_found(self, mock_run):
        """Test handling of missing samtools."""
        mock_run.side_effect = ExternalToolError(
            tool="command", exit_code=1, stderr="not found", cmd="samtools --version"
        )
        assert check_samtools_available() is False

    @patch(MOCK_TARGET)
    def test_samtools_error(self, mock_run):
        """Test handling of samtools execution error."""
        mock_run.side_effect = ExternalToolError(
            tool="command", exit_code=1, stderr="error", cmd="samtools --version"
        )
        assert check_samtools_available() is False

    @patch(MOCK_TARGET)
    def test_samtools_timeout(self, mock_run):
        """Test handling of samtools timeout."""
        mock_run.side_effect = ExternalToolError(
            tool="command", exit_code=-1, stderr="timed out", cmd="samtools --version"
        )
        assert check_samtools_available() is False


class TestValidateBam:
    """Tests for BAM validation."""

    @patch(MOCK_TARGET)
    def test_valid_bam(self, mock_run, mock_bam_file):
        """Test validation of valid BAM file."""
        mock_run.return_value = _ok()
        assert validate_bam(mock_bam_file) is True

    @patch(MOCK_TARGET)
    def test_corrupted_bam(self, mock_run, mock_bam_file):
        """Test detection of corrupted BAM file."""
        mock_run.side_effect = ExternalToolError(
            tool="command", exit_code=1, stderr="corrupted", cmd="samtools quickcheck"
        )
        assert validate_bam(mock_bam_file) is False


class TestGetReadCount:
    """Tests for read counting."""

    @patch(MOCK_TARGET)
    def test_read_count_success(self, mock_run, mock_bam_file):
        """Test successful read counting."""
        mock_run.return_value = _ok(stdout="12345\n")

        count = get_read_count(mock_bam_file)

        assert count == 12345
        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert "samtools" in cmd
        assert "view" in cmd
        assert "-c" in cmd

    @patch(MOCK_TARGET)
    def test_read_count_failure(self, mock_run, mock_bam_file):
        """Test handling of read count failure."""
        mock_run.side_effect = ExternalToolError(
            tool="command", exit_code=1, stderr="Count error", cmd="samtools view -c"
        )

        with pytest.raises(SamtoolsError, match="Failed to count reads"):
            get_read_count(mock_bam_file)

    @patch(MOCK_TARGET)
    def test_read_count_invalid_output(self, mock_run, mock_bam_file):
        """Test handling of invalid read count output."""
        mock_run.return_value = _ok(stdout="invalid\n")

        with pytest.raises(SamtoolsError, match="Invalid read count output"):
            get_read_count(mock_bam_file)


class TestSamtoolsErrorHierarchy:
    """Tests for SamtoolsError exception hierarchy and consolidation (#148)."""

    def test_samtools_error_inherits_from_external_tool_error(self):
        """SamtoolsError must inherit from ExternalToolError and MucOneUpError."""
        from muc_one_up.exceptions import MucOneUpError

        assert issubclass(SamtoolsError, ExternalToolError)
        assert issubclass(SamtoolsError, MucOneUpError)
        err = SamtoolsError("Extraction failed: error message")
        assert isinstance(err, ExternalToolError)
        assert isinstance(err, MucOneUpError)
        assert str(err) == "Extraction failed: error message"

    def test_samtools_wrapper_reexports_consolidated_functions(self):
        """samtools_wrapper re-exports consolidated utilities."""
        import muc_one_up.read_simulator.wrappers.samtools_wrapper as st_wrap

        assert hasattr(st_wrap, "SamtoolsError")
        assert hasattr(st_wrap, "merge_bam_files")
        assert hasattr(st_wrap, "downsample_bam")
        assert hasattr(st_wrap, "downsample_entire_bam")
