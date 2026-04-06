"""Comprehensive tests for samtools utility module.

Tests cover:
- Tool availability checking
- BAM validation
- Read extraction by region
- Downsampling operations
- BAM merging and sorting
- Coverage calculations
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
    calculate_mean_coverage,
    check_samtools_available,
    downsample_bam,
    extract_reads_by_region,
    get_bam_read_count,
    index_bam,
    merge_bams,
    sort_bam,
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


class TestExtractReadsByRegion:
    """Tests for read extraction by region."""

    @patch(MOCK_TARGET)
    def test_extract_reads_success(self, mock_run, mock_bam_file, mock_bed_file, temp_dir):
        """Test successful read extraction."""
        output_bam = temp_dir / "output.bam"
        mock_run.return_value = _ok()

        extract_reads_by_region(mock_bam_file, mock_bed_file, output_bam)

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert "samtools" in cmd
        assert "view" in cmd
        assert "-L" in cmd
        assert str(mock_bed_file) in cmd

    @patch(MOCK_TARGET)
    def test_extract_reads_failure(self, mock_run, mock_bam_file, mock_bed_file, temp_dir):
        """Test handling of extraction failure."""
        output_bam = temp_dir / "output.bam"
        mock_run.side_effect = ExternalToolError(
            tool="command", exit_code=1, stderr="Error message", cmd="samtools view"
        )

        with pytest.raises(SamtoolsError, match="Failed to extract reads"):
            extract_reads_by_region(mock_bam_file, mock_bed_file, output_bam)

    @patch(MOCK_TARGET)
    def test_extract_reads_custom_flags(self, mock_run, mock_bam_file, mock_bed_file, temp_dir):
        """Test read extraction with custom exclude flags."""
        output_bam = temp_dir / "output.bam"
        mock_run.return_value = _ok()

        extract_reads_by_region(mock_bam_file, mock_bed_file, output_bam, exclude_flags=256)

        cmd = mock_run.call_args[0][0]
        assert "-F" in cmd
        assert "256" in cmd


class TestDownsampleBam:
    """Tests for BAM downsampling."""

    @patch(MOCK_TARGET)
    def test_downsample_success(self, mock_run, mock_bam_file, temp_dir):
        """Test successful downsampling."""
        output_bam = temp_dir / "downsampled.bam"
        mock_run.return_value = _ok()

        downsample_bam(mock_bam_file, output_bam, 0.375, seed=42)

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert "samtools" in cmd
        assert "view" in cmd
        assert "-s" in cmd
        assert "42.375" in cmd

    def test_downsample_invalid_fraction(self, mock_bam_file, temp_dir):
        """Test validation of fraction parameter."""
        output_bam = temp_dir / "downsampled.bam"

        with pytest.raises(ValueError, match="Fraction must be in"):
            downsample_bam(mock_bam_file, output_bam, 1.5)

        with pytest.raises(ValueError, match="Fraction must be in"):
            downsample_bam(mock_bam_file, output_bam, 0.0)

        with pytest.raises(ValueError, match="Fraction must be in"):
            downsample_bam(mock_bam_file, output_bam, -0.5)

    @patch(MOCK_TARGET)
    def test_downsample_fraction_formatting(self, mock_run, mock_bam_file, temp_dir):
        """Test correct formatting of downsample parameter."""
        output_bam = temp_dir / "downsampled.bam"
        mock_run.return_value = _ok()

        # Test fraction 0.5 -> seed.5
        downsample_bam(mock_bam_file, output_bam, 0.5, seed=123)
        cmd = mock_run.call_args[0][0]
        assert "123.5" in cmd

        # Test fraction 0.375 -> seed.375
        downsample_bam(mock_bam_file, output_bam, 0.375, seed=42)
        cmd = mock_run.call_args[0][0]
        assert "42.375" in cmd

    @patch(MOCK_TARGET)
    def test_downsample_failure(self, mock_run, mock_bam_file, temp_dir):
        """Test handling of downsampling failure."""
        output_bam = temp_dir / "downsampled.bam"
        mock_run.side_effect = ExternalToolError(
            tool="command", exit_code=1, stderr="Downsample error", cmd="samtools view"
        )

        with pytest.raises(SamtoolsError, match="Failed to downsample BAM"):
            downsample_bam(mock_bam_file, output_bam, 0.5)


class TestMergeBams:
    """Tests for BAM merging."""

    @patch(MOCK_TARGET)
    def test_merge_success(self, mock_run, temp_dir):
        """Test successful BAM merging."""
        bam1 = temp_dir / "input1.bam"
        bam2 = temp_dir / "input2.bam"
        output_bam = temp_dir / "merged.bam"
        bam1.touch()
        bam2.touch()
        mock_run.return_value = _ok()

        merge_bams([bam1, bam2], output_bam, threads=4)

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert "samtools" in cmd
        assert "merge" in cmd
        assert "-@" in cmd
        assert "4" in cmd

    def test_merge_insufficient_inputs(self, temp_dir):
        """Test validation of input BAM count."""
        bam1 = temp_dir / "input1.bam"
        output_bam = temp_dir / "merged.bam"

        with pytest.raises(ValueError, match="Need at least 2 BAM files"):
            merge_bams([bam1], output_bam)

    @patch(MOCK_TARGET)
    def test_merge_failure(self, mock_run, temp_dir):
        """Test handling of merge failure."""
        bam1 = temp_dir / "input1.bam"
        bam2 = temp_dir / "input2.bam"
        output_bam = temp_dir / "merged.bam"
        bam1.touch()
        bam2.touch()
        mock_run.side_effect = ExternalToolError(
            tool="command", exit_code=1, stderr="Merge error", cmd="samtools merge"
        )

        with pytest.raises(SamtoolsError, match="Failed to merge BAMs"):
            merge_bams([bam1, bam2], output_bam)


class TestSortBam:
    """Tests for BAM sorting."""

    @patch(MOCK_TARGET)
    def test_sort_success(self, mock_run, mock_bam_file, temp_dir):
        """Test successful BAM sorting."""
        output_bam = temp_dir / "sorted.bam"
        mock_run.return_value = _ok()

        sort_bam(mock_bam_file, output_bam, threads=8)

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert "samtools" in cmd
        assert "sort" in cmd
        assert "-@" in cmd
        assert "8" in cmd

    @patch(MOCK_TARGET)
    def test_sort_failure(self, mock_run, mock_bam_file, temp_dir):
        """Test handling of sort failure."""
        output_bam = temp_dir / "sorted.bam"
        mock_run.side_effect = ExternalToolError(
            tool="command", exit_code=1, stderr="Sort error", cmd="samtools sort"
        )

        with pytest.raises(SamtoolsError, match="Failed to sort BAM"):
            sort_bam(mock_bam_file, output_bam)


class TestIndexBam:
    """Tests for BAM indexing."""

    @patch(MOCK_TARGET)
    def test_index_success(self, mock_run, mock_bam_file):
        """Test successful BAM indexing."""
        mock_run.return_value = _ok()

        index_bam(mock_bam_file)

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert "samtools" in cmd
        assert "index" in cmd
        assert str(mock_bam_file) in cmd

    @patch(MOCK_TARGET)
    def test_index_failure(self, mock_run, mock_bam_file):
        """Test handling of indexing failure."""
        mock_run.side_effect = ExternalToolError(
            tool="command", exit_code=1, stderr="Index error", cmd="samtools index"
        )

        with pytest.raises(SamtoolsError, match="Failed to index BAM"):
            index_bam(mock_bam_file)


class TestCalculateMeanCoverage:
    """Tests for coverage calculation."""

    @patch(MOCK_TARGET)
    def test_coverage_calculation(self, mock_run, mock_bam_file, mock_bed_file):
        """Test successful coverage calculation."""
        # Mock samtools depth output: chr, pos, depth
        mock_run.return_value = _ok(stdout="chr1\t1000\t10\nchr1\t1001\t20\nchr1\t1002\t30\n")

        coverage = calculate_mean_coverage(mock_bam_file, mock_bed_file)

        assert coverage == 20.0  # (10 + 20 + 30) / 3
        mock_run.assert_called_once()

    @patch(MOCK_TARGET)
    def test_coverage_no_data(self, mock_run, mock_bam_file, mock_bed_file):
        """Test coverage calculation with no data."""
        mock_run.return_value = _ok(stdout="")

        coverage = calculate_mean_coverage(mock_bam_file, mock_bed_file)

        assert coverage == 0.0

    @patch(MOCK_TARGET)
    def test_includes_all_positions_flag(self, mock_run, mock_bam_file, mock_bed_file):
        """Verify samtools depth is called with -a to include zero-coverage positions."""
        mock_run.return_value = _ok(stdout="chr1\t100\t10\nchr1\t101\t0\nchr1\t102\t5\n")

        result = calculate_mean_coverage(mock_bam_file, mock_bed_file)

        # Verify -a flag is present in the command
        cmd = mock_run.call_args[0][0]
        assert "-a" in cmd, (
            "samtools depth must be called with -a to include zero-coverage positions"
        )

        # Verify mean includes the zero-coverage position: (10 + 0 + 5) / 3 = 5.0
        assert result == 5.0

    @patch(MOCK_TARGET)
    def test_coverage_failure(self, mock_run, mock_bam_file, mock_bed_file):
        """Test handling of coverage calculation failure."""
        mock_run.side_effect = ExternalToolError(
            tool="command", exit_code=1, stderr="Coverage error", cmd="samtools depth"
        )

        with pytest.raises(SamtoolsError, match="Failed to calculate coverage"):
            calculate_mean_coverage(mock_bam_file, mock_bed_file)


class TestGetBamReadCount:
    """Tests for read counting."""

    @patch(MOCK_TARGET)
    def test_read_count_success(self, mock_run, mock_bam_file):
        """Test successful read counting."""
        mock_run.return_value = _ok(stdout="12345\n")

        count = get_bam_read_count(mock_bam_file)

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
            get_bam_read_count(mock_bam_file)

    @patch(MOCK_TARGET)
    def test_read_count_invalid_output(self, mock_run, mock_bam_file):
        """Test handling of invalid read count output."""
        mock_run.return_value = _ok(stdout="invalid\n")

        with pytest.raises(SamtoolsError, match="Invalid read count output"):
            get_bam_read_count(mock_bam_file)
