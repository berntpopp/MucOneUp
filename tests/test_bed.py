"""Comprehensive tests for BED file utility module.

Tests cover:
- BED file creation
- Region definition
- VNTR BED generation
- Flanking region creation
- BED subtraction with bedtools
- Coordinate validation
- Error handling
"""

import subprocess
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from muc_one_up.read_simulator.utils.bed import (
    BedError,
    check_bedtools_available,
    create_flanking_beds,
    create_non_vntr_bed_from_capture,
    create_region_bed,
    create_vntr_bed,
)


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def vntr_region():
    """Standard VNTR region for testing."""
    return {
        "chr": "chr1",
        "start": 155188487,
        "end": 155192239,
        "name": "MUC1_VNTR",
    }


class TestCreateRegionBed:
    """Tests for basic region BED creation."""

    def test_create_simple_region(self, temp_dir):
        """Test creation of basic BED region."""
        output_path = temp_dir / "region.bed"

        result = create_region_bed(output_path, "chr1", 1000, 2000, "test_region")

        assert result == output_path
        assert output_path.exists()

        # Verify BED format: chr, start, end, name
        content = output_path.read_text()
        assert content == "chr1\t1000\t2000\ttest_region\n"

    def test_create_region_zero_start(self, temp_dir):
        """Test region starting at position 0."""
        output_path = temp_dir / "region.bed"

        result = create_region_bed(output_path, "chr2", 0, 1000, "start_zero")

        assert result == output_path
        content = output_path.read_text()
        assert content == "chr2\t0\t1000\tstart_zero\n"

    def test_invalid_negative_start(self, temp_dir):
        """Test validation of negative start position."""
        output_path = temp_dir / "region.bed"

        with pytest.raises(ValueError, match="Start position must be >= 0"):
            create_region_bed(output_path, "chr1", -1, 1000, "invalid")

    def test_invalid_end_before_start(self, temp_dir):
        """Test validation of end position."""
        output_path = temp_dir / "region.bed"

        with pytest.raises(ValueError, match="End must be > start"):
            create_region_bed(output_path, "chr1", 1000, 1000, "invalid")

        with pytest.raises(ValueError, match="End must be > start"):
            create_region_bed(output_path, "chr1", 1000, 500, "invalid")

    def test_create_region_io_error(self, temp_dir):
        """Test handling of I/O errors."""
        # Try to write to a directory path instead of file
        output_path = temp_dir

        with pytest.raises(BedError, match="Failed to create BED file"):
            create_region_bed(output_path, "chr1", 1000, 2000, "test")


class TestCreateVntrBed:
    """Tests for VNTR-specific BED creation."""

    def test_create_vntr_bed(self, temp_dir, vntr_region):
        """Test creation of VNTR BED file."""
        result = create_vntr_bed(temp_dir, vntr_region)

        assert result == temp_dir / "vntr_region.bed"
        assert result.exists()

        content = result.read_text()
        assert "chr1" in content
        assert "155188487" in content
        assert "155192239" in content
        assert "MUC1_VNTR" in content

    def test_vntr_bed_missing_keys(self, temp_dir):
        """Test validation of VNTR region dictionary."""
        incomplete_region = {"chr": "chr1", "start": 1000}

        with pytest.raises(ValueError, match="vntr_region must contain keys"):
            create_vntr_bed(temp_dir, incomplete_region)

    def test_vntr_bed_all_required_keys(self, temp_dir):
        """Test that all required keys are validated."""
        # Missing 'name'
        region = {"chr": "chr1", "start": 1000, "end": 2000}
        with pytest.raises(ValueError, match="vntr_region must contain keys"):
            create_vntr_bed(temp_dir, region)

        # Missing 'end'
        region = {"chr": "chr1", "start": 1000, "name": "test"}
        with pytest.raises(ValueError, match="vntr_region must contain keys"):
            create_vntr_bed(temp_dir, region)


class TestCreateFlankingBeds:
    """Tests for flanking region BED creation."""

    def test_create_flanking_default_size(self, temp_dir, vntr_region):
        """Test creation of flanking regions with default size."""
        left_bed, right_bed, combined_bed = create_flanking_beds(temp_dir, vntr_region)

        # Check all files created
        assert left_bed.exists()
        assert right_bed.exists()
        assert combined_bed.exists()

        # Verify left flanking region
        left_content = left_bed.read_text()
        assert "chr1" in left_content
        assert "Flanking_Left" in left_content
        # Default flanking_size is 10000
        expected_left_start = vntr_region["start"] - 10000
        assert str(expected_left_start) in left_content

        # Verify right flanking region
        right_content = right_bed.read_text()
        assert "chr1" in right_content
        assert "Flanking_Right" in right_content
        expected_right_end = vntr_region["end"] + 10000
        assert str(expected_right_end) in right_content

        # Verify combined file has both regions
        combined_content = combined_bed.read_text()
        assert combined_content.count("\n") == 2
        assert "Flanking_Left" in combined_content
        assert "Flanking_Right" in combined_content

    def test_create_flanking_custom_size(self, temp_dir, vntr_region):
        """Test creation of flanking regions with custom size."""
        left_bed, right_bed, _combined_bed = create_flanking_beds(
            temp_dir, vntr_region, flanking_size=5000
        )

        left_content = left_bed.read_text()
        expected_left_start = vntr_region["start"] - 5000
        assert str(expected_left_start) in left_content

        right_content = right_bed.read_text()
        expected_right_end = vntr_region["end"] + 5000
        assert str(expected_right_end) in right_content

    def test_flanking_left_boundary_clipping(self, temp_dir):
        """Test that left flanking region doesn't go negative."""
        # VNTR very close to chromosome start
        vntr_near_start = {
            "chr": "chr1",
            "start": 5000,
            "end": 6000,
            "name": "near_start",
        }

        left_bed, _, _ = create_flanking_beds(temp_dir, vntr_near_start, flanking_size=10000)

        left_content = left_bed.read_text()
        # Should start at 0, not -5000
        assert left_content.startswith("chr1\t0\t5000")

    def test_flanking_combined_file_format(self, temp_dir, vntr_region):
        """Test format of combined flanking BED file."""
        _, _, combined_bed = create_flanking_beds(temp_dir, vntr_region)

        lines = combined_bed.read_text().strip().split("\n")
        assert len(lines) == 2

        # First line should be left flanking
        left_fields = lines[0].split("\t")
        assert len(left_fields) == 4
        assert left_fields[3] == "Flanking_Left"

        # Second line should be right flanking
        right_fields = lines[1].split("\t")
        assert len(right_fields) == 4
        assert right_fields[3] == "Flanking_Right"


class TestCreateNonVntrBedFromCapture:
    """Tests for non-VNTR BED creation via subtraction."""

    @patch("subprocess.run")
    def test_subtraction_success(self, mock_run, temp_dir):
        """Test successful BED subtraction."""
        capture_bed = temp_dir / "capture.bed"
        vntr_bed = temp_dir / "vntr.bed"
        output_bed = temp_dir / "non_vntr.bed"

        capture_bed.write_text("chr1\t1000\t5000\tcapture\n")
        vntr_bed.write_text("chr1\t2000\t3000\tvntr\n")

        mock_run.return_value = MagicMock(returncode=0)

        result = create_non_vntr_bed_from_capture(output_bed, capture_bed, vntr_bed)

        assert result == output_bed
        mock_run.assert_called_once()
        args = mock_run.call_args[0][0]
        assert args[0] == "bedtools"
        assert args[1] == "subtract"
        assert "-a" in args
        assert "-b" in args

    @patch("subprocess.run")
    def test_subtraction_bedtools_error(self, mock_run, temp_dir):
        """Test handling of bedtools execution error."""
        capture_bed = temp_dir / "capture.bed"
        vntr_bed = temp_dir / "vntr.bed"
        output_bed = temp_dir / "non_vntr.bed"

        capture_bed.touch()
        vntr_bed.touch()

        mock_run.side_effect = subprocess.CalledProcessError(
            1, "bedtools", stderr="Subtraction failed"
        )

        with pytest.raises(BedError, match="bedtools subtract failed"):
            create_non_vntr_bed_from_capture(output_bed, capture_bed, vntr_bed)

    @patch("subprocess.run")
    def test_subtraction_bedtools_not_found(self, mock_run, temp_dir):
        """Test handling of missing bedtools."""
        capture_bed = temp_dir / "capture.bed"
        vntr_bed = temp_dir / "vntr.bed"
        output_bed = temp_dir / "non_vntr.bed"

        capture_bed.touch()
        vntr_bed.touch()

        mock_run.side_effect = FileNotFoundError()

        with pytest.raises(BedError, match="bedtools not found in PATH"):
            create_non_vntr_bed_from_capture(output_bed, capture_bed, vntr_bed)


class TestCheckBedtoolsAvailable:
    """Tests for bedtools availability checking."""

    @patch("subprocess.run")
    def test_bedtools_available(self, mock_run):
        """Test detection of available bedtools."""
        mock_run.return_value = MagicMock(returncode=0)
        assert check_bedtools_available() is True

    @patch("subprocess.run")
    def test_bedtools_not_found(self, mock_run):
        """Test handling of missing bedtools."""
        mock_run.side_effect = FileNotFoundError()
        assert check_bedtools_available() is False

    @patch("subprocess.run")
    def test_bedtools_error(self, mock_run):
        """Test handling of bedtools execution error."""
        mock_run.side_effect = subprocess.CalledProcessError(1, "bedtools")
        assert check_bedtools_available() is False

    @patch("subprocess.run")
    def test_bedtools_timeout(self, mock_run):
        """Test handling of bedtools timeout."""
        mock_run.side_effect = subprocess.TimeoutExpired("bedtools", 5)
        assert check_bedtools_available() is False


class TestBedFileEdgeCases:
    """Tests for edge cases and boundary conditions."""

    def test_large_coordinates(self, temp_dir):
        """Test handling of large genomic coordinates."""
        output_path = temp_dir / "large.bed"

        # Test with realistic large coordinates (end of chr1)
        result = create_region_bed(output_path, "chr1", 248000000, 248956422, "chr1_end")

        assert result.exists()
        content = result.read_text()
        assert "248000000" in content
        assert "248956422" in content

    def test_multiple_chromosomes(self, temp_dir):
        """Test BED file with different chromosome formats."""
        # Numeric format
        bed1 = temp_dir / "chr_numeric.bed"
        create_region_bed(bed1, "chr1", 1000, 2000, "numeric")
        assert "chr1" in bed1.read_text()

        # Non-numeric format
        bed2 = temp_dir / "chr_x.bed"
        create_region_bed(bed2, "chrX", 1000, 2000, "x_chr")
        assert "chrX" in bed2.read_text()

        # Mitochondrial
        bed3 = temp_dir / "chr_mt.bed"
        create_region_bed(bed3, "chrM", 1000, 2000, "mito")
        assert "chrM" in bed3.read_text()

    def test_special_characters_in_name(self, temp_dir):
        """Test BED region names with special characters."""
        output_path = temp_dir / "special.bed"

        # Names with underscores, hyphens, periods
        result = create_region_bed(output_path, "chr1", 1000, 2000, "Region_1.2-test")

        assert result.exists()
        content = result.read_text()
        assert "Region_1.2-test" in content
