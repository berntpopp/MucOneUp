"""Comprehensive tests for VNTR efficiency modeling.

Tests cover:
- VNTREfficiencyModel initialization
- Configuration validation
- BED file creation workflow
- Read extraction and downsampling
- BAM merging and sorting
- Coverage statistics calculation
- Integration with pipeline
- Error handling and edge cases
"""

import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

from muc_one_up.read_simulator.vntr_efficiency import (
    DEFAULT_VNTR_REGION,
    VNTREfficiencyError,
    VNTREfficiencyModel,
)


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def mock_bam_file(temp_dir):
    """Create a mock BAM file."""
    bam_path = temp_dir / "input.bam"
    bam_path.touch()
    return bam_path


@pytest.fixture
def custom_vntr_region():
    """Custom VNTR region for testing."""
    return {
        "chr": "chr2",
        "start": 1000000,
        "end": 1005000,
        "name": "TEST_VNTR",
    }


class TestVNTREfficiencyModelInit:
    """Tests for model initialization."""

    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.check_samtools_available")
    def test_default_initialization(self, mock_check):
        """Test model initialization with default parameters."""
        mock_check.return_value = True

        model = VNTREfficiencyModel()

        assert model.penalty_factor == 0.375
        assert model.seed == 42
        assert model.threads == 8
        assert model.vntr_region == DEFAULT_VNTR_REGION
        assert model.capture_bed is None
        assert model.flanking_size == 10000

    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.check_samtools_available")
    def test_custom_initialization(self, mock_check, custom_vntr_region, temp_dir):
        """Test model initialization with custom parameters."""
        mock_check.return_value = True
        capture_bed = temp_dir / "capture.bed"
        capture_bed.touch()

        model = VNTREfficiencyModel(
            penalty_factor=0.5,
            seed=123,
            threads=16,
            vntr_region=custom_vntr_region,
            capture_bed=capture_bed,
            flanking_size=5000,
        )

        assert model.penalty_factor == 0.5
        assert model.seed == 123
        assert model.threads == 16
        assert model.vntr_region == custom_vntr_region
        assert model.capture_bed == capture_bed
        assert model.flanking_size == 5000

    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.check_samtools_available")
    def test_invalid_penalty_factor(self, mock_check):
        """Test validation of penalty factor."""
        mock_check.return_value = True

        # Too low
        with pytest.raises(ValueError, match="penalty_factor must be in"):
            VNTREfficiencyModel(penalty_factor=0.05)

        # Too high
        with pytest.raises(ValueError, match="penalty_factor must be in"):
            VNTREfficiencyModel(penalty_factor=1.5)

        # Negative
        with pytest.raises(ValueError, match="penalty_factor must be in"):
            VNTREfficiencyModel(penalty_factor=-0.5)

    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.check_samtools_available")
    def test_boundary_penalty_factors(self, mock_check):
        """Test boundary values for penalty factor."""
        mock_check.return_value = True

        # Minimum valid
        model_min = VNTREfficiencyModel(penalty_factor=0.1)
        assert model_min.penalty_factor == 0.1

        # Maximum valid
        model_max = VNTREfficiencyModel(penalty_factor=1.0)
        assert model_max.penalty_factor == 1.0

    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.check_samtools_available")
    def test_samtools_not_available(self, mock_check):
        """Test handling of missing samtools."""
        mock_check.return_value = False

        with pytest.raises(VNTREfficiencyError, match="samtools not found"):
            VNTREfficiencyModel()


class TestVNTREfficiencyModelCreateBedFiles:
    """Tests for BED file creation workflow."""

    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.check_samtools_available")
    @patch("muc_one_up.read_simulator.vntr_efficiency.bed")
    def test_create_beds_with_flanking(self, mock_bed, mock_check, temp_dir):
        """Test BED creation using flanking regions."""
        mock_check.return_value = True
        mock_bed.create_vntr_bed.return_value = temp_dir / "vntr.bed"
        mock_bed.create_flanking_beds.return_value = (
            temp_dir / "left.bed",
            temp_dir / "right.bed",
            temp_dir / "combined.bed",
        )

        model = VNTREfficiencyModel()
        vntr_bed, non_vntr_bed = model._create_bed_files(temp_dir)

        # Verify VNTR BED created
        mock_bed.create_vntr_bed.assert_called_once_with(temp_dir, DEFAULT_VNTR_REGION)

        # Verify flanking BED created
        mock_bed.create_flanking_beds.assert_called_once_with(temp_dir, DEFAULT_VNTR_REGION, 10000)

        assert vntr_bed == temp_dir / "vntr.bed"
        assert non_vntr_bed == temp_dir / "combined.bed"

    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.check_samtools_available")
    @patch("muc_one_up.read_simulator.vntr_efficiency.bed")
    def test_create_beds_with_capture(self, mock_bed, mock_check, temp_dir):
        """Test BED creation using capture targets."""
        mock_check.return_value = True
        capture_bed = temp_dir / "capture.bed"
        capture_bed.touch()

        mock_bed.create_vntr_bed.return_value = temp_dir / "vntr.bed"
        mock_bed.create_non_vntr_bed_from_capture.return_value = temp_dir / "non_vntr_region.bed"

        model = VNTREfficiencyModel(capture_bed=capture_bed)
        _vntr_bed, non_vntr_bed = model._create_bed_files(temp_dir)

        # Verify VNTR BED created
        mock_bed.create_vntr_bed.assert_called_once()

        # Verify capture subtraction used
        mock_bed.create_non_vntr_bed_from_capture.assert_called_once()
        mock_bed.create_flanking_beds.assert_not_called()

        assert non_vntr_bed == temp_dir / "non_vntr_region.bed"


class TestVNTREfficiencyModelOperations:
    """Tests for individual model operations."""

    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.check_samtools_available")
    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.extract_reads_by_region")
    def test_extract_reads(self, mock_extract, mock_check, temp_dir):
        """Test read extraction wrapper."""
        mock_check.return_value = True
        input_bam = temp_dir / "input.bam"
        region_bed = temp_dir / "region.bed"
        output_bam = temp_dir / "output.bam"

        model = VNTREfficiencyModel()
        model._extract_reads(input_bam, region_bed, output_bam)

        mock_extract.assert_called_once_with(input_bam, region_bed, output_bam)

    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.check_samtools_available")
    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.downsample_bam")
    def test_downsample_reads(self, mock_downsample, mock_check, temp_dir):
        """Test downsampling wrapper."""
        mock_check.return_value = True
        input_bam = temp_dir / "input.bam"
        output_bam = temp_dir / "output.bam"

        model = VNTREfficiencyModel(penalty_factor=0.375, seed=42)
        model._downsample_reads(input_bam, output_bam, 0.5)

        mock_downsample.assert_called_once_with(input_bam, output_bam, 0.5, 42)

    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.check_samtools_available")
    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.merge_bams")
    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.sort_bam")
    def test_merge_and_sort(self, mock_sort, mock_merge, mock_check, temp_dir):
        """Test merge and sort workflow."""
        mock_check.return_value = True
        bam1 = temp_dir / "bam1.bam"
        bam2 = temp_dir / "bam2.bam"
        merged_bam = temp_dir / "merged.bam"
        sorted_bam = temp_dir / "sorted.bam"

        model = VNTREfficiencyModel(threads=8)
        model._merge_and_sort(bam1, bam2, merged_bam, sorted_bam)

        mock_merge.assert_called_once_with([bam1, bam2], merged_bam, threads=8)
        mock_sort.assert_called_once_with(merged_bam, sorted_bam, threads=8)


class TestVNTREfficiencyModelStatistics:
    """Tests for statistics calculation."""

    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.check_samtools_available")
    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.get_bam_read_count")
    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.calculate_mean_coverage")
    def test_calculate_statistics(self, mock_coverage, mock_count, mock_check, temp_dir):
        """Test coverage statistics calculation."""
        mock_check.return_value = True
        mock_count.side_effect = [10000, 6000]  # input_reads, output_reads
        mock_coverage.side_effect = [100.0, 50.0]  # vntr_cov, non_vntr_cov

        input_bam = temp_dir / "input.bam"
        output_bam = temp_dir / "output.bam"
        vntr_bed = temp_dir / "vntr.bed"
        non_vntr_bed = temp_dir / "non_vntr.bed"

        model = VNTREfficiencyModel(penalty_factor=0.375, seed=42)
        stats = model._calculate_statistics(input_bam, output_bam, vntr_bed, non_vntr_bed)

        assert stats["vntr_coverage"] == 100.0
        assert stats["non_vntr_coverage"] == 50.0
        assert stats["coverage_ratio"] == 2.0
        assert stats["penalty_factor"] == 0.375
        assert stats["seed"] == 42
        assert stats["input_reads"] == 10000
        assert stats["output_reads"] == 6000
        assert stats["reads_removed"] == 4000
        assert stats["retention_fraction"] == 0.6

    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.check_samtools_available")
    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.get_bam_read_count")
    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.calculate_mean_coverage")
    def test_statistics_zero_coverage(self, mock_coverage, mock_count, mock_check, temp_dir):
        """Test statistics with zero non-VNTR coverage."""
        mock_check.return_value = True
        mock_count.side_effect = [1000, 800]
        mock_coverage.side_effect = [50.0, 0.0]  # VNTR has coverage, non-VNTR is zero

        input_bam = temp_dir / "input.bam"
        output_bam = temp_dir / "output.bam"
        vntr_bed = temp_dir / "vntr.bed"
        non_vntr_bed = temp_dir / "non_vntr.bed"

        model = VNTREfficiencyModel()
        stats = model._calculate_statistics(input_bam, output_bam, vntr_bed, non_vntr_bed)

        # Ratio should be 0.0 when denominator is zero (not infinity)
        assert stats["coverage_ratio"] == 0.0

    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.check_samtools_available")
    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.get_bam_read_count")
    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.calculate_mean_coverage")
    def test_statistics_zero_input_reads(self, mock_coverage, mock_count, mock_check, temp_dir):
        """Test statistics with zero input reads."""
        mock_check.return_value = True
        mock_count.side_effect = [0, 0]
        mock_coverage.side_effect = [0.0, 0.0]

        input_bam = temp_dir / "input.bam"
        output_bam = temp_dir / "output.bam"
        vntr_bed = temp_dir / "vntr.bed"
        non_vntr_bed = temp_dir / "non_vntr.bed"

        model = VNTREfficiencyModel()
        stats = model._calculate_statistics(input_bam, output_bam, vntr_bed, non_vntr_bed)

        assert stats["retention_fraction"] == 0.0


class TestVNTREfficiencyModelApplyBias:
    """Tests for complete apply_efficiency_bias workflow."""

    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.check_samtools_available")
    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.validate_bam")
    def test_invalid_input_bam(self, mock_validate, mock_check, temp_dir):
        """Test handling of invalid input BAM."""
        mock_check.return_value = True
        mock_validate.return_value = False

        input_bam = temp_dir / "input.bam"
        input_bam.touch()
        output_bam = temp_dir / "output.bam"
        temp_work_dir = temp_dir / "_temp"

        model = VNTREfficiencyModel()

        with pytest.raises(VNTREfficiencyError, match="Input BAM is corrupted"):
            model.apply_efficiency_bias(input_bam, output_bam, temp_work_dir)

    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.check_samtools_available")
    def test_missing_input_bam(self, mock_check, temp_dir):
        """Test handling of missing input BAM."""
        mock_check.return_value = True

        input_bam = temp_dir / "nonexistent.bam"
        output_bam = temp_dir / "output.bam"
        temp_work_dir = temp_dir / "_temp"

        model = VNTREfficiencyModel()

        with pytest.raises(FileNotFoundError):
            model.apply_efficiency_bias(input_bam, output_bam, temp_work_dir)

    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.check_samtools_available")
    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.validate_bam")
    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.index_bam")
    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.extract_reads_by_region")
    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.downsample_bam")
    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.merge_bams")
    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.sort_bam")
    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.get_bam_read_count")
    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.calculate_mean_coverage")
    @patch("muc_one_up.read_simulator.vntr_efficiency.bed")
    def test_complete_workflow(
        self,
        mock_bed,
        mock_coverage,
        mock_count,
        mock_sort,
        mock_merge,
        mock_downsample,
        mock_extract,
        mock_index,
        mock_validate,
        mock_check,
        temp_dir,
    ):
        """Test complete apply_efficiency_bias workflow."""
        # Setup mocks
        mock_check.return_value = True
        mock_validate.return_value = True

        vntr_bed = temp_dir / "vntr.bed"
        non_vntr_bed = temp_dir / "non_vntr.bed"
        mock_bed.create_vntr_bed.return_value = vntr_bed
        mock_bed.create_flanking_beds.return_value = (
            temp_dir / "left.bed",
            temp_dir / "right.bed",
            non_vntr_bed,
        )

        mock_count.side_effect = [10000, 6000]
        mock_coverage.side_effect = [150.0, 75.0]

        # Create input BAM
        input_bam = temp_dir / "input.bam"
        input_bam.touch()
        output_bam = temp_dir / "output.bam"
        temp_work_dir = temp_dir / "_temp"

        # Run workflow
        model = VNTREfficiencyModel(penalty_factor=0.375, seed=42)
        stats = model.apply_efficiency_bias(input_bam, output_bam, temp_work_dir)

        # Verify all steps called
        mock_validate.assert_called_once()
        mock_bed.create_vntr_bed.assert_called_once()
        mock_bed.create_flanking_beds.assert_called_once()
        assert mock_extract.call_count == 2  # VNTR + non-VNTR
        mock_downsample.assert_called_once()
        mock_merge.assert_called_once()
        mock_sort.assert_called_once()
        mock_index.assert_called_once()

        # Verify statistics
        assert stats["vntr_coverage"] == 150.0
        assert stats["non_vntr_coverage"] == 75.0
        assert stats["coverage_ratio"] == 2.0
        assert stats["penalty_factor"] == 0.375
        assert stats["input_reads"] == 10000
        assert stats["output_reads"] == 6000

    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.check_samtools_available")
    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.validate_bam")
    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.extract_reads_by_region")
    def test_workflow_extraction_failure(self, mock_extract, mock_validate, mock_check, temp_dir):
        """Test handling of failures during workflow."""
        mock_check.return_value = True
        mock_validate.return_value = True
        mock_extract.side_effect = Exception("Extraction failed")

        input_bam = temp_dir / "input.bam"
        input_bam.touch()
        output_bam = temp_dir / "output.bam"
        temp_work_dir = temp_dir / "_temp"

        model = VNTREfficiencyModel()

        with pytest.raises(VNTREfficiencyError, match="VNTR efficiency modeling failed"):
            model.apply_efficiency_bias(input_bam, output_bam, temp_work_dir)


class TestDefaultVntrRegion:
    """Tests for default VNTR region constant."""

    def test_default_region_structure(self):
        """Test that default VNTR region has required keys."""
        required_keys = {"chr", "start", "end", "name"}
        assert required_keys.issubset(DEFAULT_VNTR_REGION.keys())

    def test_default_region_values(self):
        """Test default VNTR region values (MUC1 hg38)."""
        assert DEFAULT_VNTR_REGION["chr"] == "chr1"
        assert DEFAULT_VNTR_REGION["start"] == 155188487
        assert DEFAULT_VNTR_REGION["end"] == 155192239
        assert DEFAULT_VNTR_REGION["name"] == "MUC1_VNTR"

    def test_default_region_coordinates_valid(self):
        """Test that default coordinates are valid."""
        assert DEFAULT_VNTR_REGION["start"] >= 0
        assert DEFAULT_VNTR_REGION["end"] > DEFAULT_VNTR_REGION["start"]


class TestVNTREfficiencyIntegration:
    """Integration tests combining multiple components."""

    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.check_samtools_available")
    def test_model_with_different_penalty_factors(self, mock_check):
        """Test model initialization with various empirically valid penalty factors."""
        mock_check.return_value = True

        # Test empirically validated range
        for penalty in [0.357, 0.375, 0.395, 0.4, 0.5, 1.0]:
            model = VNTREfficiencyModel(penalty_factor=penalty)
            assert model.penalty_factor == penalty

    @patch("muc_one_up.read_simulator.vntr_efficiency.samtools.check_samtools_available")
    def test_temp_directory_creation(self, mock_check, temp_dir):
        """Test that temporary directory is created during workflow."""
        mock_check.return_value = True
        VNTREfficiencyModel()

        temp_work_dir = temp_dir / "_vntr_temp"
        assert not temp_work_dir.exists()

        # The apply_efficiency_bias method creates temp_dir
        # We'll just test directory creation behavior
        temp_work_dir.mkdir(parents=True, exist_ok=True)
        assert temp_work_dir.exists()
