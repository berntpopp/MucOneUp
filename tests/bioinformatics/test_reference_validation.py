"""Tests for reference validation including assembly management (Issue #28)."""

import pytest

from muc_one_up.bioinformatics.reference_validation import (
    get_muc1_region_for_assembly,
    get_reference_path_for_assembly,
    validate_bam_file,
    validate_bed_file,
    validate_reference_for_assembly,
    validate_reference_genome,
)
from muc_one_up.exceptions import FileOperationError, ValidationError


class TestReferenceGenomeValidation:
    """Tests for basic reference genome validation."""

    def test_validate_reference_genome_exists(self, tmp_path):
        """Test validation passes when reference exists."""
        ref_file = tmp_path / "test.fa"
        ref_file.write_text(">chr1\nATCG\n")

        # Should return warnings about missing indices but not raise
        warnings = validate_reference_genome(ref_file, aligner="bwa")
        assert isinstance(warnings, list)
        assert len(warnings) > 0  # Missing indices

    def test_validate_reference_genome_missing(self, tmp_path):
        """Test validation fails when reference doesn't exist."""
        ref_file = tmp_path / "nonexistent.fa"

        with pytest.raises(FileOperationError, match="Reference genome not found"):
            validate_reference_genome(ref_file)

    def test_validate_reference_genome_unknown_aligner(self, tmp_path):
        """Test validation fails with unknown aligner."""
        ref_file = tmp_path / "test.fa"
        ref_file.write_text(">chr1\nATCG\n")

        with pytest.raises(ValidationError, match="Unknown aligner"):
            validate_reference_genome(ref_file, aligner="unknown")


class TestAssemblyManagement:
    """Test assembly-specific reference management (Issue #28)."""

    def test_get_reference_path_default_assembly(self, tmp_path):
        """Test getting reference path for default assembly."""
        ref_file = tmp_path / "hg38.fa"
        ref_file.write_text(">chr1\nATCG\n")

        config = {
            "reference_genomes": {
                "hg38": {"fasta_path": str(ref_file), "vntr_region": "chr1:155188487-155192239"}
            },
            "reference_assembly": "hg38",
        }

        path = get_reference_path_for_assembly(config)
        assert path == ref_file

    def test_get_reference_path_explicit_assembly(self, tmp_path):
        """Test getting reference path for explicitly specified assembly."""
        hg19_file = tmp_path / "hg19.fa"
        hg19_file.write_text(">chr1\nATCG\n")
        hg38_file = tmp_path / "hg38.fa"
        hg38_file.write_text(">chr1\nATCG\n")

        config = {
            "reference_genomes": {
                "hg38": {"fasta_path": str(hg38_file), "vntr_region": "chr1:155188487-155192239"},
                "hg19": {"fasta_path": str(hg19_file), "vntr_region": "chr1:155160963-155162030"},
            },
            "reference_assembly": "hg38",  # Default is hg38
        }

        # Explicitly request hg19
        path = get_reference_path_for_assembly(config, assembly="hg19")
        assert path == hg19_file

    def test_get_reference_path_missing_assembly(self):
        """Test error when assembly not configured."""
        config = {
            "reference_genomes": {
                "hg38": {"fasta_path": "/some/path.fa", "vntr_region": "chr1:1-2"}
            }
        }

        with pytest.raises(ValidationError, match="Assembly 'mm10' not found"):
            get_reference_path_for_assembly(config, assembly="mm10")

    def test_get_reference_path_no_reference_genomes(self):
        """Test error when reference_genomes section missing."""
        config = {"constants": {"hg38": {"left": "AAA", "right": "TTT"}}}

        with pytest.raises(ValidationError, match="No 'reference_genomes' section"):
            get_reference_path_for_assembly(config)

    def test_get_reference_path_file_not_found(self):
        """Test error when reference file doesn't exist."""
        config = {
            "reference_genomes": {
                "hg38": {
                    "fasta_path": "/nonexistent/file.fa",
                    "vntr_region": "chr1:1-2",
                    "source_url": "http://example.com/hg38.fa",
                }
            }
        }

        with pytest.raises(FileOperationError, match="Reference genome file not found"):
            get_reference_path_for_assembly(config, assembly="hg38")

    def test_get_reference_path_missing_fasta_path(self):
        """Test error when fasta_path not configured."""
        config = {"reference_genomes": {"hg38": {"vntr_region": "chr1:1-2"}}}

        with pytest.raises(ValidationError, match="No 'fasta_path' configured"):
            get_reference_path_for_assembly(config, assembly="hg38")

    def test_get_muc1_region_default(self):
        """Test getting MUC1 VNTR region for default assembly."""
        config = {
            "reference_genomes": {
                "hg38": {"fasta_path": "/some/path.fa", "vntr_region": "chr1:155188487-155192239"}
            },
            "reference_assembly": "hg38",
        }

        region = get_muc1_region_for_assembly(config)
        assert region == "chr1:155188487-155192239"

    def test_get_muc1_region_explicit_assembly(self):
        """Test getting MUC1 VNTR region for explicit assembly."""
        config = {
            "reference_genomes": {
                "hg38": {"fasta_path": "/some/path.fa", "vntr_region": "chr1:155188487-155192239"},
                "hg19": {"fasta_path": "/some/path.fa", "vntr_region": "chr1:155160963-155162030"},
            }
        }

        region_hg38 = get_muc1_region_for_assembly(config, "hg38")
        assert region_hg38 == "chr1:155188487-155192239"

        region_hg19 = get_muc1_region_for_assembly(config, "hg19")
        assert region_hg19 == "chr1:155160963-155162030"

    def test_get_muc1_region_missing(self):
        """Test error when vntr_region not configured."""
        config = {"reference_genomes": {"hg38": {"fasta_path": "/some/path.fa"}}}

        with pytest.raises(ValidationError, match="No 'vntr_region' configured"):
            get_muc1_region_for_assembly(config, "hg38")

    def test_validate_reference_for_assembly_calls_existing_validation(self, tmp_path):
        """Verify that validate_reference_for_assembly uses existing function."""
        ref_file = tmp_path / "hg38.fa"
        ref_file.write_text(">chr1\nATCG\n")

        config = {
            "reference_genomes": {"hg38": {"fasta_path": str(ref_file), "vntr_region": "chr1:1-2"}}
        }

        # Should call existing validation (which will warn about indices)
        warnings = validate_reference_for_assembly(config, "hg38")
        assert isinstance(warnings, list)
        # Warnings expected for missing indices
        assert len(warnings) > 0

    def test_validate_reference_for_assembly_missing_assembly(self):
        """Test validation fails for missing assembly."""
        config = {
            "reference_genomes": {
                "hg38": {"fasta_path": "/some/path.fa", "vntr_region": "chr1:1-2"}
            }
        }

        with pytest.raises(ValidationError, match="Assembly 'hg19' not found"):
            validate_reference_for_assembly(config, "hg19")

    def test_validate_reference_for_assembly_missing_file(self):
        """Test validation fails when reference file missing."""
        config = {
            "reference_genomes": {
                "hg38": {"fasta_path": "/nonexistent/file.fa", "vntr_region": "chr1:1-2"}
            }
        }

        with pytest.raises(FileOperationError):
            validate_reference_for_assembly(config, "hg38")


class TestBAMValidation:
    """Tests for BAM file validation."""

    def test_validate_bam_exists(self, tmp_path):
        """Test BAM validation when file exists."""
        bam_file = tmp_path / "test.bam"
        bam_file.write_text("dummy bam")

        # Should warn about missing index
        warnings = validate_bam_file(bam_file, require_index=False)
        assert len(warnings) == 1
        assert "BAM index missing" in warnings[0]

    def test_validate_bam_missing(self, tmp_path):
        """Test BAM validation when file doesn't exist."""
        bam_file = tmp_path / "nonexistent.bam"

        with pytest.raises(FileOperationError, match="BAM file not found"):
            validate_bam_file(bam_file)

    def test_validate_bam_require_index_missing(self, tmp_path):
        """Test BAM validation fails when index required but missing."""
        bam_file = tmp_path / "test.bam"
        bam_file.write_text("dummy bam")

        with pytest.raises(ValidationError, match="BAM index missing"):
            validate_bam_file(bam_file, require_index=True)


class TestBEDValidation:
    """Tests for BED file validation."""

    def test_validate_bed_valid(self, tmp_path):
        """Test BED validation with valid file."""
        bed_file = tmp_path / "test.bed"
        bed_file.write_text("chr1\t100\t200\n")

        # Should not raise
        validate_bed_file(bed_file)

    def test_validate_bed_missing(self, tmp_path):
        """Test BED validation when file doesn't exist."""
        bed_file = tmp_path / "nonexistent.bed"

        with pytest.raises(FileOperationError, match="BED file not found"):
            validate_bed_file(bed_file)

    def test_validate_bed_invalid_format(self, tmp_path):
        """Test BED validation with invalid format."""
        bed_file = tmp_path / "test.bed"
        bed_file.write_text("chr1\t100\n")  # Missing end position

        with pytest.raises(ValidationError, match="expected at least 3 fields"):
            validate_bed_file(bed_file)

    def test_validate_bed_invalid_coordinates(self, tmp_path):
        """Test BED validation with invalid coordinates."""
        bed_file = tmp_path / "test.bed"
        bed_file.write_text("chr1\t200\t100\n")  # start > end

        with pytest.raises(ValidationError, match="start.*>=.*end"):
            validate_bed_file(bed_file)

    def test_validate_bed_skip_comments(self, tmp_path):
        """Test BED validation skips comments and empty lines."""
        bed_file = tmp_path / "test.bed"
        bed_file.write_text("# This is a comment\n" "\n" "chr1\t100\t200\n")

        # Should not raise
        validate_bed_file(bed_file)
