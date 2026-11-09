"""Tests for reference validation including assembly management (Issue #28)."""

import pytest

from muc_one_up.bioinformatics.reference_validation import (
    get_reference_path_for_assembly,
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

    def test_minimap2_check_matches_index_the_aligner_reuses(self, tmp_path):
        """Only {ref}.{preset}.mmi counts; a generic ref.fa.mmi is ignored by the aligner (#15)."""
        ref_file = tmp_path / "test.fa"
        ref_file.write_text(">chr1\nATCG\n")
        (tmp_path / "test.fa.fai").touch()
        (tmp_path / "test.fa.mmi").touch()

        warnings = validate_reference_genome(ref_file, aligner="minimap2", preset="map-hifi")
        assert warnings == [
            f"Minimap2 index missing: {ref_file}.map-hifi.mmi. "
            f"Run: minimap2 -x map-hifi -d {ref_file}.map-hifi.mmi {ref_file}"
        ]

        (tmp_path / "test.fa.map-hifi.mmi").touch()
        assert validate_reference_genome(ref_file, aligner="minimap2", preset="map-hifi") == []


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
