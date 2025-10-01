"""Tests for muc_one_up/bioinformatics/validation.py - sequence validation."""

import pytest

from muc_one_up.bioinformatics.validation import (
    DNA_BASES,
    SNP_BASES,
    validate_dna_sequence,
    validate_fasta_format,
    validate_gc_content_range,
    validate_repeat_structure,
    validate_sequence_length,
    validate_snp_base,
    validate_snp_record,
)
from muc_one_up.exceptions import ValidationError


class TestValidateDNASequence:
    """Test DNA sequence validation."""

    def test_valid_sequence(self):
        """Valid DNA sequences should not raise errors."""
        validate_dna_sequence("ACGT")
        validate_dna_sequence("ACGTACGT")
        validate_dna_sequence("acgt")  # lowercase should work

    def test_sequence_with_n(self):
        """N for ambiguous bases should be valid when allowed."""
        validate_dna_sequence("ACGTN", allow_ambiguous=True)

    def test_sequence_with_n_not_allowed(self):
        """N should raise error when ambiguous=False."""
        with pytest.raises(ValidationError, match="Invalid DNA bases"):
            validate_dna_sequence("ACGTN", allow_ambiguous=False)

    def test_invalid_bases(self):
        """Sequences with invalid bases should raise ValidationError."""
        with pytest.raises(ValidationError, match="Invalid DNA bases"):
            validate_dna_sequence("ACGTX")

    def test_empty_sequence(self):
        """Empty sequences should raise ValidationError."""
        with pytest.raises(ValidationError, match="cannot be empty"):
            validate_dna_sequence("")

    def test_mixed_case(self):
        """Mixed case should be valid."""
        validate_dna_sequence("AcGt")


class TestValidateFastaFormat:
    """Test FASTA format validation."""

    def test_valid_fasta(self, tmp_path):
        """Valid FASTA file should not raise errors."""
        fasta = tmp_path / "test.fa"
        fasta.write_text(">sequence1\nACGT\nTGCA\n>sequence2\nGGCC\n")
        validate_fasta_format(fasta)

    def test_missing_file(self):
        """Missing file should raise FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            validate_fasta_format("/nonexistent/file.fa")

    def test_empty_file(self, tmp_path):
        """Empty FASTA file should raise ValidationError."""
        fasta = tmp_path / "empty.fa"
        fasta.write_text("")
        with pytest.raises(ValidationError, match="empty"):
            validate_fasta_format(fasta)

    def test_no_header(self, tmp_path):
        """FASTA without header should raise ValidationError."""
        fasta = tmp_path / "no_header.fa"
        fasta.write_text("ACGT\n")
        with pytest.raises(ValidationError, match="must start with header"):
            validate_fasta_format(fasta)

    def test_invalid_sequence(self, tmp_path):
        """FASTA with invalid sequence should raise ValidationError."""
        fasta = tmp_path / "invalid.fa"
        fasta.write_text(">sequence1\nACGTX\n")
        with pytest.raises(ValidationError, match="Invalid DNA sequence"):
            validate_fasta_format(fasta)


class TestValidateRepeatStructure:
    """Test repeat structure validation."""

    def test_valid_structure(self):
        """Valid repeat structures should not raise errors."""
        valid = {"1", "2", "7", "8", "9"}
        validate_repeat_structure("1-2-7-8-9", valid)

    def test_structure_with_mutation_markers(self):
        """Structures with 'm' markers should be valid."""
        valid = {"1", "2", "7"}
        validate_repeat_structure("1-2m-7", valid)

    def test_single_repeat(self):
        """Single repeat (no hyphens) should be valid."""
        valid = {"1", "2", "7"}
        validate_repeat_structure("1", valid)

    def test_invalid_symbol(self):
        """Structure with invalid symbol should raise ValidationError."""
        valid = {"1", "2", "7"}
        with pytest.raises(ValidationError, match="Invalid repeat symbol"):
            validate_repeat_structure("1-2-X", valid)

    def test_empty_structure(self):
        """Empty structure should raise ValidationError."""
        valid = {"1", "2"}
        with pytest.raises(ValidationError, match="cannot be empty"):
            validate_repeat_structure("", valid)


class TestValidateSNPBase:
    """Test SNP base validation."""

    def test_valid_bases(self):
        """A, C, G, T should be valid."""
        for base in ["A", "C", "G", "T"]:
            validate_snp_base(base)

    def test_lowercase(self):
        """Lowercase bases should be valid."""
        validate_snp_base("a")
        validate_snp_base("t")

    def test_ambiguous_base(self):
        """N should raise ValidationError."""
        with pytest.raises(ValidationError, match="Invalid SNP base"):
            validate_snp_base("N")

    def test_invalid_base(self):
        """Invalid bases should raise ValidationError."""
        with pytest.raises(ValidationError, match="Invalid SNP base"):
            validate_snp_base("X")

    def test_empty_base(self):
        """Empty base should raise ValidationError."""
        with pytest.raises(ValidationError, match="cannot be empty"):
            validate_snp_base("")


class TestValidateSNPRecord:
    """Test complete SNP record validation."""

    def test_valid_record(self):
        """Valid SNP record should not raise errors."""
        validate_snp_record(
            haplotype=1,
            position=100,
            ref="A",
            alt="G",
            num_haplotypes=2,
            sequence_length=1000,
        )

    def test_invalid_haplotype_index(self):
        """Invalid haplotype should raise ValidationError."""
        with pytest.raises(ValidationError, match="haplotype index.*out of range"):
            validate_snp_record(
                haplotype=3,
                position=100,
                ref="A",
                alt="G",
                num_haplotypes=2,
                sequence_length=1000,
            )

    def test_invalid_position(self):
        """Position out of range should raise ValidationError."""
        with pytest.raises(ValidationError, match="position.*out of range"):
            validate_snp_record(
                haplotype=1,
                position=1000,
                ref="A",
                alt="G",
                num_haplotypes=2,
                sequence_length=1000,
            )

    def test_same_ref_and_alt(self):
        """Identical ref and alt should raise ValidationError."""
        with pytest.raises(ValidationError, match="must differ"):
            validate_snp_record(
                haplotype=1,
                position=100,
                ref="A",
                alt="A",
                num_haplotypes=2,
                sequence_length=1000,
            )


class TestValidateGCContentRange:
    """Test GC content validation."""

    def test_valid_range(self):
        """GC content in [0, 100] should be valid."""
        validate_gc_content_range(0.0)
        validate_gc_content_range(50.0)
        validate_gc_content_range(100.0)

    def test_below_range(self):
        """GC content < 0 should raise ValidationError."""
        with pytest.raises(ValidationError, match="between 0.0 and 100.0"):
            validate_gc_content_range(-1.0)

    def test_above_range(self):
        """GC content > 100 should raise ValidationError."""
        with pytest.raises(ValidationError, match="between 0.0 and 100.0"):
            validate_gc_content_range(101.0)


class TestValidateSequenceLength:
    """Test sequence length validation."""

    def test_valid_length(self):
        """Length within bounds should be valid."""
        validate_sequence_length(100, min_length=10, max_length=200)
        validate_sequence_length(10, min_length=10, max_length=200)
        validate_sequence_length(200, min_length=10, max_length=200)

    def test_below_minimum(self):
        """Length below minimum should raise ValidationError."""
        with pytest.raises(ValidationError, match="below minimum"):
            validate_sequence_length(5, min_length=10)

    def test_above_maximum(self):
        """Length above maximum should raise ValidationError."""
        with pytest.raises(ValidationError, match="exceeds maximum"):
            validate_sequence_length(300, min_length=10, max_length=200)

    def test_no_maximum(self):
        """Should work without maximum constraint."""
        validate_sequence_length(1000000, min_length=10, max_length=None)


class TestDNABases:
    """Test DNA_BASES and SNP_BASES constants."""

    def test_dna_bases_content(self):
        """DNA_BASES should contain A, C, G, T, N."""
        assert {"A", "C", "G", "T", "N"} == DNA_BASES

    def test_snp_bases_content(self):
        """SNP_BASES should contain A, C, G, T (no N)."""
        assert {"A", "C", "G", "T"} == SNP_BASES
