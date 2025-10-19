"""Tests for muc_one_up.fasta_writer module.

Tests cover:
- Basic FASTA writing
- Header comment handling
- Sequence formatting (80-character lines)
- Error handling
"""

from pathlib import Path

import pytest

from muc_one_up.fasta_writer import write_fasta
from tests.utils import (
    assert_file_not_empty,
    assert_valid_fasta,
    count_fasta_records,
    get_fasta_sequences,
)


@pytest.mark.unit
class TestWriteFasta:
    """Tests for write_fasta function."""

    def test_write_fasta_single_sequence(self, tmp_path: Path):
        """Test writing single sequence to FASTA."""
        output_file = tmp_path / "output.fa"
        sequences = ["ATCGATCGATCG"]

        write_fasta(sequences, str(output_file))

        assert output_file.exists()
        assert_valid_fasta(output_file)
        assert count_fasta_records(output_file) == 1

        # Check content
        seqs = get_fasta_sequences(output_file)
        assert "haplotype_1" in seqs
        assert seqs["haplotype_1"] == "ATCGATCGATCG"

    def test_write_fasta_multiple_sequences(self, tmp_path: Path):
        """Test writing multiple sequences to FASTA."""
        output_file = tmp_path / "output.fa"
        sequences = ["ATCGATCGATCG", "GCTAGCTAGCTA", "TTTTAAAACCCC"]

        write_fasta(sequences, str(output_file))

        assert_valid_fasta(output_file)
        assert count_fasta_records(output_file) == 3

        seqs = get_fasta_sequences(output_file)
        assert len(seqs) == 3
        assert seqs["haplotype_1"] == "ATCGATCGATCG"
        assert seqs["haplotype_2"] == "GCTAGCTAGCTA"
        assert seqs["haplotype_3"] == "TTTTAAAACCCC"

    def test_write_fasta_long_sequence_wrapping(self, tmp_path: Path):
        """Test that long sequences are wrapped at 80 characters."""
        output_file = tmp_path / "output.fa"
        # Create a 200 bp sequence (should wrap to 3 lines: 80 + 80 + 40)
        long_sequence = "A" * 200

        write_fasta([long_sequence], str(output_file))

        assert_valid_fasta(output_file)

        # Read the file and check line lengths
        with output_file.open() as f:
            lines = [line.strip() for line in f if not line.startswith(">")]

        # Should have 3 sequence lines: 80 + 80 + 40
        assert len(lines) == 3
        assert len(lines[0]) == 80
        assert len(lines[1]) == 80
        assert len(lines[2]) == 40

        # Verify full sequence is preserved
        seqs = get_fasta_sequences(output_file)
        assert seqs["haplotype_1"] == long_sequence

    def test_write_fasta_with_global_comment(self, tmp_path: Path):
        """Test writing FASTA with global comment."""
        output_file = tmp_path / "output.fa"
        sequences = ["ATCGATCGATCG", "GCTAGCTAGCTA"]
        global_comment = "mutation=dupC"

        write_fasta(sequences, str(output_file), comment=global_comment)

        assert_valid_fasta(output_file)

        # Check that comment is in headers
        with output_file.open() as f:
            headers = [line for line in f if line.startswith(">")]

        assert len(headers) == 2
        assert "mutation=dupC" in headers[0]
        assert "mutation=dupC" in headers[1]

    def test_write_fasta_with_per_sequence_comments(self, tmp_path: Path):
        """Test writing FASTA with per-sequence comments."""
        output_file = tmp_path / "output.fa"
        sequences = ["ATCGATCGATCG", "GCTAGCTAGCTA"]
        comments = ["normal", "mutated at position 25"]

        write_fasta(sequences, str(output_file), comments=comments)

        assert_valid_fasta(output_file)

        # Check that specific comments are in headers
        with output_file.open() as f:
            headers = [line.strip() for line in f if line.startswith(">")]

        assert len(headers) == 2
        assert "normal" in headers[0]
        assert "mutated at position 25" in headers[1]

    def test_write_fasta_per_sequence_comments_override_global(self, tmp_path: Path):
        """Test that per-sequence comments take precedence over global comment."""
        output_file = tmp_path / "output.fa"
        sequences = ["ATCGATCGATCG", "GCTAGCTAGCTA"]
        global_comment = "global"
        per_seq_comments = ["specific1", "specific2"]

        write_fasta(sequences, str(output_file), comment=global_comment, comments=per_seq_comments)

        # Check that per-sequence comments are used, not global
        with output_file.open() as f:
            content = f.read()

        assert "specific1" in content
        assert "specific2" in content
        assert "global" not in content  # Global should not appear

    def test_write_fasta_partial_comments(self, tmp_path: Path):
        """Test writing with comments list shorter than sequences."""
        output_file = tmp_path / "output.fa"
        sequences = ["ATCGATCGATCG", "GCTAGCTAGCTA", "TTTTAAAACCCC"]
        comments = ["comment1", "comment2"]  # Missing comment for 3rd sequence

        write_fasta(sequences, str(output_file), comments=comments)

        assert_valid_fasta(output_file)

        with output_file.open() as f:
            headers = [line.strip() for line in f if line.startswith(">")]

        # First two should have comments, third should not
        assert "comment1" in headers[0]
        assert "comment2" in headers[1]
        assert len(headers[2].split()) == 1  # Only ">haplotype_3", no comment

    def test_write_fasta_none_comments_in_list(self, tmp_path: Path):
        """Test writing with None values in comments list."""
        output_file = tmp_path / "output.fa"
        sequences = ["ATCGATCGATCG", "GCTAGCTAGCTA", "TTTTAAAACCCC"]
        comments = ["comment1", None, "comment3"]  # None for 2nd sequence

        write_fasta(sequences, str(output_file), comments=comments)

        assert_valid_fasta(output_file)

        with output_file.open() as f:
            headers = [line.strip() for line in f if line.startswith(">")]

        assert "comment1" in headers[0]
        assert len(headers[1].split()) == 1  # No comment for 2nd
        assert "comment3" in headers[2]

    def test_write_fasta_custom_prefix(self, tmp_path: Path):
        """Test writing FASTA with custom header prefix."""
        output_file = tmp_path / "output.fa"
        sequences = ["ATCGATCGATCG", "GCTAGCTAGCTA"]
        custom_prefix = "chr1_allele"

        write_fasta(sequences, str(output_file), prefix=custom_prefix)

        assert_valid_fasta(output_file)

        seqs = get_fasta_sequences(output_file)
        assert "chr1_allele_1" in seqs
        assert "chr1_allele_2" in seqs

    def test_write_fasta_empty_sequence(self, tmp_path: Path):
        """Test writing empty sequence."""
        output_file = tmp_path / "output.fa"
        sequences = [""]

        write_fasta(sequences, str(output_file))

        assert output_file.exists()

        with output_file.open() as f:
            content = f.read()

        assert ">haplotype_1" in content
        # Empty sequence should result in just a header

    def test_write_fasta_sequence_exactly_80_chars(self, tmp_path: Path):
        """Test sequence that is exactly 80 characters (no wrapping needed)."""
        output_file = tmp_path / "output.fa"
        sequence = "A" * 80

        write_fasta([sequence], str(output_file))

        assert_valid_fasta(output_file)

        with output_file.open() as f:
            lines = [line.strip() for line in f if not line.startswith(">")]

        # Should be exactly one line of 80 characters
        assert len(lines) == 1
        assert len(lines[0]) == 80

    def test_write_fasta_sequence_81_chars(self, tmp_path: Path):
        """Test sequence that is 81 characters (wraps to 2 lines)."""
        output_file = tmp_path / "output.fa"
        sequence = "A" * 81

        write_fasta([sequence], str(output_file))

        with output_file.open() as f:
            lines = [line.strip() for line in f if not line.startswith(">")]

        # Should wrap to 2 lines: 80 + 1
        assert len(lines) == 2
        assert len(lines[0]) == 80
        assert len(lines[1]) == 1

    def test_write_fasta_invalid_directory(self):
        """Test writing to invalid directory raises error."""
        invalid_path = "/nonexistent/directory/output.fa"
        sequences = ["ATCGATCGATCG"]

        with pytest.raises((OSError, FileNotFoundError)):
            write_fasta(sequences, invalid_path)

    def test_write_fasta_preserves_sequence_order(self, tmp_path: Path):
        """Test that sequence order is preserved in output."""
        output_file = tmp_path / "output.fa"
        sequences = ["AAAA", "CCCC", "GGGG", "TTTT", "NNNN"]

        write_fasta(sequences, str(output_file))

        seqs = get_fasta_sequences(output_file)

        # Parse in order from file
        with output_file.open() as f:
            names_in_order = []
            for line in f:
                if line.startswith(">"):
                    names_in_order.append(line.strip()[1:].split()[0])  # Remove '>' and comments

        assert names_in_order == [
            "haplotype_1",
            "haplotype_2",
            "haplotype_3",
            "haplotype_4",
            "haplotype_5",
        ]
        assert seqs["haplotype_1"] == "AAAA"
        assert seqs["haplotype_5"] == "NNNN"

    def test_write_fasta_overwrite_existing(self, tmp_path: Path):
        """Test that writing to existing file overwrites it."""
        output_file = tmp_path / "output.fa"

        # Write first file
        write_fasta(["ATCGATCGATCG"], str(output_file))
        assert count_fasta_records(output_file) == 1

        # Overwrite with different content
        write_fasta(["GCTAGCTAGCTA", "TTTTAAAACCCC"], str(output_file))
        assert count_fasta_records(output_file) == 2

        # Verify only new content exists
        seqs = get_fasta_sequences(output_file)
        assert len(seqs) == 2
        assert "ATCGATCGATCG" not in seqs.values()
        assert "GCTAGCTAGCTA" in seqs.values()


@pytest.mark.bioinformatics
class TestFastaWriterBioinformatics:
    """Bioinformatics-specific tests for FASTA writing."""

    def test_write_fasta_realistic_sequences(
        self, tmp_path: Path, sample_haplotype_sequences: list[tuple[str, str]]
    ):
        """Test writing realistic MUC1 VNTR sequences."""
        output_file = tmp_path / "output.fa"

        # Extract just the sequences
        sequences = [seq for name, seq in sample_haplotype_sequences]

        write_fasta(sequences, str(output_file))

        assert_valid_fasta(output_file)
        assert_file_not_empty(output_file)

        seqs = get_fasta_sequences(output_file)
        assert len(seqs) == 2

        # Verify sequences are preserved exactly
        assert seqs["haplotype_1"] == sequences[0]
        assert seqs["haplotype_2"] == sequences[1]

    def test_write_fasta_with_mutation_annotations(self, tmp_path: Path):
        """Test writing FASTA with mutation annotations in comments."""
        output_file = tmp_path / "output.fa"
        sequences = ["ATCGATCGATCG", "GCTAGCTAGCTA"]
        mutation_comments = [
            "mutation=dupC position=1,25 type=insertion",
            "mutation=dupC position=2,30 type=insertion",
        ]

        write_fasta(sequences, str(output_file), comments=mutation_comments)

        assert_valid_fasta(output_file)

        with output_file.open() as f:
            headers = [line.strip() for line in f if line.startswith(">")]

        for header in headers:
            assert "mutation=" in header
            assert "position=" in header
            assert "type=" in header
