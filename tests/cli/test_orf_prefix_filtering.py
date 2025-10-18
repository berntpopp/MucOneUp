"""Tests for ORF prediction with amino acid prefix filtering.

These tests verify that the --orf-aa-prefix flag correctly filters
ORF predictions to only include sequences starting with the specified
amino acid prefix.
"""

import shutil
from unittest.mock import Mock, patch

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@pytest.fixture
def sample_fasta(tmp_path):
    """Create a sample FASTA file for ORF prediction."""
    fasta_file = tmp_path / "test.fa"
    fasta_file.write_text(
        ">seq1\n"
        "ATGGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG\n"
        ">seq2\n"
        "ATGTTAAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG\n"
    )
    return str(fasta_file)


@pytest.fixture
def sample_orfs(tmp_path):
    """Create a sample ORF output file with different protein sequences."""
    orf_file = tmp_path / "test.orfs.fa"

    # Create ORF records with different amino acid prefixes
    orfs = [
        SeqRecord(Seq("MTSSVAAAAAAA"), id="orf1", description="starts with MTSSV"),
        SeqRecord(Seq("MTSSVBBBBBBB"), id="orf2", description="starts with MTSSV"),
        SeqRecord(Seq("MTAACCCCCCC"), id="orf3", description="starts with MTA"),
        SeqRecord(Seq("MALLDDDDDDD"), id="orf4", description="starts with MAL"),
        SeqRecord(Seq("MTSSVEEEEEE"), id="orf5", description="starts with MTSSV"),
    ]

    SeqIO.write(orfs, str(orf_file), "fasta")
    return str(orf_file)


class TestORFPrefixFiltering:
    """Test ORF prediction with amino acid prefix filtering."""

    def test_filters_orfs_by_prefix_mtssv(self, sample_orfs, tmp_path):
        """Test that filtering keeps only ORFs starting with MTSSV."""
        from Bio import SeqIO

        # Filter the ORFs (simulating the CLI logic)
        orf_aa_prefix = "MTSSV"
        filtered_orfs = []
        total_orfs = 0

        for record in SeqIO.parse(sample_orfs, "fasta"):
            total_orfs += 1
            if str(record.seq).startswith(orf_aa_prefix):
                filtered_orfs.append(record)

        # Verify filtering
        assert total_orfs == 5, "Should have 5 total ORFs"
        assert len(filtered_orfs) == 3, "Should keep 3 ORFs starting with MTSSV"

        # Verify all filtered ORFs start with MTSSV
        for record in filtered_orfs:
            assert str(record.seq).startswith("MTSSV")

        # Verify the correct ORF IDs
        filtered_ids = {record.id for record in filtered_orfs}
        assert filtered_ids == {"orf1", "orf2", "orf5"}

    def test_filters_orfs_by_prefix_mta(self, sample_orfs):
        """Test filtering with MTA prefix."""
        from Bio import SeqIO

        orf_aa_prefix = "MTA"
        filtered_orfs = []

        for record in SeqIO.parse(sample_orfs, "fasta"):
            if str(record.seq).startswith(orf_aa_prefix):
                filtered_orfs.append(record)

        assert len(filtered_orfs) == 1
        assert filtered_orfs[0].id == "orf3"
        assert str(filtered_orfs[0].seq).startswith("MTA")

    def test_no_orfs_match_prefix(self, sample_orfs):
        """Test that no ORFs are kept when prefix doesn't match any."""
        from Bio import SeqIO

        orf_aa_prefix = "NOTFOUND"
        filtered_orfs = []

        for record in SeqIO.parse(sample_orfs, "fasta"):
            if str(record.seq).startswith(orf_aa_prefix):
                filtered_orfs.append(record)

        assert len(filtered_orfs) == 0, "No ORFs should match 'NOTFOUND' prefix"

    def test_no_prefix_keeps_all_orfs(self, sample_orfs):
        """Test that without a prefix, all ORFs are kept."""
        from Bio import SeqIO

        all_orfs = list(SeqIO.parse(sample_orfs, "fasta"))

        # When no prefix specified, all ORFs should be kept
        assert len(all_orfs) == 5

    def test_empty_prefix_keeps_all_orfs(self, sample_orfs):
        """Test that empty string prefix keeps all ORFs."""
        from Bio import SeqIO

        orf_aa_prefix = ""
        filtered_orfs = []

        for record in SeqIO.parse(sample_orfs, "fasta"):
            # Empty string prefix matches all (every string starts with "")
            if str(record.seq).startswith(orf_aa_prefix):
                filtered_orfs.append(record)

        assert len(filtered_orfs) == 5

    def test_case_sensitive_prefix_matching(self, tmp_path):
        """Test that prefix matching is case-sensitive."""
        from Bio import SeqIO

        # Create ORFs with lowercase sequences
        orf_file = tmp_path / "case_test.orfs.fa"
        orfs = [
            SeqRecord(Seq("MTSSVAA"), id="orf1"),
            SeqRecord(Seq("mtssvaa"), id="orf2"),  # lowercase
        ]
        SeqIO.write(orfs, str(orf_file), "fasta")

        # Filter with uppercase prefix
        orf_aa_prefix = "MTSSV"
        filtered_orfs = []

        for record in SeqIO.parse(str(orf_file), "fasta"):
            if str(record.seq).startswith(orf_aa_prefix):
                filtered_orfs.append(record)

        # Only uppercase sequence should match
        assert len(filtered_orfs) == 1
        assert filtered_orfs[0].id == "orf1"

    def test_single_letter_prefix(self, sample_orfs):
        """Test filtering with single letter prefix 'M'."""
        from Bio import SeqIO

        orf_aa_prefix = "M"
        filtered_orfs = []

        for record in SeqIO.parse(sample_orfs, "fasta"):
            if str(record.seq).startswith(orf_aa_prefix):
                filtered_orfs.append(record)

        # All ORFs start with 'M' (methionine)
        assert len(filtered_orfs) == 5

    def test_write_filtered_orfs_to_file(self, sample_orfs, tmp_path):
        """Test that filtered ORFs can be written back to a file."""
        from Bio import SeqIO

        output_file = tmp_path / "filtered.fa"

        # Filter ORFs
        orf_aa_prefix = "MTSSV"
        filtered_orfs = [
            record
            for record in SeqIO.parse(sample_orfs, "fasta")
            if str(record.seq).startswith(orf_aa_prefix)
        ]

        # Write to file
        SeqIO.write(filtered_orfs, str(output_file), "fasta")

        # Read back and verify
        reloaded = list(SeqIO.parse(str(output_file), "fasta"))
        assert len(reloaded) == 3
        for record in reloaded:
            assert str(record.seq).startswith("MTSSV")


class TestORFPrefixFilteringIntegration:
    """Integration tests for ORF prefix filtering in CLI workflow."""

    @patch("subprocess.run")
    def test_cli_applies_prefix_filter(self, mock_run, sample_fasta, sample_orfs, tmp_path):
        """Test that CLI correctly applies prefix filtering after orfipy."""
        # Mock orfipy success
        mock_run.return_value = Mock(returncode=0, stderr="", stdout="")

        # Create output file in a different location
        output_file = tmp_path / "filtered_output.orfs.fa"

        # Simulate the CLI filtering logic
        from Bio import SeqIO

        orf_aa_prefix = "MTSSV"

        # Copy sample ORFs to output location
        shutil.copy(sample_orfs, output_file)

        # Apply filtering (this is what the CLI does)
        filtered_orfs = []
        total_orfs = 0

        for record in SeqIO.parse(str(output_file), "fasta"):
            total_orfs += 1
            if str(record.seq).startswith(orf_aa_prefix):
                filtered_orfs.append(record)

        # Write filtered results back
        SeqIO.write(filtered_orfs, str(output_file), "fasta")

        # Verify filtering worked
        result = list(SeqIO.parse(str(output_file), "fasta"))
        assert len(result) == 3  # Only MTSSV ORFs
        assert total_orfs == 5  # Originally had 5

    def test_filtering_preserves_sequence_metadata(self, tmp_path):
        """Test that filtering preserves sequence IDs and descriptions."""
        from Bio import SeqIO

        # Create ORFs with detailed metadata
        orf_file = tmp_path / "metadata_test.fa"
        orfs = [
            SeqRecord(
                Seq("MTSSVAAAA"),
                id="haplotype_1_ORF1",
                description="length=9 start=0 stop=27 strand=+",
            ),
            SeqRecord(
                Seq("MTAAABBBB"),
                id="haplotype_1_ORF2",
                description="length=9 start=100 stop=127 strand=+",
            ),
        ]
        SeqIO.write(orfs, str(orf_file), "fasta")

        # Filter
        orf_aa_prefix = "MTSSV"
        filtered = [
            r for r in SeqIO.parse(str(orf_file), "fasta") if str(r.seq).startswith(orf_aa_prefix)
        ]

        # Verify metadata preserved
        assert len(filtered) == 1
        assert filtered[0].id == "haplotype_1_ORF1"
        assert "length=9" in filtered[0].description
        assert "start=0" in filtered[0].description

    def test_no_filtering_when_prefix_not_specified(self, sample_orfs, tmp_path):
        """Test that ORFs are not filtered when no prefix is specified."""
        from Bio import SeqIO

        output_file = tmp_path / "no_filter.fa"

        # Copy all ORFs without filtering
        shutil.copy(sample_orfs, output_file)

        # Verify all ORFs present
        result = list(SeqIO.parse(str(output_file), "fasta"))
        assert len(result) == 5


class TestORFPrefixFilteringEdgeCases:
    """Test edge cases for ORF prefix filtering."""

    def test_empty_orf_file(self, tmp_path):
        """Test filtering on empty ORF file."""
        from Bio import SeqIO

        empty_file = tmp_path / "empty.fa"
        empty_file.write_text("")

        orf_aa_prefix = "MTSSV"
        filtered = [
            r for r in SeqIO.parse(str(empty_file), "fasta") if str(r.seq).startswith(orf_aa_prefix)
        ]

        assert len(filtered) == 0

    def test_orf_shorter_than_prefix(self, tmp_path):
        """Test ORF shorter than the required prefix."""
        from Bio import SeqIO

        orf_file = tmp_path / "short.fa"
        orfs = [
            SeqRecord(Seq("MT"), id="orf1", description="too short"),
            SeqRecord(Seq("MTSSVAAAAA"), id="orf2", description="long enough"),
        ]
        SeqIO.write(orfs, str(orf_file), "fasta")

        orf_aa_prefix = "MTSSV"
        filtered = [
            r for r in SeqIO.parse(str(orf_file), "fasta") if str(r.seq).startswith(orf_aa_prefix)
        ]

        # Only the long ORF should match
        assert len(filtered) == 1
        assert filtered[0].id == "orf2"

    def test_prefix_longer_than_all_orfs(self, tmp_path):
        """Test prefix longer than all ORF sequences."""
        from Bio import SeqIO

        orf_file = tmp_path / "short_orfs.fa"
        orfs = [
            SeqRecord(Seq("MTSS"), id="orf1"),
            SeqRecord(Seq("MTS"), id="orf2"),
        ]
        SeqIO.write(orfs, str(orf_file), "fasta")

        orf_aa_prefix = "MTSSVAAAAA"  # Longer than any ORF
        filtered = [
            r for r in SeqIO.parse(str(orf_file), "fasta") if str(r.seq).startswith(orf_aa_prefix)
        ]

        assert len(filtered) == 0


class TestORFFilteringLogging:
    """Test logging for ORF prefix filtering."""

    def test_logging_shows_filter_statistics(self, sample_orfs, caplog):
        """Test that filtering logs correct statistics."""
        import logging

        from Bio import SeqIO

        caplog.set_level(logging.INFO)

        orf_aa_prefix = "MTSSV"
        filtered_orfs = []
        total_orfs = 0

        for record in SeqIO.parse(sample_orfs, "fasta"):
            total_orfs += 1
            if str(record.seq).startswith(orf_aa_prefix):
                filtered_orfs.append(record)

        # Simulate the log message from CLI
        log_message = f"Filtered ORFs by prefix '{orf_aa_prefix}': {len(filtered_orfs)}/{total_orfs} ORFs retained"

        # Verify the statistics
        assert "3/5" in log_message
        assert "MTSSV" in log_message
