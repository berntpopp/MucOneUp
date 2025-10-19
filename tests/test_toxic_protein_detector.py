"""Tests for toxic protein detector module.

Tests verify the toxic protein detection algorithm including:
1. Sliding window repeat analysis
2. Amino acid composition calculation
3. Composition similarity scoring
4. Overall toxic protein detection
5. FASTA file scanning
"""

import pytest

from muc_one_up.toxic_protein_detector import (
    composition_similarity,
    compute_amino_acid_composition,
    detect_toxic_protein_in_sequence,
    scan_orf_fasta,
    sliding_window_repeat_analysis,
)


class TestSlidingWindowRepeatAnalysis:
    """Test sliding window repeat detection."""

    def test_detects_perfect_repeats(self):
        """Test detection of perfect consensus repeats."""
        consensus = "ABCD"
        # Three perfect repeats
        region = "ABCDABCDABCD"

        count, avg_identity = sliding_window_repeat_analysis(
            region, consensus, identity_threshold=1.0
        )

        # Should detect 3 repeats with 100% identity
        assert count == 3
        assert avg_identity == 1.0

    def test_detects_partial_matches(self):
        """Test detection with partial matches above threshold."""
        consensus = "ABCD"
        # Partial match: 3/4 = 75% identity
        region = "ABXD"  # One mismatch

        count, avg_identity = sliding_window_repeat_analysis(
            region, consensus, identity_threshold=0.7
        )

        # Should detect 1 repeat with 75% identity
        assert count == 1
        assert avg_identity == 0.75

    def test_ignores_below_threshold(self):
        """Test that matches below threshold are ignored."""
        consensus = "ABCD"
        # Poor match: 2/4 = 50% identity
        region = "AXCX"

        count, avg_identity = sliding_window_repeat_analysis(
            region, consensus, identity_threshold=0.8
        )

        # Should detect 0 repeats
        assert count == 0
        assert avg_identity == 0.0

    def test_multiple_windows(self):
        """Test that sliding window scans all positions."""
        consensus = "ABC"
        # ABC appears at position 0 and position 3
        region = "ABCXABC"

        count, avg_identity = sliding_window_repeat_analysis(
            region, consensus, identity_threshold=1.0
        )

        # Should detect 2 matches
        assert count == 2
        assert avg_identity == 1.0

    def test_empty_region(self):
        """Test handling of empty region."""
        consensus = "ABCD"
        region = ""

        count, avg_identity = sliding_window_repeat_analysis(region, consensus)

        assert count == 0
        assert avg_identity == 0.0


class TestComputeAminoAcidComposition:
    """Test amino acid composition calculation."""

    def test_computes_frequencies(self):
        """Test basic frequency calculation."""
        seq = "AAABBC"  # 3 A's, 2 B's, 1 C = 6 total
        residues = ["A", "B", "C"]

        freq = compute_amino_acid_composition(seq, residues)

        assert freq["A"] == pytest.approx(3 / 6)
        assert freq["B"] == pytest.approx(2 / 6)
        assert freq["C"] == pytest.approx(1 / 6)

    def test_handles_missing_residues(self):
        """Test that missing residues have zero frequency."""
        seq = "AAA"
        residues = ["A", "B", "C"]

        freq = compute_amino_acid_composition(seq, residues)

        assert freq["A"] == 1.0
        assert freq["B"] == 0.0
        assert freq["C"] == 0.0

    def test_empty_sequence(self):
        """Test handling of empty sequence."""
        seq = ""
        residues = ["A", "B"]

        freq = compute_amino_acid_composition(seq, residues)

        assert freq["A"] == 0.0
        assert freq["B"] == 0.0


class TestCompositionSimilarity:
    """Test composition similarity scoring."""

    def test_identical_composition(self):
        """Test that identical sequences have similarity of 1.0."""
        mutant = "AABBCC"
        wildtype = "AABBCC"
        residues = ["A", "B", "C"]

        similarity = composition_similarity(mutant, wildtype, residues)

        assert similarity == pytest.approx(1.0)

    def test_different_composition(self):
        """Test that different compositions have lower similarity."""
        mutant = "AAAAAA"  # All A's
        wildtype = "AABBCC"  # Mixed
        residues = ["A", "B", "C"]

        similarity = composition_similarity(mutant, wildtype, residues)

        # Should be < 1.0 due to differences
        assert similarity < 1.0

    def test_similarity_score_calculation(self):
        """Test the similarity score calculation formula."""
        # Mutant: 4 A's, 2 B's = 6 total → A=4/6, B=2/6, C=0
        mutant = "AAAABB"
        # Wildtype: 3 A's, 3 B's = 6 total → A=3/6, B=3/6, C=0
        wildtype = "AAABBB"
        residues = ["A", "B", "C"]

        similarity = composition_similarity(mutant, wildtype, residues)

        # Diff = |4/6 - 3/6| + |2/6 - 3/6| + |0 - 0| = 1/6 + 1/6 = 2/6
        # Norm = 3/6 + 3/6 + 0 = 6/6 = 1.0
        # Score = 1 - 2/6 / 1.0 = 1 - 1/3 = 2/3
        assert similarity == pytest.approx(2 / 3)


class TestDetectToxicProteinInSequence:
    """Test main toxic protein detection function."""

    def test_detects_toxic_protein_high_score(self):
        """Test detection when overall score exceeds cutoff."""
        # Create a sequence with many repeats (high repeat score)
        consensus = "ABCD"
        protein_seq = consensus * 20  # Many perfect repeats

        result = detect_toxic_protein_in_sequence(
            protein_seq,
            consensus=consensus,
            expected_repeat_count=10,
            toxic_detection_cutoff=0.5,
        )

        # Should have high scores and be flagged as toxic
        assert result["repeat_count"] >= 10
        assert result["avg_repeat_identity"] > 0.8
        assert result["overall_score"] > 0.5
        assert result["toxic_flag"] == 1.0

    def test_detects_normal_protein_low_score(self):
        """Test detection when overall score is below cutoff."""
        # Create a sequence with few/no repeats
        protein_seq = "XYZXYZXYZXYZ"  # Different from consensus
        consensus = "ABCD"

        result = detect_toxic_protein_in_sequence(
            protein_seq,
            consensus=consensus,
            expected_repeat_count=10,
            toxic_detection_cutoff=0.5,
        )

        # Should have low overall score and NOT be flagged as toxic
        assert result["overall_score"] < 0.5
        assert result["toxic_flag"] == 0.0

    def test_removes_constant_regions(self):
        """Test that constant flanks are removed before analysis."""
        left_const = "LLLLL"
        right_const = "RRRRR"
        consensus = "ABCD"
        variable = consensus * 5

        # Full sequence with flanks
        protein_seq = left_const + variable + right_const

        result = detect_toxic_protein_in_sequence(
            protein_seq,
            left_const=left_const,
            right_const=right_const,
            consensus=consensus,
            identity_threshold=0.8,
        )

        # Should detect repeats only in variable region
        assert result["repeat_count"] == 5

    def test_uses_default_key_residues(self):
        """Test that default key residues [R, C, H] are used."""
        protein_seq = "RRRCCCHHHRRRCCCHHH"  # Rich in R, C, H

        result = detect_toxic_protein_in_sequence(protein_seq, consensus="RCHRCH")

        # Should compute composition similarity using R, C, H
        assert "composition_similarity" in result
        assert result["composition_similarity"] >= 0.0

    def test_weights_combination(self):
        """Test that repeat and composition weights are combined correctly."""
        protein_seq = "ABCDABCDABCD"
        consensus = "ABCD"

        result = detect_toxic_protein_in_sequence(
            protein_seq,
            consensus=consensus,
            w_repeat=0.6,
            w_composition=0.4,
            expected_repeat_count=3,
        )

        # Overall score should be weighted combination
        expected = (0.6 * result["repeat_score"]) + (0.4 * result["composition_similarity"])
        assert result["overall_score"] == pytest.approx(expected)


class TestScanOrfFasta:
    """Test FASTA file scanning."""

    def test_scans_single_orf(self, tmp_path):
        """Test scanning a FASTA file with one ORF."""
        # Create test FASTA file
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">ORF1\nABCDABCDABCD\n")

        results = scan_orf_fasta(str(fasta_file), consensus="ABCD")

        # Should have one result
        assert "ORF1" in results
        assert results["ORF1"]["repeat_count"] >= 1

    def test_scans_multiple_orfs(self, tmp_path):
        """Test scanning a FASTA file with multiple ORFs."""
        # Create test FASTA file with 3 ORFs
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">ORF1\nABCDABCD\n>ORF2\nXYZXYZXYZ\n>ORF3\nABCDABCDABCDABCD\n")

        results = scan_orf_fasta(str(fasta_file), consensus="ABCD")

        # Should have three results
        assert len(results) == 3
        assert "ORF1" in results
        assert "ORF2" in results
        assert "ORF3" in results

        # ORF1 and ORF3 should have repeats, ORF2 should not
        assert results["ORF1"]["repeat_count"] >= 1
        assert results["ORF3"]["repeat_count"] >= 2

    def test_handles_multiline_sequences(self, tmp_path):
        """Test that multiline sequences are correctly concatenated."""
        # Create test FASTA with sequence split across lines
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">ORF1\nABCD\nABCD\nABCD\n")

        results = scan_orf_fasta(str(fasta_file), consensus="ABCD")

        # Should concatenate lines and detect repeats
        assert "ORF1" in results
        assert results["ORF1"]["repeat_count"] == 3

    def test_passes_detection_kwargs(self, tmp_path):
        """Test that detection kwargs are passed through."""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">ORF1\nABCDABCD\n")

        results = scan_orf_fasta(
            str(fasta_file),
            consensus="ABCD",
            identity_threshold=0.9,
            toxic_detection_cutoff=0.3,
        )

        # Should use custom parameters (toxic_flag depends on cutoff)
        assert "ORF1" in results
        assert "toxic_flag" in results["ORF1"]

    def test_handles_empty_fasta(self, tmp_path):
        """Test handling of empty FASTA file."""
        fasta_file = tmp_path / "empty.fasta"
        fasta_file.write_text("")

        results = scan_orf_fasta(str(fasta_file))

        # Should return empty dict
        assert results == {}
