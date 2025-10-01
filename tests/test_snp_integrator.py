"""Tests for muc_one_up.snp_integrator module.

Tests cover:
- SNP file parsing
- VNTR boundary calculation
- Random SNP generation
- SNP application to sequences
- SNP file writing
"""

from pathlib import Path

import pytest

from muc_one_up.snp_integrator import (
    DNA_BASES,
    apply_snps_to_sequences,
    generate_random_snps,
    get_vntr_boundaries,
    parse_snp_file,
    write_snps_to_file,
)


@pytest.mark.unit
class TestParseSnpFile:
    """Tests for parse_snp_file function."""

    def test_parse_valid_snp_file(self, sample_snp_file: Path):
        """Test parsing valid SNP file."""
        snps = parse_snp_file(str(sample_snp_file))

        assert len(snps) == 3
        # Check first SNP
        assert snps[0]["haplotype"] == 1
        assert snps[0]["position"] == 110
        assert snps[0]["ref_base"] == "A"
        assert snps[0]["alt_base"] == "G"

    def test_parse_snp_file_with_comments(self, tmp_path: Path):
        """Test parsing SNP file with comment lines."""
        snp_file = tmp_path / "snps.tsv"
        snp_file.write_text("# This is a comment\n1\t100\tA\tG\n# Another comment\n2\t200\tC\tT\n")

        snps = parse_snp_file(str(snp_file))

        assert len(snps) == 2
        assert snps[0]["position"] == 100
        assert snps[1]["position"] == 200

    def test_parse_snp_file_with_empty_lines(self, tmp_path: Path):
        """Test parsing SNP file with empty lines."""
        snp_file = tmp_path / "snps.tsv"
        snp_file.write_text("\n1\t100\tA\tG\n\n\n2\t200\tC\tT\n\n")

        snps = parse_snp_file(str(snp_file))

        assert len(snps) == 2

    def test_parse_snp_file_lowercase_bases(self, tmp_path: Path):
        """Test parsing SNP file with lowercase bases (should be uppercased)."""
        snp_file = tmp_path / "snps.tsv"
        snp_file.write_text("1\t100\ta\tg\n2\t200\tc\tt\n")

        snps = parse_snp_file(str(snp_file))

        assert snps[0]["ref_base"] == "A"
        assert snps[0]["alt_base"] == "G"
        assert snps[1]["ref_base"] == "C"
        assert snps[1]["alt_base"] == "T"

    def test_parse_snp_file_missing_file(self):
        """Test parsing nonexistent file raises ValueError."""
        with pytest.raises(ValueError) as exc_info:
            parse_snp_file("/nonexistent/file.tsv")

        assert "SNP file not found" in str(exc_info.value)

    def test_parse_snp_file_invalid_field_count(self, tmp_path: Path):
        """Test parsing file with wrong number of fields raises ValueError."""
        snp_file = tmp_path / "snps.tsv"
        snp_file.write_text("1\t100\tA\n")  # Missing alt_base field

        with pytest.raises(ValueError) as exc_info:
            parse_snp_file(str(snp_file))

        assert "Expected 4 tab-separated fields" in str(exc_info.value)

    def test_parse_snp_file_invalid_haplotype(self, tmp_path: Path):
        """Test parsing file with invalid haplotype index raises ValueError."""
        snp_file = tmp_path / "snps.tsv"
        snp_file.write_text("3\t100\tA\tG\n")  # Haplotype 3 (only 1 or 2 allowed)

        with pytest.raises(ValueError) as exc_info:
            parse_snp_file(str(snp_file))

        assert "Haplotype index must be 1 or 2" in str(exc_info.value)

    def test_parse_snp_file_negative_position(self, tmp_path: Path):
        """Test parsing file with negative position raises ValueError."""
        snp_file = tmp_path / "snps.tsv"
        snp_file.write_text("1\t-10\tA\tG\n")

        with pytest.raises(ValueError) as exc_info:
            parse_snp_file(str(snp_file))

        assert "Position must be a non-negative integer" in str(exc_info.value)

    def test_parse_snp_file_invalid_ref_base(self, tmp_path: Path):
        """Test parsing file with invalid reference base raises ValueError."""
        snp_file = tmp_path / "snps.tsv"
        snp_file.write_text("1\t100\tX\tG\n")  # X is not a valid DNA base

        with pytest.raises(ValueError) as exc_info:
            parse_snp_file(str(snp_file))

        assert "Reference base must be one of" in str(exc_info.value)

    def test_parse_snp_file_invalid_alt_base(self, tmp_path: Path):
        """Test parsing file with invalid alternative base raises ValueError."""
        snp_file = tmp_path / "snps.tsv"
        snp_file.write_text("1\t100\tA\tZ\n")  # Z is not a valid DNA base

        with pytest.raises(ValueError) as exc_info:
            parse_snp_file(str(snp_file))

        assert "Alternative base must be one of" in str(exc_info.value)

    def test_parse_snp_file_same_ref_and_alt(self, tmp_path: Path):
        """Test parsing file where ref and alt are the same (should warn)."""
        snp_file = tmp_path / "snps.tsv"
        snp_file.write_text("1\t100\tA\tA\n")

        # Should parse successfully but log a warning
        snps = parse_snp_file(str(snp_file))

        assert len(snps) == 1
        assert snps[0]["ref_base"] == snps[0]["alt_base"]

    def test_parse_snp_file_empty_file(self, tmp_path: Path):
        """Test parsing empty file returns empty list."""
        snp_file = tmp_path / "snps.tsv"
        snp_file.write_text("")

        snps = parse_snp_file(str(snp_file))

        assert snps == []


@pytest.mark.unit
class TestGetVntrBoundaries:
    """Tests for get_vntr_boundaries function."""

    def test_get_vntr_boundaries_basic(self, minimal_config: dict):
        """Test calculating VNTR boundaries for simple sequences."""
        # Create simple simulation results
        left_const = minimal_config["constants"]["hg38"]["left"]
        right_const = minimal_config["constants"]["hg38"]["right"]
        vntr_region = "ATCG" * 50  # 200 bp VNTR

        seq1 = left_const + vntr_region + right_const
        seq2 = left_const + vntr_region + right_const

        simulation_results = [(seq1, ["1", "2", "X"]), (seq2, ["1", "2", "A"])]

        boundaries = get_vntr_boundaries(simulation_results, minimal_config)

        assert len(boundaries) == 2
        assert boundaries[0]["start"] == len(left_const)
        assert boundaries[0]["end"] == len(seq1) - len(right_const)
        assert boundaries[1]["start"] == len(left_const)
        assert boundaries[1]["end"] == len(seq2) - len(right_const)

    def test_get_vntr_boundaries_different_lengths(self, minimal_config: dict):
        """Test VNTR boundaries for haplotypes with different VNTR lengths."""
        left_const = minimal_config["constants"]["hg38"]["left"]
        right_const = minimal_config["constants"]["hg38"]["right"]

        # Different VNTR lengths
        seq1 = left_const + "ATCG" * 30 + right_const  # 120 bp VNTR
        seq2 = left_const + "ATCG" * 60 + right_const  # 240 bp VNTR

        simulation_results = [(seq1, ["1", "2"]), (seq2, ["1", "2"])]

        boundaries = get_vntr_boundaries(simulation_results, minimal_config)

        assert len(boundaries) == 2
        # Start should be same for both
        assert boundaries[0]["start"] == boundaries[1]["start"]
        # End should differ based on VNTR length
        assert boundaries[0]["end"] < boundaries[1]["end"]


@pytest.mark.unit
class TestGenerateRandomSnps:
    """Tests for generate_random_snps function."""

    def test_generate_random_snps_basic(self):
        """Test generating random SNPs with basic parameters."""
        sequences = ["ATCGATCGATCG" * 100]  # 1200 bp
        density = 1.0  # 1 SNP per kb

        snps = generate_random_snps(sequences, density)

        # Should generate approximately 1-2 SNPs (1.2 kb * 1 SNP/kb)
        assert len(snps) > 0
        assert all(snp["haplotype"] == 1 for snp in snps)
        assert all(isinstance(snp["position"], int) for snp in snps)
        assert all(snp["ref_base"] in DNA_BASES for snp in snps)
        assert all(snp["alt_base"] in DNA_BASES for snp in snps)
        # Ref and alt should be different
        assert all(snp["ref_base"] != snp["alt_base"] for snp in snps)

    def test_generate_random_snps_empty_sequences(self):
        """Test generating SNPs with empty sequence list."""
        snps = generate_random_snps([], density_per_kb=1.0)

        assert snps == []

    def test_generate_random_snps_zero_density(self):
        """Test generating SNPs with zero density."""
        sequences = ["ATCGATCGATCG" * 100]
        density = 0.0

        snps = generate_random_snps(sequences, density)

        assert snps == []

    def test_generate_random_snps_reproducible_with_seed(self):
        """Test that SNP generation is reproducible with seed."""
        import random

        sequences = ["ATCGATCGATCG" * 100]
        density = 2.0

        random.seed(42)
        snps1 = generate_random_snps(sequences, density)

        random.seed(42)
        snps2 = generate_random_snps(sequences, density)

        assert snps1 == snps2

    def test_generate_random_snps_invalid_region(self):
        """Test generating SNPs with invalid region raises ValueError."""
        sequences = ["ATCGATCGATCG"]

        with pytest.raises(ValueError) as exc_info:
            generate_random_snps(sequences, density_per_kb=1.0, region="invalid")

        assert "Invalid region parameter" in str(exc_info.value)


@pytest.mark.unit
class TestApplySnpsToSequences:
    """Tests for apply_snps_to_sequences function."""

    def test_apply_snps_basic(self):
        """Test applying SNPs to sequences."""
        sequences = ["ATCGATCGATCG", "GCTAGCTAGCTA"]
        snps = [
            {"haplotype": 1, "position": 0, "ref_base": "A", "alt_base": "G"},
            {"haplotype": 2, "position": 0, "ref_base": "G", "alt_base": "A"},
        ]

        modified_seqs, applied_snps = apply_snps_to_sequences(sequences, snps)

        assert len(modified_seqs) == 2
        # First sequence should have A→G at position 0
        assert modified_seqs[0][0] == "G"
        # Second sequence should have G→A at position 0
        assert modified_seqs[1][0] == "A"
        # Should track applied SNPs
        assert len(applied_snps[0]) == 1
        assert len(applied_snps[1]) == 1

    def test_apply_snps_no_snps(self):
        """Test applying empty SNP list leaves sequences unchanged."""
        sequences = ["ATCGATCGATCG", "GCTAGCTAGCTA"]
        snps = []

        modified_seqs, applied_snps = apply_snps_to_sequences(sequences, snps)

        assert modified_seqs == sequences
        assert len(applied_snps) == 2
        assert all(len(s) == 0 for s in applied_snps.values())

    def test_apply_snps_reference_mismatch(self):
        """Test applying SNP with mismatched reference base skips SNP."""
        sequences = ["ATCGATCGATCG"]
        snps = [
            {"haplotype": 1, "position": 0, "ref_base": "G", "alt_base": "C"}
            # Expected ref is A, but SNP says G
        ]

        modified_seqs, applied_snps = apply_snps_to_sequences(sequences, snps)

        # SNP should be skipped, sequence unchanged
        assert modified_seqs[0] == sequences[0]
        assert len(applied_snps[0]) == 0  # No SNPs applied

    def test_apply_snps_skip_reference_check(self):
        """Test applying SNPs with reference check disabled."""
        sequences = ["ATCGATCGATCG"]
        snps = [
            {"haplotype": 1, "position": 0, "ref_base": "G", "alt_base": "C"}
            # Ref mismatch, but we skip check
        ]

        modified_seqs, applied_snps = apply_snps_to_sequences(
            sequences, snps, skip_reference_check=True
        )

        # Should apply SNP despite ref mismatch
        assert modified_seqs[0][0] == "C"
        assert len(applied_snps[0]) == 1

    def test_apply_snps_position_out_of_range(self):
        """Test applying SNP at invalid position skips SNP."""
        sequences = ["ATCGATCGATCG"]  # 12 bp
        snps = [
            {"haplotype": 1, "position": 100, "ref_base": "A", "alt_base": "G"}
            # Position 100 is beyond sequence length
        ]

        modified_seqs, applied_snps = apply_snps_to_sequences(sequences, snps)

        # SNP should be skipped, sequence unchanged
        assert modified_seqs[0] == sequences[0]
        assert len(applied_snps[0]) == 0  # No SNPs applied

    def test_apply_snps_multiple_same_position(self):
        """Test applying multiple SNPs to same position."""
        sequences = ["ATCGATCGATCG"]
        snps = [
            {"haplotype": 1, "position": 0, "ref_base": "A", "alt_base": "G"},
            {"haplotype": 1, "position": 0, "ref_base": "G", "alt_base": "T"},
        ]

        # Second SNP ref should match first SNP alt
        modified_seqs, applied_snps = apply_snps_to_sequences(sequences, snps)

        # Final base should be T (last SNP applied)
        assert modified_seqs[0][0] == "T"
        assert len(applied_snps[0]) == 2


@pytest.mark.unit
class TestWriteSnpsToFile:
    """Tests for write_snps_to_file function."""

    def test_write_snps_to_file(self, tmp_path: Path):
        """Test writing SNPs to file."""
        output_file = tmp_path / "output_snps.tsv"
        snps = [
            {"haplotype": 1, "position": 100, "ref_base": "A", "alt_base": "G"},
            {"haplotype": 2, "position": 200, "ref_base": "C", "alt_base": "T"},
        ]

        write_snps_to_file(snps, str(output_file))

        assert output_file.exists()

        # Verify content
        with output_file.open() as f:
            lines = [line.strip() for line in f if not line.startswith("#")]

        assert len(lines) == 2
        assert "1\t100\tA\tG" in lines[0]
        assert "2\t200\tC\tT" in lines[1]

    def test_write_snps_empty_list(self, tmp_path: Path):
        """Test writing empty SNP list creates file with header only."""
        output_file = tmp_path / "output_snps.tsv"

        write_snps_to_file([], str(output_file))

        assert output_file.exists()

        with output_file.open() as f:
            lines = [line for line in f if not line.strip().startswith("#")]

        # Should have no data lines
        assert len([line for line in lines if line.strip()]) == 0


@pytest.mark.bioinformatics
class TestSnpIntegratorBioinformatics:
    """Bioinformatics-specific tests for SNP integration."""

    def test_snp_preserves_sequence_length(self):
        """Test that SNP application preserves sequence length."""
        sequences = ["ATCGATCGATCG" * 100]
        original_length = len(sequences[0])

        snps = [{"haplotype": 1, "position": 50, "ref_base": "A", "alt_base": "G"}]

        modified_seqs, _ = apply_snps_to_sequences(sequences, snps)

        assert len(modified_seqs[0]) == original_length

    def test_snp_changes_single_base_only(self):
        """Test that SNP modifies only the targeted base."""
        original_seq = "ATCGATCGATCG"
        sequences = [original_seq]
        snps = [{"haplotype": 1, "position": 5, "ref_base": "T", "alt_base": "C"}]

        modified_seqs, _ = apply_snps_to_sequences(sequences, snps)

        # Only position 5 should change
        for i, (orig, mod) in enumerate(zip(original_seq, modified_seqs[0], strict=True)):
            if i == 5:
                assert mod == "C"
            else:
                assert mod == orig

    def test_random_snps_no_duplicate_positions(self):
        """Test that random SNP generation doesn't create duplicates at same position."""
        sequences = ["ATCGATCGATCG" * 100]
        snps = generate_random_snps(sequences, density_per_kb=10.0)

        # Count positions per haplotype
        positions = [snp["position"] for snp in snps if snp["haplotype"] == 1]

        # Should have no duplicates
        assert len(positions) == len(set(positions))
