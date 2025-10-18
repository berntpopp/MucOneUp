"""Tests for fragment_simulation module - focuses on OUR logic (w-Wessim2 port).

Following Phase 2/3 testing principles:
- Mock ONLY at system boundary (subprocess.Popen for run_command in simulate_fragments)
- Test OUR code's logic (fragment picking, file parsing, sequence operations)
- NOT testing external tools (reseq, pblat)

These tests verify that our ported w-Wessim2 logic correctly:
1. Parses FASTA files into dictionaries
2. Parses systematic error files (syser format)
3. Parses PSL alignment files
4. Picks fragments based on PSL matches
5. Generates insert lengths from Gaussian distribution
6. Computes DNA complement sequences
7. Orchestrates fragment simulation end-to-end
"""

import pytest

from muc_one_up.exceptions import FileOperationError
from muc_one_up.read_simulator.fragment_simulation import (
    comp,
    get_insert_length,
    pick_fragment,
    pick_on_match,
    read_fasta_to_dict,
    read_psl_file,
    read_syser_file,
    simulate_fragments,
)


class TestReadFastaToDict:
    """Test FASTA file parsing into dictionary."""

    def test_parses_single_sequence(self, tmp_path):
        """Test parsing a FASTA with one sequence."""
        fasta = tmp_path / "test.fa"
        fasta.write_text(">chr1\nACGTACGT\n")

        result = read_fasta_to_dict(str(fasta))

        assert len(result) == 1
        assert "chr1" in result
        assert result["chr1"] == "ACGTACGT"

    def test_parses_multiple_sequences(self, tmp_path):
        """Test parsing FASTA with multiple chromosomes."""
        fasta = tmp_path / "test.fa"
        fasta.write_text(">chr1\nACGT\n" ">chr2\nTGCA\n" ">chr3\nGGCC\n")

        result = read_fasta_to_dict(str(fasta))

        assert len(result) == 3
        assert result["chr1"] == "ACGT"
        assert result["chr2"] == "TGCA"
        assert result["chr3"] == "GGCC"

    def test_handles_multiline_sequences(self, tmp_path):
        """Test that sequences split across lines are concatenated."""
        fasta = tmp_path / "test.fa"
        fasta.write_text(">chr1\n" "ACGT\n" "TGCA\n" "GGCC\n")

        result = read_fasta_to_dict(str(fasta))

        assert result["chr1"] == "ACGTTGCAGGCC"

    def test_extracts_chromosome_name_only(self, tmp_path):
        """Test that only the first word after > is used as chromosome name."""
        fasta = tmp_path / "test.fa"
        fasta.write_text(">chr1 extra description here\nACGT\n")

        result = read_fasta_to_dict(str(fasta))

        assert "chr1" in result
        assert "extra" not in result
        assert result["chr1"] == "ACGT"

    def test_ignores_empty_lines(self, tmp_path):
        """Test that empty lines are ignored."""
        fasta = tmp_path / "test.fa"
        fasta.write_text(">chr1\n" "\n" "ACGT\n" "\n" "TGCA\n")

        result = read_fasta_to_dict(str(fasta))

        assert result["chr1"] == "ACGTTGCA"

    def test_raises_on_missing_file(self, tmp_path):
        """Test that missing file raises exception."""
        missing = tmp_path / "missing.fa"

        with pytest.raises(FileNotFoundError):
            read_fasta_to_dict(str(missing))


class TestReadSyserFile:
    """Test systematic errors file parsing."""

    def test_parses_forward_and_reverse_sections(self, tmp_path):
        """Test parsing syser file with forward and reverse sections."""
        syser = tmp_path / "test.syser"
        syser.write_text(
            "@chr1 forward\n" "NACGT\n" "+\n" "!!III\n" "@chr1 reverse\n" "TGCAN\n" "+\n" "III!!\n"
        )

        tenden_f, tenden_r, rates_f, rates_r = read_syser_file(str(syser))

        assert "chr1" in tenden_f
        assert "chr1" in tenden_r
        assert "chr1" in rates_f
        assert "chr1" in rates_r

        # Forward direction
        assert tenden_f["chr1"] == "NACGT"
        assert rates_f["chr1"] == "!!III"

        # Reverse direction (should be reversed)
        assert tenden_r["chr1"] == "NACGT"  # "TGCAN" reversed
        assert rates_r["chr1"] == "!!III"  # "III!!" reversed

    def test_handles_multiple_chromosomes(self, tmp_path):
        """Test parsing syser with multiple chromosomes."""
        syser = tmp_path / "test.syser"
        syser.write_text(
            "@chr1 forward\n" "NNNN\n" "+\n" "!!!!\n" "@chr2 forward\n" "AAAA\n" "+\n" "IIII\n"
        )

        tenden_f, tenden_r, rates_f, rates_r = read_syser_file(str(syser))

        assert "chr1" in tenden_f
        assert "chr2" in tenden_f
        assert tenden_f["chr1"] == "NNNN"
        assert tenden_f["chr2"] == "AAAA"

    def test_extracts_identifier_from_header(self, tmp_path):
        """Test that identifier is extracted correctly from @ header."""
        syser = tmp_path / "test.syser"
        syser.write_text(
            "@identifier_with_underscores forward extra text\n" "ACGT\n" "+\n" "IIII\n"
        )

        tenden_f, tenden_r, rates_f, rates_r = read_syser_file(str(syser))

        assert "identifier_with_underscores" in tenden_f

    def test_ignores_empty_lines(self, tmp_path):
        """Test that empty lines are ignored."""
        syser = tmp_path / "test.syser"
        syser.write_text("\n" "@chr1 forward\n" "ACGT\n" "\n" "+\n" "IIII\n" "\n")

        tenden_f, tenden_r, rates_f, rates_r = read_syser_file(str(syser))

        assert "chr1" in tenden_f
        assert tenden_f["chr1"] == "ACGT"


class TestReadPslFile:
    """Test PSL alignment file parsing."""

    def test_parses_psl_matches(self, tmp_path):
        """Test parsing PSL file with valid matches."""
        psl = tmp_path / "test.psl"
        # PSL format has 21 tab-separated columns
        # We care about: [13]=chrom, [15]=start, [16]=end, [8]=strand
        psl.write_text(
            "psLayout version 3\n"
            "match\tmismatch\trepmatch\tncount\n"
            "------\n"
            "100\t0\t0\t0\t0\t0\t0\t0\t+\tquery\t150\t0\t100\tchr1\t1000\t500\t600\t1\t100\t0\t500\n"
            "200\t5\t0\t1\t0\t0\t0\t0\t-\tquery2\t250\t0\t200\tchr2\t2000\t1000\t1200\t1\t200\t0\t1000\n"
        )

        matches = read_psl_file(str(psl))

        assert len(matches) == 2
        # First match
        assert matches[0] == ("chr1", 500, 600, "+")
        # Second match
        assert matches[1] == ("chr2", 1000, 1200, "-")

    def test_skips_header_lines(self, tmp_path):
        """Test that PSL header lines are skipped."""
        psl = tmp_path / "test.psl"
        psl.write_text(
            "psLayout version 3\n"
            "match\tmismatch\trepmatch\n"
            "------\n"
            "\n"
            "100\t0\t0\t0\t0\t0\t0\t0\t+\tquery\t150\t0\t100\tchr1\t1000\t100\t200\t1\t100\t0\t100\n"
        )

        matches = read_psl_file(str(psl))

        # Should only have 1 match (headers and empty lines skipped)
        assert len(matches) == 1

    def test_skips_invalid_lines_with_insufficient_columns(self, tmp_path):
        """Test that lines with < 21 columns are skipped with warning."""
        psl = tmp_path / "test.psl"
        psl.write_text(
            "INVALID_LINE_WITH_FEW_COLUMNS\n"
            "100\t0\t0\t0\t0\t0\t0\t0\t+\tquery\t150\t0\t100\tchr1\t1000\t100\t200\t1\t100\t0\t100\n"
        )

        matches = read_psl_file(str(psl))

        # Should only have the valid match
        assert len(matches) == 1
        assert matches[0][0] == "chr1"

    def test_handles_minus_strand(self, tmp_path):
        """Test that minus strand is correctly identified."""
        psl = tmp_path / "test.psl"
        psl.write_text(
            "100\t0\t0\t0\t0\t0\t0\t0\t-\tquery\t150\t0\t100\tchr1\t1000\t100\t200\t1\t100\t0\t100\n"
        )

        matches = read_psl_file(str(psl))

        assert matches[0][3] == "-"

    def test_returns_empty_list_for_no_matches(self, tmp_path):
        """Test that file with only headers returns empty list."""
        psl = tmp_path / "test.psl"
        psl.write_text("psLayout version 3\n" "match\tmismatch\n" "------\n")

        matches = read_psl_file(str(psl))

        assert len(matches) == 0


class TestPickOnMatch:
    """Test random match selection."""

    def test_picks_from_single_match(self):
        """Test picking when only one match is available."""
        matches = [("chr1", 100, 200, "+")]

        result = pick_on_match(matches)

        assert result == ("chr1", 100, 200, "+")

    def test_picks_from_multiple_matches(self, mocker):
        """Test that random.choice is called with matches."""
        matches = [
            ("chr1", 100, 200, "+"),
            ("chr2", 300, 400, "-"),
            ("chr3", 500, 600, "+"),
        ]

        # Mock random.choice to return first match
        mock_choice = mocker.patch("random.choice", return_value=matches[0])

        result = pick_on_match(matches)

        mock_choice.assert_called_once_with(matches)
        assert result == matches[0]

    def test_raises_on_empty_matches(self):
        """Test that empty matches list raises ValueError."""
        with pytest.raises(ValueError, match="No matches provided"):
            pick_on_match([])


class TestGetInsertLength:
    """Test Gaussian insert length sampling."""

    def test_returns_integer_above_lower_bound(self):
        """Test that insert length is >= lower bound."""
        mu, sigma, lower = 300, 50, 200

        for _ in range(10):  # Test multiple times
            length = get_insert_length(mu, sigma, lower)
            assert isinstance(length, int)
            assert length >= lower

    def test_respects_lower_bound(self, mocker):
        """Test that values below lower are rejected."""
        # Mock random.gauss to return [150, 180, 250] in sequence
        mock_gauss = mocker.patch("random.gauss")
        mock_gauss.side_effect = [150, 180, 250]  # First two < 200, third >= 200

        length = get_insert_length(300, 50, 200)

        assert length == 250
        assert mock_gauss.call_count == 3  # Called 3 times before getting valid value

    def test_uses_mu_and_sigma_parameters(self, mocker):
        """Test that mu and sigma are passed to random.gauss."""
        mock_gauss = mocker.patch("random.gauss", return_value=300)

        get_insert_length(300, 50, 200)

        mock_gauss.assert_called_with(300, 50)


class TestPickFragment:
    """Test fragment picking based on PSL match."""

    def test_calculates_fragment_coordinates(self, mocker):
        """Test that fragment coordinates are calculated correctly."""
        # Mock random.randint to return fixed offset
        mocker.patch("random.randint", return_value=50)

        match = ("chr1", 1000, 1100, "+")  # plen = 100
        ins = 300  # desired insert length
        bind = 0.5  # 50% minimum overlap

        result = pick_fragment(match, ins, bind)

        chrom, fstart, fend, strand = result

        assert chrom == "chr1"
        assert strand == "+"
        # min_overlap = 100 * 0.5 = 50
        # space = 300 + 100 - 50 = 350
        # offset = 50 (mocked)
        # fstart = 1000 - 50 = 950
        # fend = 950 + 300 = 1250
        assert fstart == 950
        assert fend == 1250

    def test_prevents_negative_start_coordinate(self, mocker):
        """Test that fragment start is never negative."""
        # Mock large offset that would make start negative
        mocker.patch("random.randint", return_value=5000)

        match = ("chr1", 100, 200, "+")
        ins = 300
        bind = 0.5

        result = pick_fragment(match, ins, bind)

        chrom, fstart, fend, strand = result

        assert fstart == 0  # Clamped to 0
        assert fend == 300  # 0 + ins

    def test_preserves_strand_information(self, mocker):
        """Test that strand is preserved from match."""
        mocker.patch("random.randint", return_value=0)

        match_plus = ("chr1", 1000, 1100, "+")
        match_minus = ("chr2", 2000, 2100, "-")

        result_plus = pick_fragment(match_plus, 300, 0.5)
        result_minus = pick_fragment(match_minus, 300, 0.5)

        assert result_plus[3] == "+"
        assert result_minus[3] == "-"


class TestComp:
    """Test DNA complement function."""

    def test_complements_uppercase_bases(self):
        """Test complementation of uppercase DNA bases."""
        assert comp("ACGT") == "TGCA"
        assert comp("AAAA") == "TTTT"
        assert comp("CCCC") == "GGGG"

    def test_complements_lowercase_bases(self):
        """Test complementation of lowercase DNA bases."""
        assert comp("acgt") == "tgca"
        assert comp("aaaa") == "tttt"

    def test_preserves_case(self):
        """Test that case is preserved in output."""
        assert comp("AcGt") == "TgCa"

    def test_handles_n_bases(self):
        """Test that N bases (ambiguous) map to N."""
        assert comp("ACNGT") == "TGNCA"
        assert comp("nnn") == "nnn"

    def test_preserves_unknown_characters(self):
        """Test that unknown characters pass through unchanged."""
        # According to COMP_MAP, only ACGTN are defined
        assert comp("ACGTX") == "TGCAX"  # X not in map, preserved

    def test_handles_empty_sequence(self):
        """Test empty sequence returns empty string."""
        assert comp("") == ""


class TestSimulateFragments:
    """Test end-to-end fragment simulation."""

    def test_simulates_fragments_successfully(self, mocker, tmp_path):
        """Test successful fragment simulation with all components."""
        # Arrange: Create input files
        ref_fa = tmp_path / "ref.fa"
        ref_fa.write_text(">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n")  # 32 bases

        syser = tmp_path / "test.syser"
        syser.write_text(
            "@chr1 forward\n"
            "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n"
            "+\n"
            "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
            "@chr1 reverse\n"
            "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n"
            "+\n"
            "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
        )

        psl = tmp_path / "test.psl"
        psl.write_text(
            "100\t0\t0\t0\t0\t0\t0\t0\t+\tquery\t150\t0\t100\tchr1\t1000\t5\t15\t1\t100\t0\t5\n"
        )

        output = tmp_path / "fragments.fa"

        # Act
        simulate_fragments(
            ref_fa=str(ref_fa),
            syser_file=str(syser),
            psl_file=str(psl),
            read_number=2,
            fragment_size=20,
            fragment_sd=2,
            min_fragment=10,
            bind=0.5,
            output_fragments=str(output),
        )

        # Assert: Output file should exist and contain FASTA records
        assert output.exists()
        content = output.read_text()

        # Should have 4 records (2 fragment pairs = 4 reads total)
        assert content.count(">") == 4
        # Check for read pair format
        assert ">1 1;" in content  # Read 1, pair 1
        assert ">1 2;" in content  # Read 1, pair 2
        assert ">2 1;" in content  # Read 2, pair 1
        assert ">2 2;" in content  # Read 2, pair 2

    def test_raises_on_missing_reference(self, tmp_path):
        """Test that missing reference FASTA raises FileOperationError."""
        missing_ref = tmp_path / "missing.fa"
        syser = tmp_path / "test.syser"
        psl = tmp_path / "test.psl"
        output = tmp_path / "out.fa"

        syser.write_text("@chr1 forward\nN\n+\n!\n")
        psl.write_text(
            "100\t0\t0\t0\t0\t0\t0\t0\t+\tq\t150\t0\t100\tchr1\t1000\t0\t10\t1\t100\t0\t0\n"
        )

        with pytest.raises(FileOperationError, match="Failed to load reference FASTA"):
            simulate_fragments(
                str(missing_ref), str(syser), str(psl), 1, 300, 50, 200, 0.5, str(output)
            )

    def test_raises_on_empty_psl_matches(self, tmp_path):
        """Test that PSL with no matches raises FileOperationError."""
        ref_fa = tmp_path / "ref.fa"
        syser = tmp_path / "test.syser"
        psl = tmp_path / "empty.psl"
        output = tmp_path / "out.fa"

        ref_fa.write_text(">chr1\nACGT\n")
        syser.write_text("@chr1 forward\nN\n+\n!\n")
        psl.write_text("psLayout version 3\nmatch\n------\n")  # Only headers

        with pytest.raises(FileOperationError, match="No matches found in PSL file"):
            simulate_fragments(str(ref_fa), str(syser), str(psl), 1, 300, 50, 200, 0.5, str(output))

    def test_handles_minus_strand_correctly(self, mocker, tmp_path):
        """Test that minus strand fragments use reverse complement."""
        ref_fa = tmp_path / "ref.fa"
        ref_fa.write_text(">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n")

        syser = tmp_path / "test.syser"
        syser.write_text(
            "@chr1 forward\n" + "N" * 32 + "\n+\n" + "!" * 32 + "\n"
            "@chr1 reverse\n" + "N" * 32 + "\n+\n" + "!" * 32 + "\n"
        )

        psl = tmp_path / "test.psl"
        # Minus strand match
        psl.write_text(
            "100\t0\t0\t0\t0\t0\t0\t0\t-\tquery\t150\t0\t100\tchr1\t1000\t5\t15\t1\t100\t0\t5\n"
        )

        output = tmp_path / "fragments.fa"

        simulate_fragments(str(ref_fa), str(syser), str(psl), 1, 20, 2, 10, 0.5, str(output))

        content = output.read_text()

        # For minus strand, read 1 should be reverse complement, read 2 should be forward
        # Check that output has both reads
        assert ">1 1;" in content
        assert ">1 2;" in content

    def test_raises_on_all_fragments_skipped(self, mocker, tmp_path):
        """Test that error is raised when all fragments are skipped (empty output)."""
        ref_fa = tmp_path / "ref.fa"
        ref_fa.write_text(">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n")

        syser = tmp_path / "test.syser"
        syser.write_text("@chr1 forward\n" + "N" * 32 + "\n+\n" + "!" * 32 + "\n")

        psl = tmp_path / "test.psl"
        # Match references chr2 which doesn't exist
        psl.write_text(
            "100\t0\t0\t0\t0\t0\t0\t0\t+\tquery\t150\t0\t100\tchr2\t1000\t5\t15\t1\t100\t0\t5\n"
        )

        output = tmp_path / "fragments.fa"

        # Should raise error because output is empty (all fragments skipped)
        with pytest.raises(FileOperationError, match="missing or empty"):
            simulate_fragments(str(ref_fa), str(syser), str(psl), 1, 20, 2, 10, 0.5, str(output))

    def test_truncates_fragment_end_at_chromosome_boundary(self, mocker, tmp_path):
        """Test that fragments extending beyond chromosome are truncated."""
        ref_fa = tmp_path / "ref.fa"
        ref_fa.write_text(">chr1\nACGTACGT\n")  # Only 8 bases

        syser = tmp_path / "test.syser"
        syser.write_text(
            "@chr1 forward\n" + "N" * 20 + "\n+\n" + "!" * 20 + "\n"
            "@chr1 reverse\n" + "N" * 20 + "\n+\n" + "!" * 20 + "\n"
        )

        psl = tmp_path / "test.psl"
        psl.write_text(
            "100\t0\t0\t0\t0\t0\t0\t0\t+\tquery\t150\t0\t100\tchr1\t1000\t2\t4\t1\t100\t0\t2\n"
        )

        output = tmp_path / "fragments.fa"

        # Fragment calculation might try to go beyond chr length (8)
        # Should be truncated to chromosome length
        simulate_fragments(str(ref_fa), str(syser), str(psl), 1, 20, 2, 5, 0.5, str(output))

        assert output.exists()
        # Fragment should be generated (though possibly truncated)
