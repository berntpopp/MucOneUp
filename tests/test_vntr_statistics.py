"""Tests for muc_one_up.analysis.vntr_statistics module.

Tests cover:
- parse_vntr function with various dash characters
- analyze_vntr_sequences with header and without
- Error handling (empty files, invalid columns)
- Warning generation for unknown tokens
- Statistics calculations
- Transition probability matrix
"""

import io
import warnings

import pytest

from muc_one_up.analysis.vntr_statistics import analyze_vntr_sequences, parse_vntr


@pytest.mark.unit
class TestParseVNTR:
    """Tests for parse_vntr function."""

    def test_parse_standard_dash(self):
        """Test parsing with standard ASCII dash (U+002D)."""
        result = parse_vntr("1-2-3-X-A-B-6p-7-8-9")
        assert result == ["1", "2", "3", "X", "A", "B", "6p", "7", "8", "9"]

    def test_parse_unicode_dashes(self):
        """Test parsing with Unicode dash variants."""
        # U+2010 (hyphen)
        result = parse_vntr("1‐2‐3")  # noqa: RUF001
        assert result == ["1", "2", "3"]

        # U+2013 (en dash)
        result = parse_vntr("X–A–B")  # noqa: RUF001
        assert result == ["X", "A", "B"]

        # U+2014 (em dash)
        result = parse_vntr("6p—7—8")
        assert result == ["6p", "7", "8"]

    def test_parse_empty_string(self):
        """Test parsing empty string."""
        result = parse_vntr("")
        assert result == []

    def test_parse_whitespace_handling(self):
        """Test that leading/trailing whitespace is stripped."""
        result = parse_vntr("  1-2-3  ")
        assert result == ["1", "2", "3"]

    def test_parse_single_token(self):
        """Test parsing single repeat unit."""
        result = parse_vntr("X")
        assert result == ["X"]

    def test_parse_mixed_dashes(self):
        """Test parsing with mixed dash types."""
        # Mix of ASCII and Unicode dashes
        result = parse_vntr("1-2‐3–4")  # noqa: RUF001
        assert result == ["1", "2", "3", "4"]


@pytest.mark.unit
class TestAnalyzeVNTRSequences:
    """Tests for analyze_vntr_sequences function."""

    @pytest.fixture
    def known_repeats(self):
        """Fixture providing known repeat dictionary."""
        return {
            "1": "ACGT",
            "2": "GCTA",
            "3": "TGCA",
            "X": "ATCG",
            "A": "CGAT",
            "B": "TAGC",
            "6p": "CTAG",
            "7": "GATC",
            "8": "TCAG",
            "9": "AGCT",
        }

    def test_analyze_with_header(self, known_repeats):
        """Test analysis with header row."""
        data = """publication\tid\tallele\tvntr
PMID:123\tS1\t1\t1-2-3-X-A-B-6p-7-8-9
PMID:123\tS2\t2\tX-A-B-6p-7-8-9
"""
        file_handle = io.StringIO(data)
        result = analyze_vntr_sequences(
            file_handle, "vntr", has_header=True, known_repeats=known_repeats
        )

        assert result["min_repeats"] == 7
        assert result["max_repeats"] == 10
        assert result["mean_repeats"] == 8.5
        assert result["median_repeats"] == 8.5
        assert "probabilities" in result
        assert "END" in result["probabilities"]["9"]

    def test_analyze_without_header(self, known_repeats):
        """Test analysis without header row (column index)."""
        data = """PMID:123\tS1\t1\t1-2-3-X
PMID:123\tS2\t2\tX-A-B
"""
        file_handle = io.StringIO(data)
        result = analyze_vntr_sequences(
            file_handle, 3, has_header=False, known_repeats=known_repeats, delimiter="\t"
        )

        assert result["min_repeats"] == 3
        assert result["max_repeats"] == 4

    def test_duplicate_removal(self, known_repeats):
        """Test that duplicate VNTR structures are counted once."""
        data = """vntr
1-2-3
X-A-B
1-2-3
"""
        file_handle = io.StringIO(data)
        result = analyze_vntr_sequences(
            file_handle, "vntr", has_header=True, known_repeats=known_repeats
        )

        # Should only count 2 unique structures
        assert result["min_repeats"] == 3
        assert result["max_repeats"] == 3
        # Mean and median should be 3.0 (two sequences, both length 3)
        assert result["mean_repeats"] == 3.0
        assert result["median_repeats"] == 3.0

    def test_transition_probabilities(self, known_repeats):
        """Test transition probability calculation."""
        data = """vntr
1-2
1-3
1-2
"""
        file_handle = io.StringIO(data)
        result = analyze_vntr_sequences(
            file_handle, "vntr", has_header=True, known_repeats=known_repeats
        )

        # "1" transitions: 2 times to "2", 1 time to "3"
        # But wait - these are separate sequences, so "1" is the start of each
        # Actually "1" appears 3 times as first token, never in middle
        # Let me recalculate: only 2 unique sequences (1-2 and 1-3)
        # 1-2: "1" -> "2", "2" -> "END"
        # 1-3: "1" -> "3", "3" -> "END"
        # So "1" transitions: 1x to "2", 1x to "3" = 50% each
        assert result["probabilities"]["1"]["2"] == pytest.approx(0.5)
        assert result["probabilities"]["1"]["3"] == pytest.approx(0.5)

    def test_end_state_transitions(self, known_repeats):
        """Test that last tokens transition to END state."""
        data = """vntr
1-2-3
X-A-B
"""
        file_handle = io.StringIO(data)
        result = analyze_vntr_sequences(
            file_handle, "vntr", has_header=True, known_repeats=known_repeats
        )

        # "3" and "B" should both transition to END
        assert "END" in result["probabilities"]["3"]
        assert "END" in result["probabilities"]["B"]
        assert result["probabilities"]["3"]["END"] == 1.0
        assert result["probabilities"]["B"]["END"] == 1.0

    def test_unknown_tokens_warning(self, known_repeats):
        """Test that unknown tokens generate warnings."""
        data = """vntr
1-UNKNOWN-3
"""
        file_handle = io.StringIO(data)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            analyze_vntr_sequences(
                file_handle, "vntr", has_header=True, known_repeats=known_repeats
            )

            assert len(w) == 1
            assert "UNKNOWN" in str(w[0].message)

    def test_empty_file_error(self, known_repeats):
        """Test that empty file raises ValueError."""
        data = """vntr
"""
        file_handle = io.StringIO(data)

        with pytest.raises(ValueError, match="No valid VNTR sequences"):
            analyze_vntr_sequences(
                file_handle, "vntr", has_header=True, known_repeats=known_repeats
            )

    def test_invalid_column_index_error(self, known_repeats):
        """Test that invalid column index raises ValueError."""
        data = """1-2-3\tX-A-B"""
        file_handle = io.StringIO(data)

        with pytest.raises(ValueError, match="integer index"):
            analyze_vntr_sequences(
                file_handle, "invalid", has_header=False, known_repeats=known_repeats
            )

    def test_custom_delimiter(self, known_repeats):
        """Test analysis with custom delimiter."""
        data = """vntr
1-2-3"""
        file_handle = io.StringIO(data)

        result = analyze_vntr_sequences(
            file_handle,
            "vntr",
            has_header=True,
            known_repeats=known_repeats,
            delimiter=",",
        )

        assert "min_repeats" in result

    def test_missing_column_header(self, known_repeats):
        """Test behavior when specified column doesn't exist in header."""
        data = """publication\tid\tallele
PMID:123\tS1\t1
"""
        file_handle = io.StringIO(data)

        # Should result in empty file error since column doesn't exist
        with pytest.raises(ValueError, match="No valid VNTR sequences"):
            analyze_vntr_sequences(
                file_handle,
                "nonexistent_column",
                has_header=True,
                known_repeats=known_repeats,
            )


@pytest.mark.unit
class TestVNTRStatisticsEdgeCases:
    """Tests for edge cases and error conditions."""

    def test_single_sequence(self):
        """Test with single VNTR sequence."""
        data = """vntr
1-2-3"""
        known = {"1": "A", "2": "B", "3": "C"}
        file_handle = io.StringIO(data)

        result = analyze_vntr_sequences(file_handle, "vntr", has_header=True, known_repeats=known)

        assert result["min_repeats"] == 3
        assert result["max_repeats"] == 3
        assert result["mean_repeats"] == 3.0
        assert result["median_repeats"] == 3.0

    def test_very_long_sequence(self):
        """Test with very long VNTR sequence."""
        repeats = [f"R{i}" for i in range(100)]
        data = f"""vntr
{"-".join(repeats)}"""
        known = dict.fromkeys(repeats, "ACGT")
        file_handle = io.StringIO(data)

        result = analyze_vntr_sequences(file_handle, "vntr", has_header=True, known_repeats=known)

        assert result["min_repeats"] == 100
        assert result["max_repeats"] == 100
        assert result["mean_repeats"] == 100.0
        assert result["median_repeats"] == 100.0

    def test_all_same_transitions(self):
        """Test when all sequences have identical structure."""
        data = """vntr
X-A-B
X-A-B
X-A-B
"""
        known = {"X": "A", "A": "B", "B": "C"}
        file_handle = io.StringIO(data)

        result = analyze_vntr_sequences(file_handle, "vntr", has_header=True, known_repeats=known)

        # Only 1 unique sequence after deduplication
        assert result["min_repeats"] == 3
        assert result["max_repeats"] == 3
        # Transition probabilities should be deterministic
        assert result["probabilities"]["X"]["A"] == 1.0
        assert result["probabilities"]["A"]["B"] == 1.0
        assert result["probabilities"]["B"]["END"] == 1.0

    def test_column_index_out_of_range(self):
        """Test with column index beyond row length."""
        data = """col1\tcol2
val1\tval2
"""
        known = {"1": "A"}
        file_handle = io.StringIO(data)

        # Column index 5 doesn't exist, should skip rows and raise error
        with pytest.raises(ValueError, match="No valid VNTR sequences"):
            analyze_vntr_sequences(file_handle, 5, has_header=False, known_repeats=known)
