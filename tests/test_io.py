"""Tests for muc_one_up.io module.

Tests cover:
- VNTR structure file parsing
- Mutation information extraction from comments
- Error handling for invalid formats
- Edge cases
"""

from pathlib import Path

import pytest

from muc_one_up.io import (
    extract_mutation_info_from_comments,
    parse_vntr_structure_file,
)


@pytest.mark.unit
class TestParseVntrStructureFile:
    """Tests for parse_vntr_structure_file function."""

    def test_parse_valid_structure_file(self, tmp_path: Path, minimal_config: dict):
        """Test parsing valid structure file succeeds."""
        structure_file = tmp_path / "structure.txt"
        structure_file.write_text("haplotype_1\t1-2-X-B-6-7-8-9\nhaplotype_2\t1-2-A-B-6p-7-8-9\n")

        chains, mutation_info = parse_vntr_structure_file(str(structure_file), minimal_config)

        assert len(chains) == 2
        assert chains[0] == ["1", "2", "X", "B", "6", "7", "8", "9"]
        assert chains[1] == ["1", "2", "A", "B", "6p", "7", "8", "9"]
        assert mutation_info is None

    def test_parse_structure_file_with_mutation_markers(self, tmp_path: Path, minimal_config: dict):
        """Test parsing structure file with mutation markers (m suffix)."""
        structure_file = tmp_path / "structure.txt"
        structure_file.write_text(
            "# Mutation Applied: dupC (Targets: [(1, 25)])\n"
            "haplotype_1\t1-2-Xm-B-6-7-8-9\n"
            "haplotype_2\t1-2-A-Bm-6p-7-8-9\n"
        )

        chains, mutation_info = parse_vntr_structure_file(str(structure_file), minimal_config)

        assert len(chains) == 2
        # Mutation markers should be preserved in the chain
        assert chains[0] == ["1", "2", "Xm", "B", "6", "7", "8", "9"]
        assert chains[1] == ["1", "2", "A", "Bm", "6p", "7", "8", "9"]
        # Mutation info should be extracted
        assert mutation_info is not None
        assert mutation_info["name"] == "dupC"
        assert mutation_info["targets"] == [(1, 25)]

    def test_parse_structure_file_with_comments(self, tmp_path: Path, minimal_config: dict):
        """Test parsing structure file with comment lines."""
        structure_file = tmp_path / "structure.txt"
        structure_file.write_text(
            "# This is a comment\n"
            "# Another comment line\n"
            "haplotype_1\t1-2-X-B-6-7-8-9\n"
            "# Comment in the middle\n"
            "haplotype_2\t1-2-A-B-6p-7-8-9\n"
        )

        chains, _mutation_info = parse_vntr_structure_file(str(structure_file), minimal_config)

        assert len(chains) == 2
        assert chains[0] == ["1", "2", "X", "B", "6", "7", "8", "9"]
        assert chains[1] == ["1", "2", "A", "B", "6p", "7", "8", "9"]

    def test_parse_structure_file_with_empty_lines(self, tmp_path: Path, minimal_config: dict):
        """Test parsing structure file with empty lines."""
        structure_file = tmp_path / "structure.txt"
        structure_file.write_text(
            "\nhaplotype_1\t1-2-X-B-6-7-8-9\n\n\nhaplotype_2\t1-2-A-B-6p-7-8-9\n\n"
        )

        chains, _mutation_info = parse_vntr_structure_file(str(structure_file), minimal_config)

        assert len(chains) == 2

    def test_parse_structure_file_single_haplotype(self, tmp_path: Path, minimal_config: dict):
        """Test parsing structure file with single haplotype."""
        structure_file = tmp_path / "structure.txt"
        structure_file.write_text("haplotype_1\t1-2-X-B-6-7-8-9\n")

        chains, _mutation_info = parse_vntr_structure_file(str(structure_file), minimal_config)

        assert len(chains) == 1
        assert chains[0] == ["1", "2", "X", "B", "6", "7", "8", "9"]

    def test_parse_structure_file_many_haplotypes(self, tmp_path: Path, minimal_config: dict):
        """Test parsing structure file with many haplotypes."""
        structure_file = tmp_path / "structure.txt"
        lines = [f"haplotype_{i}\t1-2-X-B-6-7-8-9\n" for i in range(1, 11)]
        structure_file.write_text("".join(lines))

        chains, _mutation_info = parse_vntr_structure_file(str(structure_file), minimal_config)

        assert len(chains) == 10
        for chain in chains:
            assert chain == ["1", "2", "X", "B", "6", "7", "8", "9"]

    def test_parse_structure_file_missing_file(self, minimal_config: dict):
        """Test parsing nonexistent file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError) as exc_info:
            parse_vntr_structure_file("/nonexistent/file.txt", minimal_config)

        assert "Structure file not found" in str(exc_info.value)

    def test_parse_structure_file_invalid_repeat_symbol(self, tmp_path: Path, minimal_config: dict):
        """Test parsing structure with invalid repeat symbol raises ValueError."""
        structure_file = tmp_path / "structure.txt"
        structure_file.write_text("haplotype_1\t1-2-INVALID-6-7-8-9\n")

        with pytest.raises(ValueError) as exc_info:
            parse_vntr_structure_file(str(structure_file), minimal_config)

        assert "Invalid repeat symbol 'INVALID'" in str(exc_info.value)

    def test_parse_structure_file_missing_tab_separator(self, tmp_path: Path, minimal_config: dict):
        """Test parsing structure without tab separator raises ValueError."""
        structure_file = tmp_path / "structure.txt"
        structure_file.write_text("haplotype_1 1-2-X-B-6-7-8-9\n")  # Space instead of tab

        with pytest.raises(ValueError) as exc_info:
            parse_vntr_structure_file(str(structure_file), minimal_config)

        assert "Expected tab-separated columns" in str(exc_info.value)

    def test_parse_structure_file_empty_file(self, tmp_path: Path, minimal_config: dict):
        """Test parsing empty file raises ValueError."""
        structure_file = tmp_path / "structure.txt"
        structure_file.write_text("")

        with pytest.raises(ValueError) as exc_info:
            parse_vntr_structure_file(str(structure_file), minimal_config)

        assert "No valid data found" in str(exc_info.value)

    def test_parse_structure_file_only_comments(self, tmp_path: Path, minimal_config: dict):
        """Test parsing file with only comments raises ValueError."""
        structure_file = tmp_path / "structure.txt"
        structure_file.write_text("# Comment 1\n# Comment 2\n# Comment 3\n")

        with pytest.raises(ValueError) as exc_info:
            parse_vntr_structure_file(str(structure_file), minimal_config)

        assert "No valid data found" in str(exc_info.value)

    def test_parse_structure_file_complex_chains(self, tmp_path: Path, minimal_config: dict):
        """Test parsing structure with long complex chains."""
        structure_file = tmp_path / "structure.txt"
        structure_file.write_text(
            "haplotype_1\t1-2-X-X-X-A-A-A-B-B-B-C-6-6p-6-6p-7-8-9\n"
            "haplotype_2\t1-2-A-B-X-C-A-B-X-C-6-7-8-9\n"
        )

        chains, _mutation_info = parse_vntr_structure_file(str(structure_file), minimal_config)

        assert len(chains) == 2
        assert len(chains[0]) == 19  # 1-2-X-X-X-A-A-A-B-B-B-C-6-6p-6-6p-7-8-9
        assert len(chains[1]) == 14  # 1-2-A-B-X-C-A-B-X-C-6-7-8-9


@pytest.mark.unit
class TestExtractMutationInfo:
    """Tests for extract_mutation_info_from_comments function."""

    def test_extract_mutation_info_valid(self):
        """Test extracting valid mutation information."""
        comments = ["Mutation Applied: dupC (Targets: [(1, 25)])"]

        mutation_info = extract_mutation_info_from_comments(comments)

        assert mutation_info is not None
        assert mutation_info["name"] == "dupC"
        assert mutation_info["targets"] == [(1, 25)]

    def test_extract_mutation_info_multiple_targets(self):
        """Test extracting mutation info with multiple targets."""
        comments = ["Mutation Applied: delA (Targets: [(1, 10), (2, 20), (1, 30)])"]

        mutation_info = extract_mutation_info_from_comments(comments)

        assert mutation_info is not None
        assert mutation_info["name"] == "delA"
        assert mutation_info["targets"] == [(1, 10), (2, 20), (1, 30)]

    def test_extract_mutation_info_with_other_comments(self):
        """Test extracting mutation info from mixed comments."""
        comments = [
            "This is a regular comment",
            "Another comment line",
            "Mutation Applied: dupC (Targets: [(1, 25)])",
            "More comments after",
        ]

        mutation_info = extract_mutation_info_from_comments(comments)

        assert mutation_info is not None
        assert mutation_info["name"] == "dupC"
        assert mutation_info["targets"] == [(1, 25)]

    def test_extract_mutation_info_extra_whitespace(self):
        """Test extracting mutation info with extra whitespace."""
        comments = ["  Mutation Applied:   dupC   (  Targets:  [ ( 1 , 25 ) ]  )  "]

        mutation_info = extract_mutation_info_from_comments(comments)

        assert mutation_info is not None
        assert mutation_info["name"] == "dupC"

    def test_extract_mutation_info_no_mutation_comment(self):
        """Test extracting from comments with no mutation info."""
        comments = ["This is a comment", "Another comment", "No mutation info here"]

        mutation_info = extract_mutation_info_from_comments(comments)

        assert mutation_info is None

    def test_extract_mutation_info_empty_list(self):
        """Test extracting from empty comments list."""
        mutation_info = extract_mutation_info_from_comments([])

        assert mutation_info is None

    def test_extract_mutation_info_malformed_targets(self):
        """Test extracting mutation info with malformed targets."""
        # Malformed: not a list of tuples
        comments = ["Mutation Applied: dupC (Targets: [1, 25])"]

        mutation_info = extract_mutation_info_from_comments(comments)

        # Should return None for malformed targets
        assert mutation_info is None

    def test_extract_mutation_info_invalid_python_syntax(self):
        """Test extracting mutation info with invalid Python syntax."""
        comments = ["Mutation Applied: dupC (Targets: [(1, 25]"]  # Missing closing paren

        mutation_info = extract_mutation_info_from_comments(comments)

        # Should return None for syntax errors
        assert mutation_info is None

    def test_extract_mutation_info_first_valid_used(self):
        """Test that first valid mutation info is used when multiple present."""
        comments = [
            "Mutation Applied: dupC (Targets: [(1, 25)])",
            "Mutation Applied: delA (Targets: [(2, 30)])",  # This should be ignored
        ]

        mutation_info = extract_mutation_info_from_comments(comments)

        assert mutation_info is not None
        assert mutation_info["name"] == "dupC"
        assert mutation_info["targets"] == [(1, 25)]


@pytest.mark.bioinformatics
class TestStructureFileBioinformatics:
    """Bioinformatics-specific tests for structure file parsing."""

    def test_canonical_terminal_block_present(self, tmp_path: Path, minimal_config: dict):
        """Test that parsed chains can include canonical terminal blocks."""
        structure_file = tmp_path / "structure.txt"
        structure_file.write_text(
            "haplotype_1\t1-2-X-B-6-7-8-9\n"  # 6 → 7 → 8 → 9 is canonical
            "haplotype_2\t1-2-A-B-6p-7-8-9\n"  # 6p → 7 → 8 → 9 is canonical
        )

        chains, _ = parse_vntr_structure_file(str(structure_file), minimal_config)

        # Check that canonical terminal sequences are preserved
        for chain in chains:
            assert chain[-4:] == ["6", "7", "8", "9"] or chain[-4:] == ["6p", "7", "8", "9"]

    def test_repeat_symbols_match_config(self, tmp_path: Path, minimal_config: dict):
        """Test that all symbols in parsed chains are valid per config."""
        structure_file = tmp_path / "structure.txt"
        structure_file.write_text("haplotype_1\t1-2-X-A-B-C-6-7-8-9\n")

        chains, _ = parse_vntr_structure_file(str(structure_file), minimal_config)

        valid_symbols = set(minimal_config["repeats"].keys())

        for chain in chains:
            for symbol in chain:
                # Strip mutation marker if present
                pure_symbol = symbol.rstrip("m")
                assert pure_symbol in valid_symbols
