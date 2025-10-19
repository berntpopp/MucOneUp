"""Test utilities for MucOneUp tests.

This module provides reusable test utility functions following DRY principles.
These utilities help with common test assertions and validations.
"""

from collections import Counter
from pathlib import Path


def assert_valid_fasta(fasta_path: Path) -> None:
    """Assert that FASTA file is valid.

    Checks:
    - File is not empty
    - Starts with header line (>)
    - Contains only valid DNA bases (ACGTN)

    Args:
        fasta_path: Path to FASTA file

    Raises:
        AssertionError: If FASTA format is invalid
    """
    with fasta_path.open() as f:
        lines = [line.strip() for line in f if line.strip()]

    assert lines, "FASTA file is empty"
    assert lines[0].startswith(">"), "FASTA must start with header"

    # Check sequence lines contain only valid DNA bases
    for i, line in enumerate(lines):
        if line.startswith(">"):
            continue
        invalid_bases = set(line.upper()) - {"A", "C", "G", "T", "N"}
        assert not invalid_bases, f"Invalid DNA bases at line {i}: {invalid_bases}"


def assert_structure_file_valid(structure_path: Path, num_haplotypes: int = 2) -> None:
    """Assert that structure file is valid.

    Checks:
    - Correct number of haplotypes
    - Valid line format (name<TAB>structure)
    - Valid haplotype names
    - Valid structure format (dash-separated)

    Args:
        structure_path: Path to structure file
        num_haplotypes: Expected number of haplotypes

    Raises:
        AssertionError: If structure format is invalid
    """
    with structure_path.open() as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith("#")]

    assert len(lines) == num_haplotypes, f"Expected {num_haplotypes} haplotypes, got {len(lines)}"

    for line in lines:
        parts = line.split("\t")
        assert len(parts) == 2, f"Invalid structure line (expected TAB separator): {line}"
        haplotype_name, structure = parts
        assert haplotype_name.startswith("haplotype"), f"Invalid haplotype name: {haplotype_name}"
        assert "-" in structure, f"Invalid structure format (expected dashes): {structure}"


def count_sequence_bases(sequence: str) -> dict[str, int]:
    """Count bases in DNA sequence.

    Args:
        sequence: DNA sequence string

    Returns:
        Dictionary mapping base to count
    """
    return dict(Counter(sequence.upper()))


def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content of DNA sequence.

    Args:
        sequence: DNA sequence string

    Returns:
        GC content as percentage (0-100)
    """
    counts = count_sequence_bases(sequence)
    gc = counts.get("G", 0) + counts.get("C", 0)
    total = len(sequence)
    return (gc / total) * 100 if total > 0 else 0.0


def count_fasta_records(fasta_path: Path) -> int:
    """Count number of records in FASTA file.

    Args:
        fasta_path: Path to FASTA file

    Returns:
        Number of sequences in FASTA file
    """
    with fasta_path.open() as f:
        return sum(1 for line in f if line.startswith(">"))


def get_fasta_sequences(fasta_path: Path) -> dict[str, str]:
    """Parse FASTA file and return sequences.

    Args:
        fasta_path: Path to FASTA file

    Returns:
        Dictionary mapping sequence name to sequence
    """
    sequences = {}
    current_name = None
    current_seq = []

    with fasta_path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # Save previous sequence
                if current_name:
                    sequences[current_name] = "".join(current_seq)
                # Start new sequence
                current_name = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_name:
            sequences[current_name] = "".join(current_seq)

    return sequences


def assert_fasta_has_sequences(fasta_path: Path, expected_names: list[str]) -> None:
    """Assert that FASTA file contains expected sequences.

    Args:
        fasta_path: Path to FASTA file
        expected_names: List of expected sequence names

    Raises:
        AssertionError: If sequences don't match
    """
    sequences = get_fasta_sequences(fasta_path)
    actual_names = set(sequences.keys())
    expected_set = set(expected_names)

    assert actual_names == expected_set, (
        f"Sequence name mismatch.\n"
        f"Expected: {expected_set}\n"
        f"Got: {actual_names}\n"
        f"Missing: {expected_set - actual_names}\n"
        f"Extra: {actual_names - expected_set}"
    )


def assert_file_not_empty(file_path: Path) -> None:
    """Assert that file exists and is not empty.

    Args:
        file_path: Path to file

    Raises:
        AssertionError: If file doesn't exist or is empty
    """
    assert file_path.exists(), f"File does not exist: {file_path}"
    assert file_path.stat().st_size > 0, f"File is empty: {file_path}"
