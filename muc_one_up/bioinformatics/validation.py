"""Bioinformatics-specific validation for DNA sequences and structures.

This module provides validation for:
- DNA sequence validity (bases, format)
- FASTA file format
- Repeat structure strings
- SNP records
"""

import re
from pathlib import Path

from ..exceptions import ValidationError
from ..type_defs import DNASequence, RepeatStructure

# Valid DNA bases
DNA_BASES: set[str] = {"A", "C", "G", "T", "N"}
DNA_PATTERN = re.compile(r"^[ACGTN]+$", re.IGNORECASE)

# Valid SNP bases (no ambiguous bases allowed for SNPs)
SNP_BASES: set[str] = {"A", "C", "G", "T"}


def validate_dna_sequence(sequence: DNASequence, allow_ambiguous: bool = True) -> None:
    """Validate DNA sequence contains only valid bases.

    Args:
        sequence: DNA sequence string
        allow_ambiguous: Allow N for ambiguous bases (default: True)

    Raises:
        ValidationError: If sequence contains invalid characters or is empty
    """
    if not sequence:
        raise ValidationError("DNA sequence cannot be empty")

    # Convert to uppercase for validation
    bases = set(sequence.upper())
    allowed = DNA_BASES if allow_ambiguous else DNA_BASES - {"N"}
    invalid = bases - allowed

    if invalid:
        raise ValidationError(f"Invalid DNA bases: {sorted(invalid)}. Allowed: {sorted(allowed)}")


def validate_fasta_format(fasta_path: str | Path) -> None:
    """Validate FASTA file format.

    Checks that:
    - File is not empty
    - First line starts with '>'
    - Sequence lines contain only valid DNA bases

    Args:
        fasta_path: Path to FASTA file

    Raises:
        ValidationError: If FASTA format is invalid
        FileNotFoundError: If file doesn't exist
    """
    path = Path(fasta_path)

    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    with path.open() as f:
        lines = [line.strip() for line in f if line.strip()]

    if not lines:
        raise ValidationError(f"FASTA file is empty: {fasta_path}")

    if not lines[0].startswith(">"):
        raise ValidationError(f"FASTA must start with header line (>): {fasta_path}")

    # Validate sequence lines
    in_sequence = False
    for i, line in enumerate(lines, 1):
        if line.startswith(">"):
            in_sequence = True
            continue

        if in_sequence and not DNA_PATTERN.match(line):
            # Truncate long lines for error message (inline to avoid mypy unreachable warning)
            raise ValidationError(
                f"Invalid DNA sequence at line {i}: {line[:50] + '...' if len(line) > 50 else line}"
            )


def validate_repeat_structure(
    structure: RepeatStructure,
    valid_symbols: set[str],
) -> None:
    """Validate repeat structure string.

    Args:
        structure: Hyphen-separated repeat structure (e.g., "1-2-7-8-9")
        valid_symbols: Set of valid repeat symbols from config

    Raises:
        ValidationError: If structure is invalid
    """
    if not structure:
        raise ValidationError("Repeat structure cannot be empty")

    # Use ternary operator for simplicity
    repeats = [structure] if "-" not in structure else structure.split("-")

    for repeat in repeats:
        # Remove mutation marker if present
        clean_repeat = repeat.rstrip("m")

        if clean_repeat not in valid_symbols:
            raise ValidationError(
                f"Invalid repeat symbol '{repeat}' in structure. "
                f"Valid symbols: {sorted(valid_symbols)}"
            )


def validate_snp_base(base: str, allow_ref_mismatch: bool = False) -> None:
    """Validate SNP base is a valid nucleotide.

    Args:
        base: Single nucleotide base
        allow_ref_mismatch: Allow bases that don't match typical SNPs

    Raises:
        ValidationError: If base is invalid
    """
    if not base:
        raise ValidationError("SNP base cannot be empty")

    upper_base = base.upper()
    if upper_base not in SNP_BASES:
        raise ValidationError(f"Invalid SNP base: '{base}'. Allowed: {sorted(SNP_BASES)}")


def validate_snp_record(
    haplotype: int,
    position: int,
    ref: str,
    alt: str,
    num_haplotypes: int,
    sequence_length: int,
) -> None:
    """Validate a complete SNP record.

    Args:
        haplotype: 1-based haplotype index
        position: 0-based position in sequence
        ref: Reference base
        alt: Alternate base
        num_haplotypes: Total number of haplotypes
        sequence_length: Length of sequence

    Raises:
        ValidationError: If SNP record is invalid
    """
    # Validate haplotype index (1-based)
    if not 1 <= haplotype <= num_haplotypes:
        raise ValidationError(
            f"SNP haplotype index {haplotype} out of range. "
            f"Must be 1-{num_haplotypes} (1-based indexing)."
        )

    # Validate position (0-based)
    if not 0 <= position < sequence_length:
        raise ValidationError(
            f"SNP position {position} out of range. "
            f"Must be 0-{sequence_length - 1} (0-based indexing)."
        )

    # Validate bases
    validate_snp_base(ref)
    validate_snp_base(alt)

    # Validate ref != alt
    if ref.upper() == alt.upper():
        raise ValidationError(f"SNP reference and alternate bases must differ: {ref} == {alt}")


def validate_gc_content_range(gc_content: float) -> None:
    """Validate GC content is in valid range.

    Args:
        gc_content: GC content percentage (0.0-100.0)

    Raises:
        ValidationError: If GC content is invalid
    """
    if not 0.0 <= gc_content <= 100.0:
        raise ValidationError(f"GC content must be between 0.0 and 100.0, got: {gc_content}")


def validate_sequence_length(
    length: int, min_length: int = 0, max_length: int | None = None
) -> None:
    """Validate sequence length is within bounds.

    Args:
        length: Sequence length
        min_length: Minimum allowed length (default: 0)
        max_length: Maximum allowed length (optional)

    Raises:
        ValidationError: If length is invalid
    """
    if length < min_length:
        raise ValidationError(f"Sequence length {length} below minimum: {min_length}")

    if max_length is not None and length > max_length:
        raise ValidationError(f"Sequence length {length} exceeds maximum: {max_length}")
