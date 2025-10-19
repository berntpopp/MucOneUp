#!/usr/bin/env python3
"""
Reference file utilities for diploid simulation.

This module provides utilities for analyzing and manipulating FASTA reference files,
with special support for diploid references (2 haplotypes).
"""

import logging
from pathlib import Path
from typing import NamedTuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ...exceptions import ValidationError


class ReferenceInfo(NamedTuple):
    """Information about a reference FASTA file.

    Attributes:
        num_sequences: Number of sequences in the FASTA file.
        sequence_ids: List of sequence IDs.
        sequence_lengths: List of sequence lengths.
        total_length: Total length of all sequences combined.
        is_diploid: True if file contains exactly 2 sequences.
    """
    num_sequences: int
    sequence_ids: list[str]
    sequence_lengths: list[int]
    total_length: int
    is_diploid: bool


def get_reference_info(fasta_path: str | Path) -> ReferenceInfo:
    """
    Get information about a FASTA reference file.

    Args:
        fasta_path: Path to FASTA file.

    Returns:
        ReferenceInfo object containing reference metadata.

    Raises:
        ValidationError: If file doesn't exist or is not a valid FASTA.

    Example:
        >>> info = get_reference_info("diploid.fa")
        >>> print(f"Diploid: {info.is_diploid}, Sequences: {info.num_sequences}")
        Diploid: True, Sequences: 2
    """
    fasta_path = Path(fasta_path)

    if not fasta_path.exists():
        raise ValidationError(f"Reference file not found: {fasta_path}")

    if not fasta_path.is_file():
        raise ValidationError(f"Reference path is not a file: {fasta_path}")

    try:
        sequences = list(SeqIO.parse(str(fasta_path), "fasta"))
    except Exception as e:
        raise ValidationError(
            f"Failed to parse FASTA file {fasta_path}: {e}"
        ) from e

    if not sequences:
        raise ValidationError(f"FASTA file is empty: {fasta_path}")

    sequence_ids = [seq.id for seq in sequences]
    sequence_lengths = [len(seq.seq) for seq in sequences]
    total_length = sum(sequence_lengths)
    num_sequences = len(sequences)
    is_diploid = num_sequences == 2

    logging.debug(
        "Reference info for %s: %d sequences, %d total bp, diploid=%s",
        fasta_path.name,
        num_sequences,
        total_length,
        is_diploid,
    )

    return ReferenceInfo(
        num_sequences=num_sequences,
        sequence_ids=sequence_ids,
        sequence_lengths=sequence_lengths,
        total_length=total_length,
        is_diploid=is_diploid,
    )


def is_diploid_reference(fasta_path: str | Path) -> bool:
    """
    Determine if a FASTA file contains a diploid reference (exactly 2 sequences).

    This is a convenience function that returns only the diploid status.
    For complete reference information, use get_reference_info().

    Args:
        fasta_path: Path to FASTA file.

    Returns:
        True if file contains exactly 2 sequences, False otherwise.

    Raises:
        ValidationError: If file doesn't exist or is not a valid FASTA.

    Example:
        >>> if is_diploid_reference("sample.fa"):
        ...     print("Diploid reference detected")
    """
    info = get_reference_info(fasta_path)
    return info.is_diploid


def extract_haplotypes(
    diploid_fasta: str | Path,
    output_dir: str | Path,
    base_name: str | None = None,
) -> tuple[str, str]:
    """
    Extract individual haplotypes from a diploid FASTA file.

    Splits a diploid reference (2 sequences) into two separate FASTA files,
    one per haplotype. This is required for split-simulation to eliminate
    length-proportional sampling bias.

    Args:
        diploid_fasta: Path to diploid FASTA file (must contain exactly 2 sequences).
        output_dir: Directory where haplotype files will be created.
        base_name: Base name for output files (default: use input filename).
                  Output files will be named {base_name}_hap1.fa and {base_name}_hap2.fa.

    Returns:
        Tuple of (hap1_path, hap2_path) as strings.

    Raises:
        ValidationError: If input is not a valid diploid reference.

    Example:
        >>> hap1, hap2 = extract_haplotypes("diploid.fa", "output", "sample")
        >>> print(hap1)
        output/sample_hap1.fa

    Notes:
        - Haplotypes are numbered in the order they appear in the input file
        - Original sequence IDs are preserved in output files
        - Output directory is created if it doesn't exist
    """
    diploid_fasta = Path(diploid_fasta)
    output_dir = Path(output_dir)

    # Validate diploid reference
    info = get_reference_info(diploid_fasta)
    if not info.is_diploid:
        raise ValidationError(
            f"Expected diploid reference (2 sequences), found {info.num_sequences} "
            f"sequences in {diploid_fasta}"
        )

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine base name
    if base_name is None:
        base_name = diploid_fasta.stem

    # Parse sequences
    sequences = list(SeqIO.parse(str(diploid_fasta), "fasta"))
    hap1_seq = sequences[0]
    hap2_seq = sequences[1]

    # Create output paths
    hap1_path = output_dir / f"{base_name}_hap1.fa"
    hap2_path = output_dir / f"{base_name}_hap2.fa"

    # Write haplotype files
    try:
        SeqIO.write([hap1_seq], str(hap1_path), "fasta")
        SeqIO.write([hap2_seq], str(hap2_path), "fasta")
    except Exception as e:
        raise ValidationError(
            f"Failed to write haplotype files: {e}"
        ) from e

    logging.info(
        "Extracted haplotypes from %s:",
        diploid_fasta.name,
    )
    logging.info(
        "  Haplotype 1 (%s): %d bp -> %s",
        hap1_seq.id,
        len(hap1_seq.seq),
        hap1_path.name,
    )
    logging.info(
        "  Haplotype 2 (%s): %d bp -> %s",
        hap2_seq.id,
        len(hap2_seq.seq),
        hap2_path.name,
    )

    return str(hap1_path), str(hap2_path)


def validate_reference_compatibility(
    reference_path: str | Path,
    min_sequences: int = 1,
    max_sequences: int | None = None,
) -> None:
    """
    Validate that a reference file meets sequence count requirements.

    Args:
        reference_path: Path to FASTA file.
        min_sequences: Minimum number of sequences required.
        max_sequences: Maximum number of sequences allowed (None = no limit).

    Raises:
        ValidationError: If reference doesn't meet requirements.

    Example:
        >>> # Ensure reference is diploid
        >>> validate_reference_compatibility("sample.fa", min_sequences=2, max_sequences=2)
    """
    info = get_reference_info(reference_path)

    if info.num_sequences < min_sequences:
        raise ValidationError(
            f"Reference {reference_path} has {info.num_sequences} sequences, "
            f"but at least {min_sequences} required"
        )

    if max_sequences is not None and info.num_sequences > max_sequences:
        raise ValidationError(
            f"Reference {reference_path} has {info.num_sequences} sequences, "
            f"but maximum {max_sequences} allowed"
        )

    logging.debug(
        "Reference validation passed: %s (%d sequences)",
        Path(reference_path).name,
        info.num_sequences,
    )
