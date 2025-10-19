# muc_one_up/fasta_writer.py
"""FASTA file writing utilities for MucOneUp.

This module provides FASTA output functionality with support for custom
headers, per-sequence comments, and proper line wrapping (80 characters).

Key Functions:
    write_fasta: Write DNA/protein sequences to FASTA file with formatting

Features:
    - Automatic 80-character line wrapping for readability
    - Per-sequence comment support (for mutation annotations)
    - Global comment fallback for all sequences
    - Safe file handling with proper error reporting

Example:
    Write haplotype sequences with mutation annotations::

        from muc_one_up.fasta_writer import write_fasta

        sequences = [seq1, seq2]
        comments = [
            "mutation=dupC targets=[(1, 25)]",
            "mutation=normal"
        ]
        write_fasta(
            sequences=sequences,
            filename="output.fa",
            prefix="haplotype",
            comments=comments
        )

Output Format:
    >haplotype_1 mutation=dupC targets=[(1, 25)]
    ATCGATCGATCGATCG...
    >haplotype_2 mutation=normal
    GCTAGCTAGCTAGCTA...
"""

import logging
from pathlib import Path


def write_fasta(
    sequences: list[str],
    filename: str,
    prefix: str = "haplotype",
    comment: str | None = None,
    comments: list[str | None] | None = None,
) -> None:
    """Write DNA/protein sequences to FASTA file with formatting.

    Writes sequences with auto-generated sequential IDs (prefix_1, prefix_2, ...)
    and optional comments. Sequences are wrapped at 80 characters for readability.

    Args:
        sequences: List of DNA or protein sequence strings
        filename: Output FASTA filename
        prefix: Header prefix for sequential IDs, defaults to "haplotype"
        comment: Optional comment appended to all headers (if comments not provided)
        comments: Optional list of per-sequence comments (overrides comment parameter)

    Raises:
        Exception: If file writing fails (logged and re-raised)

    Example:
        Write haplotypes with mutation annotations::

            sequences = [seq1, seq2]
            comments = ["mutation=dupC", "mutation=normal"]
            write_fasta(sequences, "output.fa", comments=comments)

        Output::

            >haplotype_1 mutation=dupC
            ATCGATCGATCG...
            >haplotype_2 mutation=normal
            GCTAGCTAGCTA...

    Note:
        Per-sequence comments take precedence over global comment parameter.
    """
    try:
        with Path(filename).open("w") as fh:
            for i, seq in enumerate(sequences, start=1):
                header = f">{prefix}_{i}"

                # Use sequence-specific comment if available
                if comments and i - 1 < len(comments) and comments[i - 1]:
                    header += f" {comments[i - 1]}"
                elif comment:  # Fall back to global comment
                    header += f" {comment}"

                # Write the sequence in chunks of 80 characters for readability
                fh.write(f"{header}\n")
                for j in range(0, len(seq), 80):
                    fh.write(f"{seq[j : j + 80]}\n")
        logging.info("FASTA file written to %s", filename)
    except Exception as e:
        logging.error("Error writing FASTA file: %s", e)
        raise
