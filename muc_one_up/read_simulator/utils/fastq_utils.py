#!/usr/bin/env python3
"""
FASTQ file utilities for read simulation.

This module provides utilities for manipulating and validating FASTQ files,
including merging files from split-simulation workflows.
"""

import gzip
import logging
from pathlib import Path

from ...exceptions import ValidationError


def validate_fastq(fastq_path: str | Path, check_quality: bool = True) -> bool:
    """
    Validate that a file is a properly formatted FASTQ file.

    Performs basic validation:
    - File exists and is readable
    - Has .fastq or .fq extension (optionally .gz)
    - Contains valid FASTQ records (4 lines per record)

    Args:
        fastq_path: Path to FASTQ file (can be gzipped).
        check_quality: If True, verify quality scores match sequence length.

    Returns:
        True if file is valid.

    Raises:
        ValidationError: If file is invalid or malformed.

    Example:
        >>> validate_fastq("reads.fastq.gz")
        True
    """
    fastq_path = Path(fastq_path)

    if not fastq_path.exists():
        raise ValidationError(f"FASTQ file not found: {fastq_path}")

    if not fastq_path.is_file():
        raise ValidationError(f"FASTQ path is not a file: {fastq_path}")

    # Check extension
    valid_extensions = {".fastq", ".fq", ".fastq.gz", ".fq.gz"}
    if not any(str(fastq_path).endswith(ext) for ext in valid_extensions):
        logging.warning(
            "File %s does not have standard FASTQ extension (.fastq, .fq, or .gz variants)",
            fastq_path.name,
        )

    # Open file (handle gzip if needed)
    try:
        if str(fastq_path).endswith(".gz"):
            file_handle = gzip.open(fastq_path, "rt")
        else:
            file_handle = open(fastq_path)
    except Exception as e:
        raise ValidationError(f"Cannot open FASTQ file {fastq_path}: {e}") from e

    try:
        # Read first record to validate format
        line_num = 0
        for i in range(4):  # FASTQ records are 4 lines
            line = file_handle.readline()
            line_num += 1
            if not line:
                raise ValidationError(f"FASTQ file {fastq_path} is too short (fewer than 4 lines)")

            if i == 0 and not line.startswith("@"):
                raise ValidationError(
                    f"FASTQ file {fastq_path} line {line_num}: "
                    f"Expected header starting with '@', got: {line[:20]}"
                )

            if i == 2 and not line.startswith("+"):
                raise ValidationError(
                    f"FASTQ file {fastq_path} line {line_num}: "
                    f"Expected separator '+', got: {line[:20]}"
                )

        if check_quality:
            # The first 4 lines we just read are: header, seq, +, qual
            # Reopen to check properly
            file_handle.close()
            if str(fastq_path).endswith(".gz"):
                file_handle = gzip.open(fastq_path, "rt")
            else:
                file_handle = open(fastq_path)

            file_handle.readline().strip()
            seq = file_handle.readline().strip()
            file_handle.readline().strip()
            qual = file_handle.readline().strip()

            if len(seq) != len(qual):
                raise ValidationError(
                    f"FASTQ file {fastq_path}: "
                    f"Sequence length ({len(seq)}) != quality length ({len(qual)}) "
                    f"in first record"
                )

    except ValidationError:
        raise
    except Exception as e:
        raise ValidationError(f"Error validating FASTQ file {fastq_path}: {e}") from e
    finally:
        file_handle.close()

    logging.debug("FASTQ validation passed: %s", fastq_path.name)
    return True


def count_fastq_reads(fastq_path: str | Path) -> int:
    """
    Count the number of reads in a FASTQ file.

    Args:
        fastq_path: Path to FASTQ file (can be gzipped).

    Returns:
        Number of reads in the file.

    Raises:
        ValidationError: If file cannot be read or is malformed.

    Example:
        >>> count = count_fastq_reads("reads.fastq.gz")
        >>> print(f"Found {count} reads")
        Found 1000 reads
    """
    fastq_path = Path(fastq_path)

    if not fastq_path.exists():
        raise ValidationError(f"FASTQ file not found: {fastq_path}")

    try:
        if str(fastq_path).endswith(".gz"):
            file_handle = gzip.open(fastq_path, "rt")
        else:
            file_handle = open(fastq_path)
    except Exception as e:
        raise ValidationError(f"Cannot open FASTQ file {fastq_path}: {e}") from e

    try:
        # Count lines and divide by 4 (FASTQ has 4 lines per read)
        line_count = sum(1 for _ in file_handle)
        read_count = line_count // 4

        if line_count % 4 != 0:
            logging.warning(
                "FASTQ file %s has %d lines (not divisible by 4), may be truncated or malformed",
                Path(fastq_path).name,
                line_count,
            )

        logging.debug(
            "Counted %d reads in %s (%d lines)",
            read_count,
            fastq_path.name,
            line_count,
        )
        return read_count

    except Exception as e:
        raise ValidationError(f"Error counting reads in {fastq_path}: {e}") from e
    finally:
        file_handle.close()


def merge_fastq_files(
    input_files: list[str | Path],
    output_file: str | Path,
    validate_inputs: bool = True,
) -> str:
    """
    Merge multiple FASTQ files into a single output file.

    This is used to combine reads from split-simulation of diploid haplotypes.
    Reads are concatenated in the order provided. Output compression is
    automatically determined from the output filename extension.

    Args:
        input_files: List of input FASTQ files to merge (can be gzipped).
        output_file: Path to merged output file (.gz extension for compression).
        validate_inputs: If True, validate input files before merging.

    Returns:
        Path to the merged output file (as string).

    Raises:
        ValidationError: If any input file is invalid or merge fails.

    Example:
        >>> merged = merge_fastq_files(
        ...     ["hap1_reads.fastq", "hap2_reads.fastq"],
        ...     "merged_reads.fastq.gz"
        ... )
        >>> print(f"Merged to: {merged}")
        Merged to: merged_reads.fastq.gz

    Notes:
        - Input files can be mixed (some gzipped, some not)
        - Output will be gzipped if filename ends with .gz
        - Original files are not modified
        - No de-duplication is performed
    """
    if not input_files:
        raise ValidationError("No input files provided for merging")

    input_files = [Path(f) for f in input_files]
    output_file = Path(output_file)

    # Validate inputs
    if validate_inputs:
        for input_file in input_files:
            validate_fastq(input_file, check_quality=False)

    # Create output directory if needed
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Determine if output should be compressed
    compress_output = str(output_file).endswith(".gz")

    try:
        # Open output file
        out_handle = gzip.open(output_file, "wt") if compress_output else open(output_file, "w")

        total_reads = 0

        # Merge files
        for input_file in input_files:
            # Open input file
            if str(input_file).endswith(".gz"):
                in_handle = gzip.open(input_file, "rt")
            else:
                in_handle = open(input_file)

            # Copy lines
            read_count = 0
            for line in in_handle:
                out_handle.write(line)
                # Count reads (every 4th line is end of a read)
                if line.startswith("@"):
                    read_count += 1

            in_handle.close()

            total_reads += read_count
            logging.info(
                "Merged %d reads from %s",
                read_count,
                Path(input_file).name,
            )

        out_handle.close()

        logging.info(
            "Successfully merged %d FASTQ files (%d total reads) -> %s",
            len(input_files),
            total_reads,
            output_file.name,
        )

        return str(output_file)

    except Exception as e:
        # Clean up partial output file on error
        if output_file.exists():
            output_file.unlink()
        raise ValidationError(f"Failed to merge FASTQ files: {e}") from e
