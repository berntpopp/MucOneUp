#!/usr/bin/env python3
"""
Wrapper module for BWA operations.

This module provides wrapper functions for BWA operations:
- align_reads: Align paired-end reads using bwa mem
"""

import logging
from pathlib import Path

from ...exceptions import FileOperationError
from ..utils import run_command


def align_reads(
    read1: str,
    read2: str,
    human_reference: str,
    output_bam: str,
    tools: dict[str, str],
    threads: int = 4,
    timeout: int | None = 300,
) -> None:
    """
    Align paired-end reads with bwa mem, sort with samtools, and index the BAM file.

    Args:
        read1: FASTQ filename for read1.
        read2: FASTQ filename for read2.
        human_reference: Human reference FASTA.
        output_bam: Output BAM filename.
        tools: Dictionary of tool commands.
        threads: Number of threads.

    Raises:
        FileOperationError: If input files not found or output files missing
    """
    # Check input files exist
    for file in [read1, read2, human_reference]:
        if not Path(file).exists():
            raise FileOperationError(f"Input file not found for BWA alignment: {file}")

    # Get output directory and filename base for intermediate files
    output_bam_path = Path(output_bam)
    output_dir = output_bam_path.parent
    output_base = output_bam_path.stem

    # Create intermediate filenames for unsorted BAM
    unsorted_bam = str(output_dir / f"{output_base}_unsorted.bam")

    # Run alignment with BWA mem and pipe to samtools to create BAM
    cmd = (
        f"{tools['bwa']} mem -t {threads} {human_reference} {read1} {read2} | "
        f"{tools['samtools']} view -bS - > {unsorted_bam}"
    )

    run_command(
        cmd,
        shell=True,
        timeout=300,
        stderr_prefix="[bwa+samtools] ",
        stderr_log_level=logging.INFO,
    )

    # Sort the BAM file
    cmd = [
        tools["samtools"],
        "sort",
        "-@",
        str(threads),
        "-o",
        output_bam,
        unsorted_bam,
    ]
    run_command(cmd, timeout=60, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    # Index the BAM file
    cmd = [tools["samtools"], "index", output_bam]
    run_command(cmd, timeout=60, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    # Check output files exist
    for file in [output_bam, f"{output_bam}.bai"]:
        file_path = Path(file)
        if not file_path.exists() or file_path.stat().st_size == 0:
            raise FileOperationError(f"BWA alignment output file missing or empty: {file}")

    # Clean up intermediate files
    unsorted_bam_path = Path(unsorted_bam)
    if unsorted_bam_path.exists():
        try:
            unsorted_bam_path.unlink()
            logging.info(f"Removed intermediate file: {unsorted_bam}")
        except Exception as e:
            logging.warning(f"Could not remove intermediate file {unsorted_bam}: {e!s}")

    logging.info(f"Successfully aligned reads to {output_bam}")
