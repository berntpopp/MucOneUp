#!/usr/bin/env python3
"""
Wrapper module for BWA operations.

This module provides wrapper functions for BWA operations:
- align_reads: Align paired-end reads using bwa mem
"""

import logging
from pathlib import Path

from ...exceptions import FileOperationError
from ..command_utils import build_tool_command
from ..utils.common_utils import run_command, run_pipeline


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
        timeout: Timeout in seconds for the alignment pipeline.

    Raises:
        FileOperationError: If input files not found or output files missing
        ExternalToolError: If any external tool fails
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
    # Use run_pipeline to connect bwa mem stdout to samtools view stdin.
    # samtools writes binary BAM directly to file via -o to avoid text decoding issues.
    bwa_cmd = build_tool_command(tools["bwa"], "mem", "-t", threads, human_reference, read1, read2)
    samtools_cmd = build_tool_command(tools["samtools"], "view", "-bS", "-o", unsorted_bam, "-")

    logging.info("Running BWA mem | samtools view pipeline")
    run_pipeline([bwa_cmd, samtools_cmd], capture=True, timeout=timeout or 300)
    logging.info("[bwa+samtools] Alignment and BAM conversion completed")

    # Sort the BAM file
    sort_cmd = build_tool_command(
        tools["samtools"], "sort", "-@", threads, "-o", output_bam, unsorted_bam
    )
    logging.info("[samtools] Sorting BAM: %s", " ".join(sort_cmd))
    run_command(sort_cmd, capture=True, timeout=60)

    # Index the BAM file
    index_cmd = build_tool_command(tools["samtools"], "index", output_bam)
    logging.info("[samtools] Indexing BAM: %s", " ".join(index_cmd))
    run_command(index_cmd, capture=True, timeout=60)

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
            logging.info("Removed intermediate file: %s", unsorted_bam)
        except Exception as e:
            logging.warning("Could not remove intermediate file %s: %s", unsorted_bam, str(e))

    logging.info("Successfully aligned reads to %s", output_bam)
