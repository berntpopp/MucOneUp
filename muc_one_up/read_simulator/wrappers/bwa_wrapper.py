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

    # Run alignment with BWA mem and stream directly to samtools sort
    bwa_cmd = build_tool_command(tools["bwa"], "mem", "-t", threads, human_reference, read1, read2)
    sort_cmd = build_tool_command(tools["samtools"], "sort", "-@", threads, "-o", output_bam, "-")

    logging.info("Running BWA mem | samtools sort pipeline")
    run_pipeline([bwa_cmd, sort_cmd], capture=True, timeout=timeout or 300)
    logging.info("[bwa+samtools] Alignment and BAM sorting completed")

    # Index the BAM file
    index_cmd = build_tool_command(tools["samtools"], "index", output_bam)
    logging.info("[samtools] Indexing BAM: %s", " ".join(index_cmd))
    run_command(index_cmd, capture=True, timeout=60)

    # Check output files exist
    for file in [output_bam, f"{output_bam}.bai"]:
        file_path = Path(file)
        if not file_path.exists() or file_path.stat().st_size == 0:
            raise FileOperationError(f"BWA alignment output file missing or empty: {file}")

    logging.info("Successfully aligned reads to %s", output_bam)
