#!/usr/bin/env python3
"""
Wrapper module for BWA operations.

This module provides wrapper functions for BWA operations:
- align_reads: Align paired-end reads using bwa mem
"""

import logging
import os
import sys
from typing import Dict, Optional

from ..utils import run_command


def align_reads(
    read1: str,
    read2: str,
    human_reference: str,
    output_bam: str,
    tools: Dict[str, str],
    threads: int = 4,
    timeout: Optional[int] = 300,
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
        SystemExit: If any alignment or sorting command fails.
    """
    # Check input files exist
    for file in [read1, read2, human_reference]:
        if not os.path.exists(file):
            logging.error(f"Input file not found: {file}")
            sys.exit(1)

    # Get output directory and filename base for intermediate files
    output_dir = os.path.dirname(output_bam)
    output_base = os.path.basename(output_bam).replace(".bam", "")

    # Create intermediate filenames for unsorted BAM
    unsorted_bam = os.path.join(output_dir, f"{output_base}_unsorted.bam")

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
    run_command(
        cmd, timeout=60, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO
    )

    # Index the BAM file
    cmd = [tools["samtools"], "index", output_bam]
    run_command(
        cmd, timeout=60, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO
    )

    # Check output files exist
    for file in [output_bam, f"{output_bam}.bai"]:
        if not os.path.exists(file) or os.path.getsize(file) == 0:
            logging.error(f"Expected output file missing or empty: {file}")
            sys.exit(1)

    # Clean up intermediate files
    if os.path.exists(unsorted_bam):
        try:
            os.remove(unsorted_bam)
            logging.info(f"Removed intermediate file: {unsorted_bam}")
        except Exception as e:
            logging.warning(
                f"Could not remove intermediate file {unsorted_bam}: {str(e)}"
            )

    logging.info(f"Successfully aligned reads to {output_bam}")
