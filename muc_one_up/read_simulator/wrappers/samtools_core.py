#!/usr/bin/env python3
"""
Samtools core operations: reference extraction and merging.

This module provides wrapper functions for fundamental samtools operations:
- extract_subset_reference: Extract reference sequences from BAM files
- merge_bam_files: Merge multiple BAM files into a single BAM
"""

import logging
from pathlib import Path

from ...exceptions import FileOperationError
from ..command_utils import build_tool_command
from ..utils import run_command
from ..utils.samtools import build_merge_args


def extract_subset_reference(sample_bam: str, output_fa: str, tools: dict[str, str]) -> str:
    """
    Extract a subset reference from a BAM file.

    Args:
        sample_bam: Input BAM filename.
        output_fa: Output FASTA filename.
        tools: Dictionary of tool commands.

    Returns:
        Intermediate collated BAM filename.

    Raises:
        SystemExit: If any samtools command fails.
    """
    # Get the base directory and name for intermediate files
    output_path = Path(output_fa)
    output_dir = output_path.parent
    output_base = output_path.name.removesuffix(".fasta").removesuffix(".fa")

    # Create intermediate filenames
    collated_bam = str(output_dir / f"{output_base}_collated.bam")

    # Collate the BAM file (group by read name)
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    cmd = build_tool_command(tools["samtools"], "collate", "-o", collated_bam, sample_bam)
    run_command(cmd, timeout=60, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    # Extract sequences into FASTA format
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    cmd = build_tool_command(tools["samtools"], "fasta", "-n", "-F", "0x900", collated_bam)

    logging.info(f"[samtools] Running: {' '.join(cmd)} > {output_fa}")
    run_command(cmd, timeout=60, stdout_path=Path(output_fa))

    # Verify the output FASTA exists and is non-empty
    if not output_path.exists() or output_path.stat().st_size == 0:
        raise FileOperationError(
            f"Failed to extract subset reference: Output FASTA {output_fa} missing or empty"
        )

    return collated_bam


def merge_bam_files(
    samtools_cmd: str,
    input_bams: list[str],
    output_bam: str,
    threads: int = 4,
    timeout: int = 1800,
) -> str:
    """
    Merge multiple BAM files into a single BAM file.
    """
    if len(input_bams) < 2:
        raise FileOperationError(
            f"merge_bam_files requires at least 2 input BAMs, got {len(input_bams)}"
        )

    for input_bam in input_bams:
        input_bam_path = Path(input_bam)
        if not input_bam_path.exists():
            raise FileOperationError(f"Input BAM file not found: {input_bam}")

    cmd = build_tool_command(
        samtools_cmd,
        *build_merge_args(output_bam, input_bams, threads=threads, force=True),
    )

    logging.info(f"Merging {len(input_bams)} BAM files into: {output_bam}")
    for i, bam in enumerate(input_bams, 1):
        logging.info(f"  Input {i}: {bam}")

    run_command(cmd, timeout=timeout, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    output_bam_path = Path(output_bam)
    if not output_bam_path.exists() or output_bam_path.stat().st_size == 0:
        raise FileOperationError(
            f"Failed to merge BAM files: Output BAM {output_bam} missing or empty"
        )

    logging.info(f"BAM merge complete: {output_bam}")
    logging.info(f"  Output size: {output_bam_path.stat().st_size / (1024 * 1024):.2f} MB")

    return output_bam
