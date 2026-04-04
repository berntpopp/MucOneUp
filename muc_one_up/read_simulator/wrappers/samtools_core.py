#!/usr/bin/env python3
"""
Samtools core operations: reference extraction, sorting/indexing, and merging.

This module provides wrapper functions for fundamental samtools operations:
- extract_subset_reference: Extract reference sequences from BAM files
- sort_and_index_bam: Sort and index BAM files
- merge_bam_files: Merge multiple BAM files into a single BAM
"""

import logging
from pathlib import Path

from ...exceptions import FileOperationError
from ..command_utils import build_tool_command
from ..utils import run_command


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
    output_base = output_path.name.replace(".fa", "").replace(".fasta", "")

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


def sort_and_index_bam(
    samtools_exe: str,
    input_bam: str,
    output_bam: str | None = None,
    threads: int = 4,
    by_name: bool = False,
) -> str:
    """
    Sort and index a BAM file.

    Args:
        samtools_exe: Path to the samtools executable.
        input_bam: Input BAM file to sort.
        output_bam: Output sorted BAM file. If None, replaces input_bam with sorted version.
        threads: Number of threads to use.
        by_name: If True, sort by read name instead of coordinate.

    Returns:
        Path to the sorted BAM file.

    Raises:
        SystemExit: If any samtools command fails.
    """
    # Determine output filename
    if output_bam is None:
        output_bam = input_bam

    temp_bam = f"{output_bam}.sorting.bam"

    # Sort the BAM file
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    sort_extra_args = [
        "sort",
        "-@",
        threads,  # build_tool_command handles conversion
    ]

    if by_name:
        sort_extra_args.append("-n")  # Sort by name

    sort_extra_args.extend(["-o", temp_bam, input_bam])

    sort_cmd = build_tool_command(samtools_exe, *sort_extra_args)

    run_command(sort_cmd, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    # Replace original file with sorted version
    temp_bam_path = Path(temp_bam)
    output_bam_path = Path(output_bam)
    if temp_bam_path.exists():
        if output_bam_path.exists():
            output_bam_path.unlink()
        temp_bam_path.rename(output_bam_path)

    # Index the BAM file (only if sorted by coordinate)
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    if not by_name:
        cmd = build_tool_command(samtools_exe, "index", output_bam)
        run_command(cmd, timeout=60, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    return output_bam


def merge_bam_files(
    samtools_cmd: str,
    input_bams: list[str],
    output_bam: str,
    threads: int = 4,
    timeout: int = 1800,
) -> str:
    """
    Merge multiple BAM files into a single BAM file.

    This function is used for combining haplotype-specific BAM files from
    diploid simulations (e.g., pbsim3 creates separate BAMs for each haplotype).

    Args:
        samtools_cmd: Path to the samtools executable.
        input_bams: List of input BAM file paths to merge.
        output_bam: Output merged BAM file path.
        threads: Number of threads to use (default: 4).
        timeout: Timeout in seconds (default: 1800).

    Returns:
        Path to the output merged BAM file.

    Raises:
        ExternalToolError: If samtools merge command fails (propagated from run_command).
        FileOperationError: If input BAMs don't exist or output BAM creation fails.

    Example:
        Merge diploid CLR BAMs from pbsim3::

            from muc_one_up.read_simulator.wrappers.samtools_wrapper import merge_bam_files

            merged_bam = merge_bam_files(
                samtools_cmd="samtools",
                input_bams=["sim_0001.bam", "sim_0002.bam"],
                output_bam="sim_merged.bam",
                threads=8
            )

    Notes:
        - Validates all input BAMs exist before merging
        - Validates output BAM is non-empty after merging
        - Uses build_tool_command for safe command construction
        - Lets ExternalToolError propagate (no catching/re-raising)
        - Requires at least 2 input BAMs (raises error for single BAM)
    """
    # Validate we have at least 2 BAMs to merge
    if len(input_bams) < 2:
        raise FileOperationError(
            f"merge_bam_files requires at least 2 input BAMs, got {len(input_bams)}"
        )

    # Validate all input BAMs exist
    for input_bam in input_bams:
        input_bam_path = Path(input_bam)
        if not input_bam_path.exists():
            raise FileOperationError(f"Input BAM file not found: {input_bam}")

    # Merge BAMs using samtools merge
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    cmd = build_tool_command(
        samtools_cmd,
        "merge",
        "-f",  # Force overwrite if output exists
        "-@",
        threads,  # Threads (build_tool_command handles conversion)
        output_bam,
        *input_bams,  # Unpack list of input BAMs
    )

    logging.info(f"Merging {len(input_bams)} BAM files into: {output_bam}")
    for i, bam in enumerate(input_bams, 1):
        logging.info(f"  Input {i}: {bam}")

    run_command(cmd, timeout=timeout, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    # Validate output BAM exists and is non-empty
    output_bam_path = Path(output_bam)
    if not output_bam_path.exists() or output_bam_path.stat().st_size == 0:
        raise FileOperationError(
            f"Failed to merge BAM files: Output BAM {output_bam} missing or empty"
        )

    logging.info(f"BAM merge complete: {output_bam}")
    logging.info(f"  Output size: {output_bam_path.stat().st_size / (1024 * 1024):.2f} MB")

    return output_bam
