#!/usr/bin/env python3
"""
Wrapper module for samtools operations.

This module provides wrapper functions for samtools operations:
- extract_subset_reference: Extract reference sequences from BAM files
- calculate_vntr_coverage: Calculate coverage over VNTR regions
- calculate_target_coverage: Calculate coverage over target regions
- downsample_bam: Downsample BAM files to target coverage
"""

import logging
from pathlib import Path

from ...exceptions import FileOperationError
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
    cmd = [
        tools["samtools"],
        "collate",
        "-o",
        collated_bam,
        sample_bam,
    ]
    run_command(cmd, timeout=60, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    # Extract sequences into FASTA format
    # Using shell=True and string command to handle redirection properly
    cmd = f"{tools['samtools']} fasta -n -F 0x900 {collated_bam} > {output_fa}"  # type: ignore[assignment]
    run_command(
        cmd,
        shell=True,
        timeout=60,
        stderr_prefix="[samtools] ",
        stderr_log_level=logging.INFO,
    )

    # Verify the output FASTA exists and is non-empty
    if not output_path.exists() or output_path.stat().st_size == 0:
        raise FileOperationError(
            f"Failed to extract subset reference: Output FASTA {output_fa} missing or empty"
        )

    return collated_bam


def calculate_vntr_coverage(
    samtools_exe: str,
    bam_file: str,
    region: str,
    threads: int,
    output_dir: str,
    output_name: str,
) -> float:
    """
    Calculate mean coverage over the VNTR region using samtools depth.

    Args:
        samtools_exe: Path to the samtools executable.
        bam_file: BAM file for which coverage is calculated.
        region: Genomic region in format "chr:start-end".
        threads: Number of threads to use.
        output_dir: Directory for output files.
        output_name: Base name for the coverage output file.

    Returns:
        Mean coverage (float).

    Raises:
        SystemExit: If the samtools command fails.
    """
    # Create output filename for depth data
    depth_file = str(Path(output_dir) / f"{output_name}_vntr_depth.txt")

    # Run samtools depth to get per-base coverage
    cmd = [
        samtools_exe,
        "depth",
        "-a",  # Output all positions including zero coverage
        "-r",
        region,  # Target region
        "-@",
        str(threads),  # Threads
        bam_file,
        ">",
        depth_file,
    ]
    run_command(
        cmd,
        shell=True,
        timeout=60,
        stderr_prefix="[samtools] ",
        stderr_log_level=logging.INFO,
    )

    # Calculate mean coverage from depth file
    total_depth = 0
    num_positions = 0

    try:
        with Path(depth_file).open() as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:  # chr, pos, depth
                    total_depth += int(parts[2])
                    num_positions += 1
    except Exception as e:
        raise FileOperationError(f"Error reading depth file for VNTR coverage: {e}") from e

    # Calculate mean coverage (handle empty file case)
    if num_positions == 0:
        logging.warning("No positions found in VNTR region for coverage calculation")
        return 0.0

    mean_coverage = total_depth / num_positions
    logging.info(f"Mean coverage in VNTR region: {mean_coverage:.2f}x")

    return mean_coverage


def calculate_target_coverage(
    samtools_exe: str,
    bam_file: str,
    bed_file: str,
    threads: int,
    output_dir: str,
    output_name: str,
) -> float:
    """
    Calculate mean coverage over the target (non-VNTR) regions using samtools depth and a BED file.

    Args:
        samtools_exe: Path to the samtools executable.
        bam_file: BAM file for which coverage is calculated.
        bed_file: BED file specifying target regions.
        threads: Number of threads to use.
        output_dir: Directory for output files.
        output_name: Base name for the coverage output file.

    Returns:
        Mean coverage (float).

    Raises:
        SystemExit: If the samtools command fails or bed file is invalid.
    """
    # Verify bed file exists
    if not Path(bed_file).exists():
        raise FileOperationError(f"BED file not found: {bed_file}")

    # Create output filename for depth data
    depth_file = str(Path(output_dir) / f"{output_name}_target_depth.txt")

    # Run samtools depth with BED file to get per-base coverage
    cmd = [
        samtools_exe,
        "depth",
        "-a",  # Output all positions including zero coverage
        "-b",
        bed_file,  # Target regions from BED file
        "-@",
        str(threads),  # Threads
        bam_file,
        ">",
        depth_file,
    ]
    run_command(cmd, shell=True, timeout=60)

    # Calculate mean coverage from depth file
    total_depth = 0
    num_positions = 0

    try:
        with Path(depth_file).open() as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:  # chr, pos, depth
                    total_depth += int(parts[2])
                    num_positions += 1
    except Exception as e:
        raise FileOperationError(f"Error reading depth file for target coverage: {e}") from e

    # Calculate mean coverage (handle empty file case)
    if num_positions == 0:
        logging.warning("No positions found in target regions for coverage calculation")
        return 0.0

    mean_coverage = total_depth / num_positions
    logging.info(f"Mean coverage in target regions: {mean_coverage:.2f}x")

    return mean_coverage


def downsample_bam(
    samtools_exe: str,
    input_bam: str,
    output_bam: str,
    region: str,
    fraction: float,
    seed: int,
    threads: int,
) -> None:
    """
    Downsample the BAM file to the specified fraction in the given region.

    Args:
        samtools_exe: Path to the samtools executable.
        input_bam: Input BAM file.
        output_bam: Output downsampled BAM file.
        region: Genomic region to downsample (format: "chr:start-end").
        fraction: Downsampling fraction (0.0-1.0).
        seed: Random seed for reproducibility.
        threads: Number of threads to use.

    Raises:
        SystemExit: If any samtools command fails.
    """
    # Format the fraction string for samtools view -s option (seed.fraction)
    fraction_str = f"{seed}.{int(fraction * 10000):04d}"

    # Extract reads from the target region, downsample, and save to temporary BAM
    region_bam = output_bam.replace(".bam", "_region_only.bam")
    cmd = [
        samtools_exe,
        "view",
        "-b",  # Output BAM
        "-s",
        fraction_str,  # Downsampling fraction with seed
        "-@",
        str(threads),  # Threads
        input_bam,
        region,
        "-o",
        region_bam,
    ]
    run_command(cmd, timeout=60)

    # Extract reads from outside the target region (keep all)
    other_bam = output_bam.replace(".bam", "_other_regions.bam")
    cmd = [
        samtools_exe,
        "view",
        "-b",  # Output BAM
        "-@",
        str(threads),  # Threads
        input_bam,
        f"^{region}",  # Exclude target region
        "-o",
        other_bam,
    ]
    run_command(cmd, timeout=60, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    # Merge the downsampled region BAM with the unchanged "other" BAM
    cmd = [
        samtools_exe,
        "merge",
        "-f",  # Force overwrite output
        "-@",
        str(threads),  # Threads
        output_bam,
        region_bam,
        other_bam,
    ]
    run_command(cmd, timeout=60, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    # Index the final BAM
    cmd = [samtools_exe, "index", output_bam]
    run_command(cmd, timeout=60, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    # Clean up intermediate files
    for tmp_file in [region_bam, other_bam]:
        tmp_path = Path(tmp_file)
        if tmp_path.exists():
            try:
                tmp_path.unlink()
            except Exception as e:
                logging.warning(f"Could not remove temporary file {tmp_file}: {e!s}")

    logging.info(
        f"Downsampled BAM file in region {region} to fraction {fraction:.4f} (seed: {seed})"
    )


def downsample_entire_bam(
    samtools_exe: str,
    input_bam: str,
    output_bam: str,
    fraction: float,
    seed: int,
    threads: int,
) -> None:
    """
    Downsample the entire BAM file to the specified fraction.

    Args:
        samtools_exe: Path to the samtools executable.
        input_bam: Input BAM file.
        output_bam: Output downsampled BAM file.
        fraction: Downsampling fraction (0.0-1.0).
        seed: Random seed for reproducibility.
        threads: Number of threads to use.

    Raises:
        SystemExit: If any samtools command fails.
    """
    # Format the fraction string for samtools view -s option (seed.fraction)
    fraction_str = f"{seed}.{int(fraction * 10000):04d}"

    # Downsample the entire BAM file
    cmd = [
        samtools_exe,
        "view",
        "-b",  # Output BAM
        "-s",
        fraction_str,  # Downsampling fraction with seed
        "-@",
        str(threads),  # Threads
        "-o",
        output_bam,
        input_bam,
    ]
    run_command(cmd, timeout=60, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    # Index the output BAM
    cmd = [samtools_exe, "index", output_bam]
    run_command(cmd, timeout=60, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    logging.info(f"Downsampled entire BAM file to fraction {fraction:.4f} (seed: {seed})")


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
    sort_cmd = [
        samtools_exe,
        "sort",
        "-@",
        str(threads),
    ]

    if by_name:
        sort_cmd.append("-n")  # Sort by name

    sort_cmd.extend(["-o", temp_bam, input_bam])

    run_command(sort_cmd, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    # Replace original file with sorted version
    temp_bam_path = Path(temp_bam)
    output_bam_path = Path(output_bam)
    if temp_bam_path.exists():
        if output_bam_path.exists():
            output_bam_path.unlink()
        temp_bam_path.rename(output_bam_path)

    # Index the BAM file (only if sorted by coordinate)
    if not by_name:
        cmd = [samtools_exe, "index", output_bam]
        run_command(cmd, timeout=60, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    return output_bam
