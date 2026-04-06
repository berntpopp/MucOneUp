#!/usr/bin/env python3
"""
Samtools coverage and downsampling operations.

This module provides wrapper functions for samtools coverage and downsampling:
- calculate_vntr_coverage: Calculate coverage over VNTR regions
- calculate_target_coverage: Calculate coverage over target regions
- downsample_bam: Downsample BAM files to target coverage in a specific region
- downsample_entire_bam: Downsample the entire BAM file
"""

import logging
from pathlib import Path

from ...exceptions import FileOperationError
from ..command_utils import build_tool_command
from ..utils import run_command


def calculate_vntr_coverage(
    samtools_exe: str,
    bam_file: str,
    region: str,
    threads: int,
    output_dir: str,
    output_name: str,
) -> tuple[float, str]:
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
        Tuple of (mean_coverage: float, depth_file: str).

    Raises:
        SystemExit: If the samtools command fails.
    """
    # Create output filename for depth data
    depth_file = str(Path(output_dir) / f"{output_name}_vntr_depth.txt")

    # Run samtools depth to get per-base coverage
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    cmd = build_tool_command(
        samtools_exe,
        "depth",
        "-a",  # Output all positions including zero coverage
        "-r",
        region,  # Target region
        "-@",
        threads,  # Threads (build_tool_command handles conversion)
        bam_file,
    )

    logging.info(f"[samtools] Running: {' '.join(cmd)} > {depth_file}")
    run_command(cmd, timeout=60, stdout_path=Path(depth_file))

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
        return 0.0, depth_file

    mean_coverage = total_depth / num_positions
    logging.info(f"Mean coverage in VNTR region: {mean_coverage:.2f}x")

    return mean_coverage, depth_file


def calculate_target_coverage(
    samtools_exe: str,
    bam_file: str,
    bed_file: str,
    threads: int,
    output_dir: str,
    output_name: str,
) -> tuple[float, str]:
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
        Tuple of (mean_coverage: float, depth_file: str).

    Raises:
        SystemExit: If the samtools command fails or bed file is invalid.
    """
    # Verify bed file exists
    if not Path(bed_file).exists():
        raise FileOperationError(f"BED file not found: {bed_file}")

    # Create output filename for depth data
    depth_file = str(Path(output_dir) / f"{output_name}_target_depth.txt")

    # Run samtools depth with BED file to get per-base coverage
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    cmd = build_tool_command(
        samtools_exe,
        "depth",
        "-a",  # Output all positions including zero coverage
        "-b",
        bed_file,  # Target regions from BED file
        "-@",
        threads,  # Threads (build_tool_command handles conversion)
        bam_file,
    )

    logging.info(f"[samtools] Running: {' '.join(cmd)} > {depth_file}")
    run_command(cmd, timeout=60, stdout_path=Path(depth_file))

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
        return 0.0, depth_file

    mean_coverage = total_depth / num_positions
    logging.info(f"Mean coverage in target regions: {mean_coverage:.2f}x")

    return mean_coverage, depth_file


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

    # Create a temporary BED file for the target region
    region_bed = output_bam.replace(".bam", "_region.bed")
    chrom, coords = region.split(":")
    start, end = coords.split("-")
    with Path(region_bed).open("w") as f:
        f.write(f"{chrom}\t{start}\t{end}\n")

    # Two-pass approach: extract non-region reads (full) then region reads (downsampled).
    # Cannot use -L/-U/-s together because -U captures reads failing ANY filter,
    # so VNTR reads dropped by -s would end up in -U and get merged back.
    region_bam = output_bam.replace(".bam", "_region_only.bam")
    other_bam = output_bam.replace(".bam", "_other_regions.bam")

    # Pass 1: Extract non-region reads at full coverage using -L + -U
    cmd = build_tool_command(
        samtools_exe,
        "view",
        "-b",  # Output BAM
        "-L",
        region_bed,  # Target region BED
        "-U",
        other_bam,  # Non-region reads written here (full coverage)
        "-@",
        threads,
        "-o",
        "/dev/null",  # Discard the region reads from this pass
        input_bam,
    )
    run_command(cmd, timeout=60)

    # Pass 2: Extract and downsample region reads
    cmd = build_tool_command(
        samtools_exe,
        "view",
        "-b",  # Output BAM
        "-s",
        fraction_str,  # Downsampling fraction with seed
        "-L",
        region_bed,  # Target region BED
        "-@",
        threads,
        "-o",
        region_bam,
        input_bam,
    )
    run_command(cmd, timeout=60)

    # Merge the downsampled region BAM with the unchanged "other" BAM
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    cmd = build_tool_command(
        samtools_exe,
        "merge",
        "-f",  # Force overwrite output
        "-@",
        threads,  # Threads (build_tool_command handles conversion)
        output_bam,
        region_bam,
        other_bam,
    )
    run_command(cmd, timeout=60, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    # Index the final BAM
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    cmd = build_tool_command(samtools_exe, "index", output_bam)
    run_command(cmd, timeout=60, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    # Clean up intermediate files
    for tmp_file in [region_bam, other_bam, region_bed]:
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
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    cmd = build_tool_command(
        samtools_exe,
        "view",
        "-b",  # Output BAM
        "-s",
        fraction_str,  # Downsampling fraction with seed
        "-@",
        threads,  # Threads (build_tool_command handles conversion)
        "-o",
        output_bam,
        input_bam,
    )
    run_command(cmd, timeout=60, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    # Index the output BAM
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    cmd = build_tool_command(samtools_exe, "index", output_bam)
    run_command(cmd, timeout=60, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    logging.info(f"Downsampled entire BAM file to fraction {fraction:.4f} (seed: {seed})")
