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
from ..utils.samtools import (
    build_downsample_args,
    build_merge_args,
    mean_depth_from_lines,
)


def _read_mean_coverage(depth_file: str, label: str) -> float:
    """Read a samtools depth file and return mean coverage."""
    try:
        with Path(depth_file).open() as f:
            mean_coverage, num_positions = mean_depth_from_lines(f)
    except Exception as e:
        raise FileOperationError(f"Error reading depth file for {label} coverage: {e}") from e

    if num_positions == 0:
        logging.warning("No positions found in %s for coverage calculation", label)
        return 0.0

    logging.info("Mean coverage in %s: %.2fx", label, mean_coverage)
    return mean_coverage


def calculate_vntr_coverage(
    samtools_exe: str,
    bam_file: str,
    region: str,
    threads: int,
    output_dir: str,
    output_name: str,
) -> tuple[float, str]:
    """Calculate mean coverage over the VNTR region using samtools depth."""
    depth_file = str(Path(output_dir) / f"{output_name}_vntr_depth.txt")
    cmd = build_tool_command(
        samtools_exe,
        "depth",
        "-a",
        "-r",
        region,
        "-@",
        threads,
        bam_file,
    )
    logging.info(f"[samtools] Running: {' '.join(cmd)} > {depth_file}")
    run_command(cmd, timeout=60, stdout_path=Path(depth_file))
    return _read_mean_coverage(depth_file, "VNTR region"), depth_file


def calculate_target_coverage(
    samtools_exe: str,
    bam_file: str,
    bed_file: str,
    threads: int,
    output_dir: str,
    output_name: str,
) -> tuple[float, str]:
    """Calculate mean coverage over the target (non-VNTR) regions using samtools depth and a BED file."""
    if not Path(bed_file).exists():
        raise FileOperationError(f"BED file not found: {bed_file}")

    depth_file = str(Path(output_dir) / f"{output_name}_target_depth.txt")
    cmd = build_tool_command(
        samtools_exe,
        "depth",
        "-a",
        "-b",
        bed_file,
        "-@",
        threads,
        bam_file,
    )
    logging.info(f"[samtools] Running: {' '.join(cmd)} > {depth_file}")
    run_command(cmd, timeout=60, stdout_path=Path(depth_file))
    return _read_mean_coverage(depth_file, "target regions"), depth_file


def downsample_bam(
    samtools_exe: str,
    input_bam: str,
    output_bam: str,
    region: str,
    fraction: float,
    seed: int,
    threads: int,
) -> None:
    """Downsample the BAM file to the specified fraction in the given region."""
    region_bed = f"{output_bam.removesuffix('.bam')}_region.bed"
    chrom, coords = region.split(":")
    start, end = coords.split("-")
    bed_start = int(start) - 1
    with Path(region_bed).open("w") as f:
        f.write(f"{chrom}\t{bed_start}\t{end}\n")

    region_bam = f"{output_bam.removesuffix('.bam')}_region_only.bam"
    other_bam = f"{output_bam.removesuffix('.bam')}_other_regions.bam"

    # Pass 1: Extract non-region reads at full coverage using -L + -U
    cmd = build_tool_command(
        samtools_exe,
        "view",
        "-b",
        "-L",
        region_bed,
        "-U",
        other_bam,
        "-@",
        threads,
        "-o",
        "/dev/null",
        input_bam,
    )
    run_command(cmd, timeout=60)

    # Pass 2: Extract and downsample region reads
    cmd = build_tool_command(
        samtools_exe,
        *build_downsample_args(
            input_bam, region_bam, fraction, seed, threads=threads, region_bed=region_bed
        ),
    )
    run_command(cmd, timeout=60)

    # Merge the downsampled region BAM with the unchanged "other" BAM
    cmd = build_tool_command(
        samtools_exe,
        *build_merge_args(output_bam, [region_bam, other_bam], threads=threads, force=True),
    )
    run_command(cmd, timeout=60, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    cmd = build_tool_command(samtools_exe, "index", output_bam)
    run_command(cmd, timeout=60, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

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
    """Downsample the entire BAM file to the specified fraction."""
    cmd = build_tool_command(
        samtools_exe,
        *build_downsample_args(input_bam, output_bam, fraction, seed, threads=threads),
    )
    run_command(cmd, timeout=60, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    cmd = build_tool_command(samtools_exe, "index", output_bam)
    run_command(cmd, timeout=60, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    logging.info(f"Downsampled entire BAM file to fraction {fraction:.4f} (seed: {seed})")
