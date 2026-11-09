#!/usr/bin/env python3
"""
Samtools coverage and downsampling operations.

This module provides wrapper functions for samtools coverage and downsampling:
- calculate_vntr_coverage: Calculate coverage over VNTR regions
- calculate_target_coverage: Calculate coverage over target regions
- downsample_bam: Pair-aware downsampling of the templates in a BED region
- downsample_entire_bam: Downsample the entire BAM file
"""

import logging
from pathlib import Path

from ...exceptions import FileOperationError
from ..command_utils import build_tool_command
from ..utils import run_command
from ..utils.samtools import build_downsample_args, mean_depth_from_lines
from .samtools_core import merge_bam_files


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
    run_command(cmd, stdout_path=Path(depth_file))  # no cap: capture BEDs can span an exome
    return _read_mean_coverage(depth_file, "target regions"), depth_file


def _template_names(samtools_exe: str, bam: str, bed: str, names_file: str) -> None:
    """Write the unique names of reads in *bam* that overlap *bed* to *names_file*."""
    sam_file = f"{names_file}.sam"
    run_command(
        build_tool_command(samtools_exe, "view", "-L", bed, bam), stdout_path=Path(sam_file)
    )
    with Path(sam_file).open() as sam:
        names = dict.fromkeys(line.split("\t", 1)[0] for line in sam)
    Path(names_file).write_text("".join(f"{n}\n" for n in names))
    Path(sam_file).unlink()


def downsample_bam(
    samtools_exe: str,
    input_bam: str | Path,
    output_bam: str | Path,
    region_bed: str | Path,
    fraction: float,
    seed: int,
    threads: int,
    keep_bed: str | Path | None = None,
) -> None:
    """Pair-aware downsampling of the templates that touch *region_bed*.

    Reads are split by template name, so a pair with one mate in the region and
    the other outside is treated as one unit: every template with any read
    overlapping *region_bed* is subsampled (``--subsample`` hashes the read
    name, so mates share their fate); all other templates are kept in full, or,
    when *keep_bed* is given, only those with a read overlapping *keep_bed*.
    The coordinate-sorted parts are merged (merge keeps the sort order) and the
    result is indexed.
    """
    if not 0.0 < fraction <= 1.0:
        raise ValueError(f"Fraction must be in (0.0, 1.0], got {fraction}")

    tmp = str(output_bam).removesuffix(".bam")
    names = f"{tmp}_region_names.txt"
    region_bam = f"{tmp}_region_only.bam"
    other_bam = f"{tmp}_other_regions.bam"
    kept_names = f"{tmp}_kept_names.txt"
    kept_bam = f"{tmp}_kept_regions.bam"
    sampled_bam = f"{tmp}_region_sampled.bam"
    input_bam, region_bed = str(input_bam), str(region_bed)

    try:
        _template_names(samtools_exe, input_bam, region_bed, names)
        run_command(
            build_tool_command(
                samtools_exe,
                "view",
                "-b",
                "-@",
                threads,
                "-N",
                names,
                "-o",
                region_bam,
                "-U",
                other_bam,
                input_bam,
            )
        )
        if keep_bed is not None:
            _template_names(samtools_exe, other_bam, str(keep_bed), kept_names)
            run_command(
                build_tool_command(
                    samtools_exe,
                    "view",
                    "-b",
                    "-@",
                    threads,
                    "-N",
                    kept_names,
                    "-o",
                    kept_bam,
                    other_bam,
                )
            )
            other_bam = kept_bam
        run_command(
            build_tool_command(
                samtools_exe,
                *build_downsample_args(region_bam, sampled_bam, fraction, seed, threads=threads),
            )
        )
        merge_bam_files(samtools_exe, [sampled_bam, other_bam], str(output_bam), threads=threads)
        run_command(build_tool_command(samtools_exe, "index", str(output_bam)))
    finally:
        for tmp_file in (
            names,
            region_bam,
            f"{tmp}_other_regions.bam",
            kept_names,
            kept_bam,
            sampled_bam,
        ):
            Path(tmp_file).unlink(missing_ok=True)

    logging.info(
        "Downsampled templates in %s to fraction %.4f (seed: %d)", region_bed, fraction, seed
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
