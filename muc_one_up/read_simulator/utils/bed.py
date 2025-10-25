"""
BED file utilities for genomic region management.

Provides functions for creating and manipulating BED files for
VNTR efficiency modeling and other genomic analyses.
"""

import logging
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)


class BedError(Exception):
    """Exception raised for BED file operation failures."""

    pass


def create_region_bed(
    output_path: Path, chr: str, start: int, end: int, name: str = "region"
) -> Path:
    """
    Create a simple BED file for a genomic region.

    Args:
        output_path: Path to output BED file
        chr: Chromosome name (e.g., "chr1")
        start: Start position (0-based)
        end: End position (exclusive)
        name: Region name

    Returns:
        Path to created BED file

    Raises:
        BedError: If BED creation fails
        ValueError: If coordinates are invalid
    """
    if start < 0:
        raise ValueError(f"Start position must be >= 0, got {start}")
    if end <= start:
        raise ValueError(f"End must be > start, got start={start}, end={end}")

    logger.debug(f"Creating BED file: {output_path}")
    logger.debug(f"  Region: {chr}:{start}-{end} ({name})")

    try:
        with open(output_path, "w") as f:
            f.write(f"{chr}\t{start}\t{end}\t{name}\n")

        logger.debug(f"BED file created: {output_path}")
        return output_path

    except OSError as e:
        raise BedError(f"Failed to create BED file {output_path}: {e}") from e


def create_vntr_bed(output_dir: Path, vntr_region: dict) -> Path:
    """
    Create BED file for VNTR region from config.

    Args:
        output_dir: Directory for BED file
        vntr_region: Dict with keys: chr, start, end, name

    Returns:
        Path to VNTR BED file

    Raises:
        BedError: If creation fails
        ValueError: If vntr_region is invalid
    """
    required_keys = {"chr", "start", "end", "name"}
    if not required_keys.issubset(vntr_region.keys()):
        raise ValueError(f"vntr_region must contain keys: {required_keys}")

    bed_file = output_dir / "vntr_region.bed"

    return create_region_bed(
        bed_file, vntr_region["chr"], vntr_region["start"], vntr_region["end"], vntr_region["name"]
    )


def create_flanking_beds(
    output_dir: Path, vntr_region: dict, flanking_size: int = 10000
) -> tuple[Path, Path, Path]:
    """
    Create BED files for flanking regions around VNTR.

    Args:
        output_dir: Directory for BED files
        vntr_region: Dict with keys: chr, start, end
        flanking_size: Size of flanking regions in bp

    Returns:
        Tuple of (left_bed, right_bed, combined_bed)

    Raises:
        BedError: If creation fails
    """
    chr = vntr_region["chr"]
    vntr_start = vntr_region["start"]
    vntr_end = vntr_region["end"]

    # Left flanking region
    left_start = max(0, vntr_start - flanking_size)
    left_end = vntr_start
    left_bed = output_dir / "flanking_left.bed"

    create_region_bed(left_bed, chr, left_start, left_end, "Flanking_Left")

    # Right flanking region
    right_start = vntr_end
    right_end = vntr_end + flanking_size
    right_bed = output_dir / "flanking_right.bed"

    create_region_bed(right_bed, chr, right_start, right_end, "Flanking_Right")

    # Combined flanking regions
    combined_bed = output_dir / "flanking_combined.bed"

    try:
        with open(combined_bed, "w") as f:
            f.write(f"{chr}\t{left_start}\t{left_end}\tFlanking_Left\n")
            f.write(f"{chr}\t{right_start}\t{right_end}\tFlanking_Right\n")

        logger.debug("Created flanking BED files:")
        logger.debug(f"  Left: {left_bed}")
        logger.debug(f"  Right: {right_bed}")
        logger.debug(f"  Combined: {combined_bed}")

        return left_bed, right_bed, combined_bed

    except OSError as e:
        raise BedError(f"Failed to create combined flanking BED: {e}") from e


def create_non_vntr_bed_from_capture(output_path: Path, capture_bed: Path, vntr_bed: Path) -> Path:
    """
    Create non-VNTR BED by subtracting VNTR from capture targets.

    Uses bedtools subtract to create mutually exclusive regions.

    Args:
        output_path: Path to output non-VNTR BED
        capture_bed: Capture targets BED file
        vntr_bed: VNTR region BED file

    Returns:
        Path to non-VNTR BED file

    Raises:
        BedError: If bedtools subtract fails
    """
    logger.debug("Creating non-VNTR BED from capture targets")
    logger.debug(f"  Capture BED: {capture_bed}")
    logger.debug(f"  VNTR BED: {vntr_bed}")

    try:
        with open(output_path, "w") as out_f:
            subprocess.run(
                ["bedtools", "subtract", "-a", str(capture_bed), "-b", str(vntr_bed)],
                stdout=out_f,
                check=True,
                stderr=subprocess.PIPE,
                text=True,
            )

        logger.debug(f"Non-VNTR BED created: {output_path}")
        return output_path

    except subprocess.CalledProcessError as e:
        raise BedError(f"bedtools subtract failed: {e.stderr}") from e
    except FileNotFoundError as e:
        raise BedError("bedtools not found in PATH. Please install bedtools.") from e


def check_bedtools_available() -> bool:
    """
    Check if bedtools is available in PATH.

    Returns:
        bool: True if bedtools is available
    """
    try:
        subprocess.run(["bedtools", "--version"], capture_output=True, check=True, timeout=5)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        return False
