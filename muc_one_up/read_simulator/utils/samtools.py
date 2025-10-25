"""
Samtools utility wrappers for read simulation pipeline.

Provides consistent interface to samtools operations with proper
error handling and logging.
"""

import logging
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)


class SamtoolsError(Exception):
    """Exception raised for samtools operation failures."""

    pass


def check_samtools_available() -> bool:
    """
    Check if samtools is available in PATH.

    Returns:
        bool: True if samtools is available
    """
    try:
        subprocess.run(["samtools", "--version"], capture_output=True, check=True, timeout=5)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        return False


def extract_reads_by_region(
    input_bam: Path,
    region_bed: Path,
    output_bam: Path,
    exclude_flags: int = 2304,  # Secondary + supplementary
) -> None:
    """
    Extract reads overlapping a BED region using samtools view.

    Args:
        input_bam: Input BAM file
        region_bed: BED file defining region
        output_bam: Output BAM file
        exclude_flags: SAM flags to exclude (default: 2304)

    Raises:
        SamtoolsError: If extraction fails
    """
    logger.debug(f"Extracting reads from {input_bam} using region {region_bed}")

    try:
        subprocess.run(
            [
                "samtools",
                "view",
                "-h",
                "-b",
                "-L",
                str(region_bed),
                "-F",
                str(exclude_flags),
                str(input_bam),
                "-o",
                str(output_bam),
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        logger.debug(f"Reads extracted to {output_bam}")

    except subprocess.CalledProcessError as e:
        raise SamtoolsError(f"Failed to extract reads: {e.stderr}") from e


def downsample_bam(input_bam: Path, output_bam: Path, fraction: float, seed: int = 42) -> None:
    """
    Downsample BAM file using samtools view -s.

    Uses hash-based sampling which preserves read pairs automatically.

    Args:
        input_bam: Input BAM file
        output_bam: Output BAM file
        fraction: Fraction of reads to keep (0.0-1.0)
        seed: Random seed for reproducibility

    Raises:
        SamtoolsError: If downsampling fails
        ValueError: If fraction not in valid range
    """
    if not 0.0 < fraction <= 1.0:
        raise ValueError(f"Fraction must be in (0.0, 1.0], got {fraction}")

    logger.debug(f"Downsampling {input_bam} to {fraction:.3f} (seed={seed})")

    # Format: SEED.FRACTION (e.g., 42.375)
    fraction_str = str(fraction).lstrip("0.")
    downsample_param = f"{seed}.{fraction_str}"

    try:
        subprocess.run(
            [
                "samtools",
                "view",
                "-s",
                downsample_param,
                "-b",
                str(input_bam),
                "-o",
                str(output_bam),
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        logger.debug(f"Downsampled BAM written to {output_bam}")

    except subprocess.CalledProcessError as e:
        raise SamtoolsError(f"Failed to downsample BAM: {e.stderr}") from e


def merge_bams(
    input_bams: list[Path], output_bam: Path, threads: int = 1, force: bool = True
) -> None:
    """
    Merge multiple BAM files using samtools merge.

    Args:
        input_bams: List of input BAM files
        output_bam: Output merged BAM file
        threads: Number of threads for compression
        force: Overwrite output if exists

    Raises:
        SamtoolsError: If merging fails
        ValueError: If less than 2 input BAMs provided
    """
    if len(input_bams) < 2:
        raise ValueError("Need at least 2 BAM files to merge")

    logger.debug(f"Merging {len(input_bams)} BAM files into {output_bam}")

    cmd = ["samtools", "merge"]

    if force:
        cmd.append("-f")

    cmd.extend(["-@", str(threads)])
    cmd.append(str(output_bam))
    cmd.extend([str(bam) for bam in input_bams])

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.debug(f"Merged BAM written to {output_bam}")

    except subprocess.CalledProcessError as e:
        raise SamtoolsError(f"Failed to merge BAMs: {e.stderr}") from e


def sort_bam(input_bam: Path, output_bam: Path, threads: int = 1) -> None:
    """
    Sort BAM file by coordinate using samtools sort.

    Args:
        input_bam: Input BAM file
        output_bam: Output sorted BAM file
        threads: Number of threads for sorting

    Raises:
        SamtoolsError: If sorting fails
    """
    logger.debug(f"Sorting {input_bam}")

    try:
        subprocess.run(
            [
                "samtools",
                "sort",
                "-@",
                str(threads),
                "-o",
                str(output_bam),
                str(input_bam),
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        logger.debug(f"Sorted BAM written to {output_bam}")

    except subprocess.CalledProcessError as e:
        raise SamtoolsError(f"Failed to sort BAM: {e.stderr}") from e


def index_bam(bam_file: Path) -> None:
    """
    Index BAM file using samtools index.

    Args:
        bam_file: BAM file to index

    Raises:
        SamtoolsError: If indexing fails
    """
    logger.debug(f"Indexing {bam_file}")

    try:
        subprocess.run(
            ["samtools", "index", str(bam_file)],
            check=True,
            capture_output=True,
            text=True,
        )

        logger.debug(f"Index created: {bam_file}.bai")

    except subprocess.CalledProcessError as e:
        raise SamtoolsError(f"Failed to index BAM: {e.stderr}") from e


def calculate_mean_coverage(bam_file: Path, region_bed: Path) -> float:
    """
    Calculate mean coverage for a region using samtools depth.

    Args:
        bam_file: BAM file
        region_bed: BED file defining region

    Returns:
        Mean coverage (float), 0.0 if no coverage

    Raises:
        SamtoolsError: If coverage calculation fails
    """
    logger.debug(f"Calculating coverage for {bam_file} in region {region_bed}")

    try:
        result = subprocess.run(
            ["samtools", "depth", "-b", str(region_bed), str(bam_file)],
            check=True,
            capture_output=True,
            text=True,
        )

        if not result.stdout.strip():
            logger.warning(f"No coverage data for {region_bed}")
            return 0.0

        # Parse depth output (chr, pos, depth)
        depths = [int(line.split()[2]) for line in result.stdout.strip().split("\n")]

        mean_cov = sum(depths) / len(depths) if depths else 0.0
        logger.debug(f"Mean coverage: {mean_cov:.2f}x")

        return mean_cov

    except subprocess.CalledProcessError as e:
        raise SamtoolsError(f"Failed to calculate coverage: {e.stderr}") from e


def validate_bam(bam_file: Path) -> bool:
    """
    Validate BAM file using samtools quickcheck.

    Args:
        bam_file: BAM file to validate

    Returns:
        bool: True if valid, False otherwise
    """
    logger.debug(f"Validating {bam_file}")

    try:
        subprocess.run(
            ["samtools", "quickcheck", str(bam_file)],
            check=True,
            capture_output=True,
        )

        logger.debug(f"BAM file is valid: {bam_file}")
        return True

    except subprocess.CalledProcessError:
        logger.error(f"BAM file is corrupted: {bam_file}")
        return False


def get_bam_read_count(bam_file: Path) -> int:
    """
    Get total read count from BAM file using samtools view -c.

    Args:
        bam_file: BAM file

    Returns:
        int: Number of reads

    Raises:
        SamtoolsError: If counting fails
    """
    logger.debug(f"Counting reads in {bam_file}")

    try:
        result = subprocess.run(
            ["samtools", "view", "-c", str(bam_file)],
            check=True,
            capture_output=True,
            text=True,
        )

        count = int(result.stdout.strip())
        logger.debug(f"Read count: {count:,}")

        return count

    except subprocess.CalledProcessError as e:
        raise SamtoolsError(f"Failed to count reads: {e.stderr}") from e
    except ValueError as e:
        raise SamtoolsError(f"Invalid read count output: {result.stdout}") from e
