"""
Samtools utility wrappers for read simulation pipeline.

Provides consistent interface to samtools operations with proper
error handling and logging.
"""

import logging
from pathlib import Path
from typing import Any

from muc_one_up.read_simulator.command_utils import build_tool_command

from ...exceptions import ExternalToolError
from .common_utils import run_command

logger = logging.getLogger(__name__)


class SamtoolsError(ExternalToolError):
    """Exception raised for samtools operation failures."""

    def __init__(
        self,
        message: str,
        tool: str = "samtools",
        exit_code: int = 1,
        stderr: str = "",
        stdout: str = "",
        cmd: str = "",
    ) -> None:
        super().__init__(tool=tool, exit_code=exit_code, stderr=stderr, stdout=stdout, cmd=cmd)
        self.message = message

    def __str__(self) -> str:
        return self.message


def subsample_args(fraction: float, seed: int) -> list[str]:
    """Build ``samtools view`` arguments for seeded fractional subsampling.

    Uses ``--subsample``/``--subsample-seed`` (samtools >= 1.13) instead of the
    deprecated ``-s SEED.FRAC`` form, which newer samtools versions misread.

    Args:
        fraction: Fraction of templates to keep (0.0-1.0).
        seed: Random seed for reproducible subsampling.

    Returns:
        Argument list to splice into a ``samtools view`` command.
    """
    return ["--subsample", str(fraction), "--subsample-seed", str(seed)]


def run_samtools(
    args: list[str],
    samtools_cmd: str = "samtools",
    **kwargs: Any,
) -> Any:
    """Run a samtools command using build_tool_command."""
    cmd = build_tool_command(samtools_cmd, *args)
    return run_command(cmd, **kwargs)


def check_samtools_available(samtools_cmd: str = "samtools") -> bool:
    """Check if samtools is available in PATH or at specified executable.

    Args:
        samtools_cmd: Samtools executable or command.

    Returns:
        bool: True if samtools is available.
    """
    try:
        run_samtools(["--version"], samtools_cmd=samtools_cmd, capture=True, timeout=5)
        return True
    except (ExternalToolError, FileNotFoundError):
        return False


def extract_reads_by_region(
    input_bam: Path,
    region_bed: Path,
    output_bam: Path,
    exclude_flags: int | None = 2304,  # Secondary + supplementary
    unselected_bam: Path | None = None,
    samtools_cmd: str = "samtools",
) -> None:
    """
    Extract reads overlapping a BED region using samtools view.

    Args:
        input_bam: Input BAM file
        region_bed: BED file defining region
        output_bam: Output BAM file
        exclude_flags: SAM flags to exclude (default: 2304)
        unselected_bam: Optional path for unselected reads
        samtools_cmd: Samtools executable or command (default: "samtools")

    Raises:
        SamtoolsError: If extraction fails
    """
    logger.debug(f"Extracting reads from {input_bam} using region {region_bed}")

    cmd_args = ["view", "-b"]
    if exclude_flags is not None:
        cmd_args.extend(["-F", str(exclude_flags)])
    cmd_args.extend(["-L", str(region_bed), "-o", str(output_bam)])
    if unselected_bam is not None:
        cmd_args.extend(["-U", str(unselected_bam)])
    cmd_args.append(str(input_bam))

    try:
        run_samtools(cmd_args, samtools_cmd=samtools_cmd, capture=True)
        logger.debug(f"Reads extracted to {output_bam}")
    except ExternalToolError as e:
        raise SamtoolsError(f"Failed to extract reads: {e}") from e


def build_downsample_args(
    input_bam: str | Path,
    output_bam: str | Path,
    fraction: float,
    seed: int,
    *,
    threads: int | None = None,
    region_bed: str | Path | None = None,
) -> list[str | int]:
    """Build ``samtools view`` arguments for seeded BAM downsampling."""
    args: list[str | int] = ["view", "-b", *subsample_args(fraction, seed)]
    if region_bed is not None:
        args.extend(["-L", str(region_bed)])
    if threads is not None:
        args.extend(["-@", threads])
    args.extend(["-o", str(output_bam), str(input_bam)])
    return args


def build_merge_args(
    output_bam: str | Path,
    input_bams: list[str] | list[Path],
    *,
    threads: int = 1,
    force: bool = True,
) -> list[str | int]:
    """Build ``samtools merge`` arguments."""
    args: list[str | int] = ["merge"]
    if force:
        args.append("-f")
    args.extend(["-@", threads, str(output_bam), *(str(bam) for bam in input_bams)])
    return args


def build_sort_args(
    input_bam: str | Path,
    output_bam: str | Path,
    *,
    threads: int = 1,
    by_name: bool = False,
) -> list[str | int]:
    """Build ``samtools sort`` arguments."""
    args: list[str | int] = ["sort", "-@", threads]
    if by_name:
        args.append("-n")
    args.extend(["-o", str(output_bam), str(input_bam)])
    return args


def mean_depth_from_lines(lines: Any) -> tuple[float, int]:
    """Compute (mean_coverage, num_positions) from ``samtools depth`` output lines."""
    total_depth = 0
    num_positions = 0
    for line in lines:
        parts = line.strip().split()
        if len(parts) >= 3:
            total_depth += int(parts[2])
            num_positions += 1
    if num_positions == 0:
        return 0.0, 0
    return total_depth / num_positions, num_positions


def downsample_bam(
    input_bam: Path,
    output_bam: Path,
    fraction: float,
    seed: int = 42,
    samtools_cmd: str = "samtools",
) -> None:
    """
    Downsample BAM file using samtools view --subsample.

    Uses hash-based sampling which preserves read pairs automatically.

    Args:
        input_bam: Input BAM file
        output_bam: Output BAM file
        fraction: Fraction of reads to keep (0.0-1.0)
        seed: Random seed for reproducibility
        samtools_cmd: Samtools executable or command (default: "samtools")

    Raises:
        SamtoolsError: If downsampling fails
        ValueError: If fraction not in valid range
    """
    if not 0.0 < fraction <= 1.0:
        raise ValueError(f"Fraction must be in (0.0, 1.0], got {fraction}")

    logger.debug(f"Downsampling {input_bam} to {fraction:.3f} (seed={seed})")

    try:
        run_samtools(
            [str(a) for a in build_downsample_args(input_bam, output_bam, fraction, seed)],
            samtools_cmd=samtools_cmd,
            capture=True,
        )

        logger.debug(f"Downsampled BAM written to {output_bam}")

    except ExternalToolError as e:
        raise SamtoolsError(f"Failed to downsample BAM: {e}") from e


def merge_bams(
    input_bams: list[Path],
    output_bam: Path,
    threads: int = 1,
    force: bool = True,
    samtools_cmd: str = "samtools",
) -> None:
    """
    Merge multiple BAM files using samtools merge.

    Args:
        input_bams: List of input BAM files
        output_bam: Output merged BAM file
        threads: Number of threads for compression
        force: Overwrite output if exists
        samtools_cmd: Samtools executable or command (default: "samtools")

    Raises:
        SamtoolsError: If merging fails
        ValueError: If less than 2 input BAMs provided
    """
    if len(input_bams) < 2:
        raise ValueError("Need at least 2 BAM files to merge")

    logger.debug(f"Merging {len(input_bams)} BAM files into {output_bam}")

    try:
        run_samtools(
            [
                str(a)
                for a in build_merge_args(output_bam, input_bams, threads=threads, force=force)
            ],
            samtools_cmd=samtools_cmd,
            capture=True,
        )
        logger.debug(f"Merged BAM written to {output_bam}")

    except ExternalToolError as e:
        raise SamtoolsError(f"Failed to merge BAMs: {e}") from e


def sort_bam(
    input_bam: Path,
    output_bam: Path,
    threads: int = 1,
    samtools_cmd: str = "samtools",
) -> None:
    """
    Sort BAM file by coordinate using samtools sort.

    Args:
        input_bam: Input BAM file
        output_bam: Output sorted BAM file
        threads: Number of threads for sorting
        samtools_cmd: Samtools executable or command (default: "samtools")

    Raises:
        SamtoolsError: If sorting fails
    """
    logger.debug(f"Sorting {input_bam}")

    try:
        run_samtools(
            [str(a) for a in build_sort_args(input_bam, output_bam, threads=threads)],
            samtools_cmd=samtools_cmd,
            capture=True,
        )

        logger.debug(f"Sorted BAM written to {output_bam}")

    except ExternalToolError as e:
        raise SamtoolsError(f"Failed to sort BAM: {e}") from e


def index_bam(bam_file: Path, samtools_cmd: str = "samtools") -> None:
    """
    Index BAM file using samtools index.

    Args:
        bam_file: BAM file to index
        samtools_cmd: Samtools executable or command (default: "samtools")

    Raises:
        SamtoolsError: If indexing fails
    """
    logger.debug(f"Indexing {bam_file}")

    try:
        run_samtools(
            ["index", str(bam_file)],
            samtools_cmd=samtools_cmd,
            capture=True,
        )

        logger.debug(f"Index created: {bam_file}.bai")

    except ExternalToolError as e:
        raise SamtoolsError(f"Failed to index BAM: {e}") from e


def calculate_mean_coverage(
    bam_file: Path,
    region_bed: Path,
    samtools_cmd: str = "samtools",
) -> float:
    """
    Calculate mean coverage for a region using samtools depth.

    Args:
        bam_file: BAM file
        region_bed: BED file defining region
        samtools_cmd: Samtools executable or command (default: "samtools")

    Returns:
        Mean coverage (float), 0.0 if no coverage

    Raises:
        SamtoolsError: If coverage calculation fails
    """
    logger.debug(f"Calculating coverage for {bam_file} in region {region_bed}")

    try:
        result = run_samtools(
            ["depth", "-a", "-b", str(region_bed), str(bam_file)],
            samtools_cmd=samtools_cmd,
            capture=True,
        )

        stdout = result.stdout or ""
        if not stdout.strip():
            logger.warning(f"No coverage data for {region_bed}")
            return 0.0

        mean_cov, _ = mean_depth_from_lines(stdout.strip().split("\n"))
        logger.debug(f"Mean coverage: {mean_cov:.2f}x")

        return mean_cov

    except ExternalToolError as e:
        raise SamtoolsError(f"Failed to calculate coverage: {e}") from e


def validate_bam(bam_file: Path, samtools_cmd: str = "samtools") -> bool:
    """
    Validate BAM file using samtools quickcheck.

    Args:
        bam_file: BAM file to validate
        samtools_cmd: Samtools executable or command (default: "samtools")

    Returns:
        bool: True if valid, False otherwise
    """
    logger.debug(f"Validating {bam_file}")

    try:
        run_samtools(
            ["quickcheck", str(bam_file)],
            samtools_cmd=samtools_cmd,
            capture=True,
        )

        logger.debug(f"BAM file is valid: {bam_file}")
        return True

    except ExternalToolError:
        logger.error(f"BAM file is corrupted: {bam_file}")
        return False


def get_read_count(bam_file: Path, samtools_cmd: str = "samtools") -> int:
    """
    Get total read count from BAM file using samtools view -c.

    Args:
        bam_file: BAM file
        samtools_cmd: Samtools executable or command (default: "samtools")

    Returns:
        int: Number of reads

    Raises:
        SamtoolsError: If counting fails
    """
    logger.debug(f"Counting reads in {bam_file}")

    try:
        result = run_samtools(
            ["view", "-c", str(bam_file)],
            samtools_cmd=samtools_cmd,
            capture=True,
        )

        stdout = result.stdout
        if stdout is None:
            stdout = "0"
        elif not stdout.strip():
            raise SamtoolsError("Invalid read count output: empty stdout")
        count = int(stdout.strip())
        logger.debug(f"Read count: {count:,}")

        return count

    except ExternalToolError as e:
        raise SamtoolsError(f"Failed to count reads: {e}") from e
    except ValueError as e:
        raise SamtoolsError(f"Invalid read count output: {result.stdout}") from e
