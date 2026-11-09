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
