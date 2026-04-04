#!/usr/bin/env python3
"""
Samtools format conversion operations.

This module provides wrapper functions for samtools format conversions:
- convert_sam_to_bam: Convert SAM format to BAM format
- convert_bam_to_fastq: Convert BAM to FASTQ format (single-end)
- convert_bam_to_paired_fastq: Convert BAM to paired FASTQ with validation
- FastqConversionOptions: Dataclass for paired FASTQ conversion options
"""

import logging
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from ...exceptions import ExternalToolError, FileOperationError
from ..command_utils import build_tool_command
from ..utils import run_command, run_pipeline


def convert_sam_to_bam(
    samtools_cmd: str,
    input_sam: str,
    output_bam: str,
    threads: int = 4,
    timeout: int = 1800,
) -> str:
    """
    Convert SAM file to BAM format using samtools view.

    This function provides a reusable SAM→BAM conversion for read simulation
    pipelines (PacBio, ONT, Illumina). Used when tools output SAM format but
    downstream processing requires BAM.

    Args:
        samtools_cmd: Path to the samtools executable.
        input_sam: Input SAM file path.
        output_bam: Output BAM file path.
        threads: Number of threads to use (default: 4).
        timeout: Timeout in seconds (default: 1800).

    Returns:
        Path to the output BAM file.

    Raises:
        ExternalToolError: If samtools view command fails (propagated from run_command).
        FileOperationError: If input SAM doesn't exist or output BAM creation fails.

    Example:
        Convert pbsim3 SAM output to BAM::

            from muc_one_up.read_simulator.wrappers.samtools_wrapper import convert_sam_to_bam

            bam_file = convert_sam_to_bam(
                samtools_cmd="samtools",
                input_sam="simulation.sam",
                output_bam="simulation.bam",
                threads=8
            )

    Notes:
        - Validates input SAM exists before conversion
        - Validates output BAM is non-empty after conversion
        - Uses build_tool_command for safe command construction
        - Lets ExternalToolError propagate (no catching/re-raising)
    """
    # Validate input SAM exists
    input_sam_path = Path(input_sam)
    if not input_sam_path.exists():
        raise FileOperationError(f"Input SAM file not found: {input_sam}")

    # Convert SAM to BAM using samtools view
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    cmd = build_tool_command(
        samtools_cmd,
        "view",
        "-b",  # Output BAM format
        "-@",
        threads,  # Threads (build_tool_command handles conversion)
        "-o",
        output_bam,
        input_sam,
    )

    logging.info(f"Converting SAM to BAM: {input_sam} → {output_bam}")
    run_command(cmd, timeout=timeout, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    # Validate output BAM exists and is non-empty
    output_bam_path = Path(output_bam)
    if not output_bam_path.exists() or output_bam_path.stat().st_size == 0:
        raise FileOperationError(
            f"Failed to convert SAM to BAM: Output BAM {output_bam} missing or empty"
        )

    logging.info(f"SAM→BAM conversion complete: {output_bam}")
    return output_bam


def convert_bam_to_fastq(
    samtools_cmd: str,
    input_bam: str,
    output_fastq: str,
    threads: int = 4,
    timeout: int = 1800,
) -> str:
    """
    Convert BAM file to FASTQ format using samtools fastq.

    This function provides a reusable BAM→FASTQ conversion for read simulation
    pipelines (PacBio, ONT, Illumina). Used when downstream tools require FASTQ
    format for alignment.

    Args:
        samtools_cmd: Path to the samtools executable.
        input_bam: Input BAM file path.
        output_fastq: Output FASTQ file path.
        threads: Number of threads to use (default: 4).
        timeout: Timeout in seconds (default: 1800).

    Returns:
        Path to the output FASTQ file.

    Raises:
        ExternalToolError: If samtools fastq command fails (propagated from run_command).
        FileOperationError: If input BAM doesn't exist or output FASTQ creation fails.

    Example:
        Convert CCS HiFi BAM to FASTQ for alignment::

            from muc_one_up.read_simulator.wrappers.samtools_wrapper import convert_bam_to_fastq

            fastq_file = convert_bam_to_fastq(
                samtools_cmd="samtools",
                input_bam="hifi_reads.bam",
                output_fastq="hifi_reads.fastq",
                threads=8
            )

    Notes:
        - Validates input BAM exists before conversion
        - Validates output FASTQ is non-empty after conversion
        - Uses build_tool_command for safe command construction
        - Lets ExternalToolError propagate (no catching/re-raising)
        - Output is uncompressed FASTQ (compress with gzip if needed)
    """
    # Validate input BAM exists
    input_bam_path = Path(input_bam)
    if not input_bam_path.exists():
        raise FileOperationError(f"Input BAM file not found: {input_bam}")

    # Convert BAM to FASTQ using samtools fastq
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    cmd = build_tool_command(
        samtools_cmd,
        "fastq",
        "-@",
        threads,  # Threads (build_tool_command handles conversion)
        "-0",
        output_fastq,  # Output file for reads without pairing information
        input_bam,
    )

    logging.info(f"Converting BAM to FASTQ: {input_bam} → {output_fastq}")
    run_command(cmd, timeout=timeout, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

    # Validate output FASTQ exists and is non-empty
    output_fastq_path = Path(output_fastq)
    if not output_fastq_path.exists() or output_fastq_path.stat().st_size == 0:
        raise FileOperationError(
            f"Failed to convert BAM to FASTQ: Output FASTQ {output_fastq} missing or empty"
        )

    logging.info(f"BAM→FASTQ conversion complete: {output_fastq}")
    return output_fastq


@dataclass
class FastqConversionOptions:
    """Optional parameters for BAM to FASTQ conversion.

    Attributes:
        output_singleton: Optional path (str or Path) for singleton reads (one mate unmapped).
                         If None, singletons are discarded. (default: None)
        preserve_read_names: If True, preserve original read names without /1 /2 suffixes.
                            If False, add /1 /2 suffixes for legacy tools. (default: True)
        validate_pairs: If True, verify R1 and R2 have same read count.
                       Prevents silent data corruption. (default: True)
        collate_before_conversion: If True, collate BAM by read name before conversion.
                                  Required for proper paired-end handling. (default: True)
        discard_singletons: If True, discard singleton reads (one mate unmapped).
                           Ensures only proper pairs in output. (default: True)
        threads: Number of threads for samtools operations (default: 4)
        timeout: Timeout in seconds for samtools command (default: 1800)
    """

    output_singleton: str | Path | None = None
    preserve_read_names: bool = True
    validate_pairs: bool = True
    collate_before_conversion: bool = True
    discard_singletons: bool = True
    threads: int = 4
    timeout: int = 1800


def convert_bam_to_paired_fastq(
    samtools_cmd: str | Sequence[str],
    input_bam: str | Path,
    output_fq1: str | Path,
    output_fq2: str | Path,
    options: FastqConversionOptions | None = None,
) -> tuple[str, str]:
    """Convert paired-end BAM to FASTQ format with integrity validation.

    Extracts paired-end reads from a BAM file to separate R1/R2 FASTQ files.
    Uses two-step pipeline (collate → fastq) to ensure proper pairing.

    Args:
        samtools_cmd: Path to samtools executable (e.g., "samtools")
        input_bam: Input BAM file path (must exist)
        output_fq1: Output R1 FASTQ path (.gz extension enables auto-compression)
        output_fq2: Output R2 FASTQ path (.gz extension enables auto-compression)
        options: Optional conversion parameters (default: FastqConversionOptions())

    Returns:
        Tuple of (output_fq1, output_fq2) paths as strings

    Raises:
        FileOperationError: Input missing, outputs invalid, or read count mismatch
        ExternalToolError: samtools command failed (propagated from run_command)

    Example:
        Basic usage (most common)::

            fq1, fq2 = convert_bam_to_paired_fastq(
                "samtools",
                "sample.bam",
                "sample_R1.fastq.gz",
                "sample_R2.fastq.gz"
            )

        With custom options::

            opts = FastqConversionOptions(threads=8, validate_pairs=True)
            fq1, fq2 = convert_bam_to_paired_fastq(
                "samtools", "sample.bam", "R1.fq.gz", "R2.fq.gz", opts
            )

    Notes:
        - Uses samtools collate to sort by read name before conversion (default)
        - Singletons (unpaired reads) discarded by default via -s /dev/null
        - Output files auto-compress if filenames end with .gz
        - Read pair integrity validated by default (R1/R2 must have same read count)
        - Only primary alignments extracted (secondary/supplementary filtered: -F 0x900)

    References:
        samtools collate: http://www.htslib.org/doc/samtools-collate.html
        samtools fastq: http://www.htslib.org/doc/samtools-fastq.html
    """
    import shlex
    import subprocess

    opts = options or FastqConversionOptions()

    # Normalize samtools_cmd into both string and list forms:
    # - String form: used for conda wrapper detection and shell commands
    # - List form: used for subprocess.Popen (preserves token boundaries)
    if isinstance(samtools_cmd, str):
        samtools_cmd_str = samtools_cmd
        samtools_cmd_list = shlex.split(samtools_cmd)
    else:
        samtools_cmd_list = list(samtools_cmd)
        samtools_cmd_str = shlex.join(samtools_cmd_list)

    # ========== VALIDATION: Input BAM ==========
    input_bam_path = Path(input_bam)

    if not input_bam_path.exists():
        raise FileOperationError(f"Input BAM file not found: {input_bam}")

    if input_bam_path.stat().st_size == 0:
        raise FileOperationError(f"Input BAM file is empty: {input_bam}")

    # ========== PATH PREPARATION ==========
    output_fq1_path = Path(output_fq1)
    output_fq2_path = Path(output_fq2)

    # Ensure parent directories exist
    output_fq1_path.parent.mkdir(parents=True, exist_ok=True)
    output_fq2_path.parent.mkdir(parents=True, exist_ok=True)

    # ========== COMMAND CONSTRUCTION ==========
    # Detect if samtools_cmd uses conda/mamba wrapping
    use_conda_wrapper = "mamba run" in samtools_cmd_str or "conda run" in samtools_cmd_str

    # Build singleton handling option
    singleton_opt = ""
    if opts.discard_singletons:
        if opts.output_singleton:
            singleton_opt = f"-s {shlex.quote(str(opts.output_singleton))}"
        else:
            singleton_opt = "-s /dev/null"
    elif opts.output_singleton:
        singleton_opt = f"-s {shlex.quote(str(opts.output_singleton))}"

    # Build read name formatting option
    read_name_opt = "-n" if opts.preserve_read_names else "-N"

    # ========== EXECUTION ==========
    if opts.collate_before_conversion:
        # Two-step pipeline: collate | fastq
        logging.info(
            f"Converting paired-end BAM to FASTQ (with collation): "
            f"{input_bam_path.name} → {output_fq1_path.name}, {output_fq2_path.name}"
        )

        if use_conda_wrapper:
            # Use shell=True with single conda/mamba environment for piping
            # This is required because two separate conda/mamba processes
            # cannot communicate via pipe
            pipe_cmd = (
                f"samtools collate -u -O -@ {opts.threads} {shlex.quote(str(input_bam))} | "
                f"samtools fastq "
                f"-1 {shlex.quote(str(output_fq1))} "
                f"-2 {shlex.quote(str(output_fq2))} "
                f"-F 0x900 "
                f"{singleton_opt} "
                f"-0 /dev/null "
                f"{read_name_opt} "
                f"-@ {opts.threads} -"
            )

            # Extract the mamba/conda wrapper prefix (e.g., "mamba run -n env_wessim")
            # and use bash to run the piped command inside the environment
            wrapper_prefix = samtools_cmd_str.rsplit("samtools", 1)[0].strip()
            shell_cmd = f"{wrapper_prefix} bash -c {shlex.quote(pipe_cmd)}"

            logging.debug(f"  Shell command: {shell_cmd}")

            try:
                # shell=True is required for piping within conda/mamba environment.
                # All user inputs are properly quoted with shlex.quote() for safety.
                result = subprocess.run(  # nosec B602
                    shell_cmd,
                    shell=True,
                    check=True,
                    timeout=opts.timeout,
                    capture_output=True,
                    text=True,
                )

                # Log any warnings from samtools
                if result.stderr:
                    stderr_text = result.stderr.strip()
                    if stderr_text:
                        logging.debug(f"  Samtools stderr: {stderr_text}")

            except subprocess.TimeoutExpired as e:
                raise ExternalToolError(
                    tool="samtools collate|fastq",
                    exit_code=-1,
                    stderr=f"Command timed out after {opts.timeout}s",
                    cmd=shell_cmd,
                ) from e
            except subprocess.CalledProcessError as e:
                raise ExternalToolError(
                    tool="samtools collate|fastq",
                    exit_code=e.returncode,
                    stderr=e.stderr if e.stderr else str(e),
                    cmd=shell_cmd,
                ) from e

        else:
            # Use run_pipeline for non-wrapped commands

            collate_cmd = [
                *samtools_cmd_list,
                "collate",
                "-u",  # Uncompressed output for faster piping
                "-O",  # Output to stdout
                "-@",
                str(opts.threads),
                str(input_bam),
            ]

            fastq_cmd = [
                *samtools_cmd_list,
                "fastq",
                "-1",
                str(output_fq1),
                "-2",
                str(output_fq2),
                "-F",
                "0x900",
            ]

            # Add singleton handling
            if singleton_opt:
                fastq_cmd.extend(singleton_opt.split())

            # Add other options
            fastq_cmd.extend(["-0", "/dev/null", read_name_opt, "-@", str(opts.threads), "-"])

            logging.debug(f"  Collate command: {' '.join(collate_cmd)}")
            logging.debug(f"  Fastq command: {' '.join(fastq_cmd)}")

            pipeline_result = run_pipeline(
                [collate_cmd, fastq_cmd],
                capture=True,
                timeout=opts.timeout,
            )

            # Log any warnings from samtools
            if pipeline_result.stderr:
                stderr_text = pipeline_result.stderr.strip()
                if stderr_text:
                    logging.debug(f"  Samtools stderr: {stderr_text}")

    else:
        # Single-step: direct fastq conversion (no collation)
        logging.info(
            f"Converting paired-end BAM to FASTQ (no collation): "
            f"{input_bam_path.name} → {output_fq1_path.name}, {output_fq2_path.name}"
        )

        # Build fastq command for non-collation path

        fastq_cmd = [
            *samtools_cmd_list,
            "fastq",
            "-1",
            str(output_fq1),
            "-2",
            str(output_fq2),
            "-F",
            "0x900",
        ]

        # Add singleton handling
        if singleton_opt:
            fastq_cmd.extend(singleton_opt.split())

        # Add other options
        fastq_cmd.extend(["-0", "/dev/null", read_name_opt, "-@", str(opts.threads)])

        # Input BAM (direct input, no stdin)
        fastq_cmd.append(str(input_bam))

        logging.debug(f"  Fastq command: {' '.join(fastq_cmd)}")

        run_command(
            fastq_cmd,
            timeout=opts.timeout,
            stderr_prefix="[samtools] ",
            stderr_log_level=logging.INFO,
        )

    # ========== VALIDATION: Output Files ==========
    # Check R1 output exists and non-empty
    if not output_fq1_path.exists():
        raise FileOperationError(
            f"Failed to convert BAM to FASTQ: R1 output {output_fq1} not created"
        )

    if output_fq1_path.stat().st_size == 0:
        raise FileOperationError(f"Failed to convert BAM to FASTQ: R1 output {output_fq1} is empty")

    # Check R2 output exists and non-empty
    if not output_fq2_path.exists():
        raise FileOperationError(
            f"Failed to convert BAM to FASTQ: R2 output {output_fq2} not created"
        )

    if output_fq2_path.stat().st_size == 0:
        raise FileOperationError(f"Failed to convert BAM to FASTQ: R2 output {output_fq2} is empty")

    # ========== VALIDATION: Read Pair Integrity ==========
    if opts.validate_pairs:
        logging.info("  Validating read pair integrity...")

        r1_reads = _count_fastq_reads(output_fq1_path)
        r2_reads = _count_fastq_reads(output_fq2_path)

        if r1_reads != r2_reads:
            raise FileOperationError(
                f"Read pair count mismatch: R1 has {r1_reads} reads, "
                f"R2 has {r2_reads} reads. This indicates unpaired reads in the BAM file. "
                f"Enable collate_before_conversion=True and discard_singletons=True to fix."
            )

        logging.info(f"  Validated {r1_reads} read pairs")

    # ========== SUCCESS ==========
    fq1_size_mb = output_fq1_path.stat().st_size / (1024 * 1024)
    fq2_size_mb = output_fq2_path.stat().st_size / (1024 * 1024)

    logging.info(
        f"Paired-end FASTQ conversion complete: R1={fq1_size_mb:.2f} MB, R2={fq2_size_mb:.2f} MB"
    )

    return str(output_fq1), str(output_fq2)


def _count_fastq_reads(fastq_path: Path) -> int:
    """Count reads in FASTQ file (handles gzip compression).

    Args:
        fastq_path: Path to FASTQ file (.fq, .fastq, .fq.gz, .fastq.gz)

    Returns:
        Number of reads in file

    Raises:
        FileOperationError: If file cannot be read or is invalid

    Notes:
        - Automatically detects gzip compression based on file extension
        - Assumes standard FASTQ format (4 lines per read)
        - Fast implementation using line counting
    """
    import gzip

    # Determine opener based on file extension
    opener = gzip.open if str(fastq_path).endswith(".gz") else open

    try:
        with opener(fastq_path, "rt") as f:
            # Count lines and divide by 4 (FASTQ format: header, seq, +, qual)
            line_count = sum(1 for _ in f)
            return line_count // 4
    except Exception as e:
        raise FileOperationError(f"Failed to count reads in {fastq_path}: {e}") from e
