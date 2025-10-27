#!/usr/bin/env python3
"""
Wrapper module for samtools operations.

This module provides wrapper functions for samtools operations:
- extract_subset_reference: Extract reference sequences from BAM files
- calculate_vntr_coverage: Calculate coverage over VNTR regions
- calculate_target_coverage: Calculate coverage over target regions
- downsample_bam: Downsample BAM files to target coverage
- convert_sam_to_bam: Convert SAM format to BAM format
- convert_bam_to_fastq: Convert BAM to FASTQ format (single-end)
- convert_bam_to_paired_fastq: Convert BAM to paired FASTQ with validation
- merge_bam_files: Merge multiple BAM files into a single BAM
- sort_and_index_bam: Sort and index BAM files
"""

import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path

from ...exceptions import ExternalToolError, FileOperationError
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
    # SECURITY: Use subprocess.run with file handle redirection, never shell=True
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    cmd = build_tool_command(tools["samtools"], "fasta", "-n", "-F", "0x900", collated_bam)

    try:
        logging.info(f"[samtools] Running: {' '.join(cmd)} > {output_fa}")
        with Path(output_fa).open("w") as fasta_out:
            result = subprocess.run(
                cmd, stdout=fasta_out, stderr=subprocess.PIPE, timeout=60, check=True
            )
        if result.stderr:
            logging.info(f"[samtools] {result.stderr.decode('utf-8', errors='ignore')}")
    except subprocess.CalledProcessError as e:
        raise ExternalToolError(
            tool="samtools fasta",
            exit_code=e.returncode,
            stderr=e.stderr.decode("utf-8", errors="ignore") if e.stderr else "Unknown error",
            cmd=" ".join(cmd),
        ) from e
    except subprocess.TimeoutExpired as e:
        raise ExternalToolError(
            tool="samtools fasta",
            exit_code=-1,
            stderr="Timed out after 60s",
            cmd=" ".join(cmd),
        ) from e

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
    # SECURITY: Use subprocess.run with file handle redirection, never shell=True
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

    try:
        logging.info(f"[samtools] Running: {' '.join(cmd)} > {depth_file}")
        with Path(depth_file).open("w") as depth_out:
            result = subprocess.run(
                cmd, stdout=depth_out, stderr=subprocess.PIPE, timeout=60, check=True
            )
        if result.stderr:
            logging.info(f"[samtools] {result.stderr.decode('utf-8', errors='ignore')}")
    except subprocess.CalledProcessError as e:
        raise ExternalToolError(
            tool="samtools depth",
            exit_code=e.returncode,
            stderr=e.stderr.decode("utf-8", errors="ignore") if e.stderr else "Unknown error",
            cmd=" ".join(cmd),
        ) from e
    except subprocess.TimeoutExpired as e:
        raise ExternalToolError(
            tool="samtools depth",
            exit_code=-1,
            stderr="Timed out after 60s",
            cmd=" ".join(cmd),
        ) from e

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
    # SECURITY: Use subprocess.run with file handle redirection, never shell=True
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

    try:
        logging.info(f"[samtools] Running: {' '.join(cmd)} > {depth_file}")
        with Path(depth_file).open("w") as depth_out:
            result = subprocess.run(
                cmd, stdout=depth_out, stderr=subprocess.PIPE, timeout=60, check=True
            )
        if result.stderr:
            logging.info(f"[samtools] {result.stderr.decode('utf-8', errors='ignore')}")
    except subprocess.CalledProcessError as e:
        raise ExternalToolError(
            tool="samtools depth",
            exit_code=e.returncode,
            stderr=e.stderr.decode("utf-8", errors="ignore") if e.stderr else "Unknown error",
            cmd=" ".join(cmd),
        ) from e
    except subprocess.TimeoutExpired as e:
        raise ExternalToolError(
            tool="samtools depth",
            exit_code=-1,
            stderr="Timed out after 60s",
            cmd=" ".join(cmd),
        ) from e

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

    # Extract reads from the target region, downsample, and save to temporary BAM
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    region_bam = output_bam.replace(".bam", "_region_only.bam")
    cmd = build_tool_command(
        samtools_exe,
        "view",
        "-b",  # Output BAM
        "-s",
        fraction_str,  # Downsampling fraction with seed
        "-@",
        threads,  # Threads (build_tool_command handles conversion)
        input_bam,
        region,
        "-o",
        region_bam,
    )
    run_command(cmd, timeout=60)

    # Extract reads from outside the target region (keep all)
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    other_bam = output_bam.replace(".bam", "_other_regions.bam")
    cmd = build_tool_command(
        samtools_exe,
        "view",
        "-b",  # Output BAM
        "-@",
        threads,  # Threads (build_tool_command handles conversion)
        input_bam,
        f"^{region}",  # Exclude target region
        "-o",
        other_bam,
    )
    run_command(cmd, timeout=60, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO)

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
        threads: Number of threads for samtools operations (default: 4)
        timeout: Timeout in seconds for samtools command (default: 1800)
    """

    output_singleton: str | Path | None = None
    preserve_read_names: bool = True
    validate_pairs: bool = True
    threads: int = 4
    timeout: int = 1800


def convert_bam_to_paired_fastq(
    samtools_cmd: str,
    input_bam: str | Path,
    output_fq1: str | Path,
    output_fq2: str | Path,
    options: FastqConversionOptions | None = None,
) -> tuple[str, str]:
    """Convert paired-end BAM to FASTQ format with integrity validation.

    Extracts paired-end reads from a BAM file to separate R1/R2 FASTQ files.
    Filters secondary/supplementary alignments and validates read pair integrity.

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
        - Output files auto-compress if filenames end with .gz (handled by samtools)
        - Read pair integrity validated by default (R1/R2 must have same read count)
        - Singletons (one mate unmapped) discarded by default
        - Only primary alignments extracted (secondary/supplementary filtered: -F 0x900)

    References:
        samtools fastq: http://www.htslib.org/doc/samtools-fastq.html
    """
    opts = options or FastqConversionOptions()

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
    # Build samtools fastq command
    # Note: samtools auto-compresses if output filename ends with .gz
    cmd_args = [
        "fastq",
        "-1",
        str(output_fq1),  # READ1 output
        "-2",
        str(output_fq2),  # READ2 output
    ]

    # Handle singleton reads
    if opts.output_singleton:
        cmd_args.extend(["-s", str(opts.output_singleton)])
    else:
        # Discard singletons (standard for paired-end workflows)
        cmd_args.extend(["-0", "/dev/null"])

    # Read name formatting
    if opts.preserve_read_names:
        cmd_args.append("-n")  # Preserve original names
    else:
        cmd_args.append("-N")  # Force /1 /2 suffixes

    # Threading
    cmd_args.extend(["-@", str(opts.threads)])

    # Input BAM (must be last)
    cmd_args.append(str(input_bam))

    # Build full command using utility function
    cmd = build_tool_command(samtools_cmd, *cmd_args)

    # ========== EXECUTION ==========
    logging.info(
        f"Converting paired-end BAM to FASTQ: {input_bam_path.name} → "
        f"{output_fq1_path.name}, {output_fq2_path.name}"
    )
    logging.debug(f"  Command: {' '.join(cmd)}")

    # Execute with timeout and error handling
    run_command(
        cmd, timeout=opts.timeout, stderr_prefix="[samtools] ", stderr_log_level=logging.INFO
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
                f"R2 has {r2_reads} reads. This indicates incomplete conversion "
                f"or corrupted output."
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
