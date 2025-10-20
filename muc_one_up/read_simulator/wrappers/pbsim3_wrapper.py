#!/usr/bin/env python3
"""
Wrapper module for pbsim3 PacBio read simulation.

This module provides a wrapper for pbsim3 multi-pass CLR (Continuous Long Read)
simulation, which generates raw PacBio reads that can be processed with CCS
to produce HiFi reads.

Key Functions:
    run_pbsim3_simulation: Multi-pass CLR simulation with model-based approach

Key Features:
    - Model-based simulation (QSHMM/ERRHMM)
    - Multi-pass read generation for CCS input
    - Configurable read length and accuracy distributions
    - Automatic SAM→BAM conversion if needed
    - Validation of model files and parameters

Example:
    Simulate multi-pass CLR reads::

        from muc_one_up.read_simulator.wrappers.pbsim3_wrapper import run_pbsim3_simulation

        bam_file = run_pbsim3_simulation(
            pbsim3_cmd="pbsim",
            samtools_cmd="samtools",
            reference="template.fa",
            model_type="qshmm",
            model_file="/path/to/QSHMM-RSII.model",
            coverage=30,
            pass_num=3,
            output_prefix="sim"
        )

Notes:
    - Multi-pass reads are required for CCS consensus generation
    - pass_num ≥2 produces sufficient data for HiFi (≥3 recommended)
    - Model files are technology and chemistry specific
    - Automatic validation prevents silent failures

References:
    - pbsim3: https://github.com/yukiteruono/pbsim3
    - Model files: https://github.com/yukiteruono/pbsim3/tree/master/data
    - CCS: https://ccs.how/
"""

import logging
from pathlib import Path

from ...exceptions import FileOperationError
from ..command_utils import build_tool_command
from ..constants import (
    DEFAULT_PBSIM3_ACCURACY_MEAN,
    DEFAULT_PBSIM3_ACCURACY_MIN,
    DEFAULT_PBSIM3_ACCURACY_SD,
    DEFAULT_PBSIM3_LENGTH_MAX,
    DEFAULT_PBSIM3_LENGTH_MEAN,
    DEFAULT_PBSIM3_LENGTH_MIN,
    DEFAULT_PBSIM3_LENGTH_SD,
    DEFAULT_PBSIM3_PASS_NUM,
    DEFAULT_PBSIM3_TIMEOUT,
    VALID_PBSIM3_MODEL_TYPES,
)
from ..utils import run_command
from .samtools_wrapper import convert_sam_to_bam


def run_pbsim3_simulation(
    pbsim3_cmd: str,
    samtools_cmd: str,
    reference: str,
    model_type: str,
    model_file: str,
    coverage: float,
    output_prefix: str,
    pass_num: int = DEFAULT_PBSIM3_PASS_NUM,
    accuracy_mean: float = DEFAULT_PBSIM3_ACCURACY_MEAN,
    accuracy_sd: float = DEFAULT_PBSIM3_ACCURACY_SD,
    accuracy_min: float = DEFAULT_PBSIM3_ACCURACY_MIN,
    length_mean: int = DEFAULT_PBSIM3_LENGTH_MEAN,
    length_sd: int = DEFAULT_PBSIM3_LENGTH_SD,
    length_min: int = DEFAULT_PBSIM3_LENGTH_MIN,
    length_max: int = DEFAULT_PBSIM3_LENGTH_MAX,
    seed: int | None = None,
    timeout: int = DEFAULT_PBSIM3_TIMEOUT,
) -> list[str]:
    """
    Simulate multi-pass PacBio CLR reads using pbsim3 model-based approach.

    Generates multiple passes of the same molecule (simulating circular consensus
    sequencing) which can be processed with CCS to produce high-accuracy HiFi reads.

    Args:
        pbsim3_cmd: Path to the pbsim3 executable.
        samtools_cmd: Path to the samtools executable (for SAM→BAM conversion).
        reference: Reference sequence FASTA file path.
        model_type: Simulation model type ("qshmm" or "errhmm").
        model_file: Path to pbsim3 model file (.model extension).
        coverage: Target sequencing coverage (e.g., 30 for 30x coverage).
        output_prefix: Prefix for output files (e.g., "sim" → "sim.bam", "sim.ref").
        pass_num: Number of passes per molecule (default: 3, minimum: 2 for HiFi).
        accuracy_mean: Mean per-base accuracy before consensus (default: 0.85).
        accuracy_sd: Standard deviation of accuracy (default: 0.05).
        accuracy_min: Minimum per-base accuracy (default: 0.75).
        length_mean: Mean read length in base pairs (default: 15000).
        length_sd: Standard deviation of read length (default: 5000).
        length_min: Minimum read length in base pairs (default: 5000).
        length_max: Maximum read length in base pairs (default: 30000).
        seed: Random seed for reproducibility (default: None = random).
        timeout: Timeout in seconds (default: 7200 = 2 hours).

    Returns:
        List of paths to output BAM files containing multi-pass CLR reads.
        For diploid references, returns multiple BAMs (one per haplotype).
        For haploid references, returns a single-element list.

    Raises:
        FileOperationError: If model file or reference doesn't exist, or if
                           model type is invalid.
        ExternalToolError: If pbsim3 or samtools command fails (propagated).

    Example:
        Simulate with QSHMM model for Sequel II chemistry::

            from muc_one_up.read_simulator.wrappers.pbsim3_wrapper import run_pbsim3_simulation

            bam_file = run_pbsim3_simulation(
                pbsim3_cmd="pbsim",
                samtools_cmd="samtools",
                reference="vntr_diploid.fa",
                model_type="qshmm",
                model_file="/models/QSHMM-SEQUEL.model",
                coverage=30,
                output_prefix="sim",
                pass_num=3,
                seed=42
            )

        Simulate with custom length distribution::

            bam_file = run_pbsim3_simulation(
                pbsim3_cmd="pbsim",
                samtools_cmd="samtools",
                reference="target_region.fa",
                model_type="qshmm",
                model_file="/models/QSHMM-RSII.model",
                coverage=50,
                output_prefix="deep_sim",
                pass_num=5,
                length_mean=20000,
                length_sd=8000,
                length_min=10000,
                length_max=40000
            )

    Notes:
        - Model files are technology-specific (RSII, Sequel, Sequel II)
        - pass_num ≥2 required for CCS; ≥3 recommended for high-quality HiFi
        - accuracy_* parameters affect pre-consensus error rate
        - Final HiFi accuracy depends on CCS settings (min_passes, min_rq)
        - Seed enables reproducible simulations for benchmarking
        - Output includes: .bam (reads), .maf (alignment), .ref (reference)
        - Automatically detects SAM output and converts to BAM
    """
    # Validate model type against allowed values from constants
    if model_type not in VALID_PBSIM3_MODEL_TYPES:
        valid_types = ", ".join(sorted(VALID_PBSIM3_MODEL_TYPES))
        raise FileOperationError(
            f"Invalid pbsim3 model type: '{model_type}'. Valid options: {valid_types}"
        )

    # Validate model file exists
    model_file_path = Path(model_file)
    if not model_file_path.exists():
        raise FileOperationError(f"pbsim3 model file not found: {model_file}")

    # Validate reference exists
    reference_path = Path(reference)
    if not reference_path.exists():
        raise FileOperationError(f"Reference FASTA file not found: {reference}")

    # Validate pass_num is at least 2 (required for multi-pass)
    if pass_num < 2:
        raise FileOperationError(
            f"Invalid pass_num: {pass_num}. Multi-pass simulation requires pass_num ≥ 2"
        )

    # Build pbsim3 command with all parameters
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    # Note: ERRHMM only supports --accuracy-mean and --length-* parameters
    # QSHMM supports additional --accuracy-sd/--accuracy-min parameters
    cmd_args = [
        "--strategy",
        "wgs",  # Whole-genome shotgun strategy
        "--method",
        model_type,  # qshmm or errhmm
        "--qshmm" if model_type == "qshmm" else "--errhmm",
        model_file,  # Model file path
        "--genome",
        reference,  # Reference FASTA (--genome required for wgs strategy)
        "--depth",
        coverage,  # Target coverage (build_tool_command handles conversion)
        "--pass-num",
        pass_num,  # Number of passes per molecule
        "--accuracy-mean",
        accuracy_mean,  # Mean accuracy
        "--length-mean",
        length_mean,  # Mean read length
        "--length-sd",
        length_sd,  # Read length standard deviation
        "--length-min",
        length_min,  # Minimum read length
        "--length-max",
        length_max,  # Maximum read length
        "--prefix",
        output_prefix,  # Output file prefix
    ]

    # Add seed if specified (for reproducible simulations)
    if seed is not None:
        cmd_args.extend(["--seed", seed])

    cmd = build_tool_command(pbsim3_cmd, *cmd_args)

    # Log simulation parameters for debugging
    logging.info("Running pbsim3 multi-pass CLR simulation:")
    logging.info(f"  Reference: {reference}")
    logging.info(f"  Model: {model_type} ({model_file})")
    logging.info(f"  Coverage: {coverage}x")
    logging.info(f"  Pass number: {pass_num}")
    logging.info(f"  Length: {length_mean} ± {length_sd} bp (min={length_min}, max={length_max})")
    logging.info(f"  Accuracy: {accuracy_mean} ± {accuracy_sd} (min={accuracy_min})")
    logging.info(f"  Output prefix: {output_prefix}")
    if seed is not None:
        logging.info(f"  Seed: {seed}")

    # Run pbsim3 simulation
    # Let ExternalToolError propagate from run_command (don't catch/re-raise)
    run_command(cmd, timeout=timeout, stderr_prefix="[pbsim3] ", stderr_log_level=logging.INFO)

    # Determine output file format and handle multiple BAM files (diploid simulations)
    # pbsim3 can produce:
    # 1. Single file: {prefix}.bam or {prefix}.sam (haploid reference)
    # 2. Multiple files: {prefix}_0001.bam, {prefix}_0002.bam, ... (diploid reference)

    # Check for multiple BAM files (diploid simulation)
    # Pattern: {output_prefix}_*.bam (e.g., sim_0001.bam, sim_0002.bam)
    output_prefix_path = Path(output_prefix)
    multi_bam_pattern = f"{output_prefix_path.name}_*.bam"
    multi_bam_files = sorted(str(p) for p in output_prefix_path.parent.glob(multi_bam_pattern))

    output_bams = []

    if multi_bam_files:
        # Multiple BAM files detected (diploid/polyploid reference)
        logging.info(f"pbsim3 produced {len(multi_bam_files)} BAM files (diploid simulation)")

        # Validate each BAM file exists and is non-empty
        for i, bam_file in enumerate(multi_bam_files, 1):
            bam_path = Path(bam_file)
            if not bam_path.exists():
                raise FileOperationError(f"pbsim3 output BAM {bam_file} not found after simulation")
            if bam_path.stat().st_size == 0:
                raise FileOperationError(f"pbsim3 produced empty BAM file: {bam_file}")

            logging.info(
                f"  Haplotype {i} CLR BAM: {bam_file} ({bam_path.stat().st_size / (1024 * 1024):.2f} MB)"
            )
            output_bams.append(bam_file)

    else:
        # Single haploid reference - check for .bam or .sam output
        output_bam = f"{output_prefix}.bam"
        output_sam = f"{output_prefix}.sam"

        output_bam_path = Path(output_bam)
        output_sam_path = Path(output_sam)

        if output_sam_path.exists() and not output_bam_path.exists():
            # Single SAM file detected (needs conversion)
            logging.info(
                f"pbsim3 produced SAM output, converting to BAM: {output_sam} → {output_bam}"
            )

            # Convert SAM to BAM using reusable samtools wrapper
            convert_sam_to_bam(
                samtools_cmd=samtools_cmd,
                input_sam=output_sam,
                output_bam=output_bam,
                threads=4,  # Use moderate threads for conversion
                timeout=timeout,
            )

            # Clean up intermediate SAM file
            try:
                if output_sam_path.exists():
                    output_sam_path.unlink()
                    logging.debug(f"Removed intermediate SAM file: {output_sam}")
            except Exception as e:
                logging.warning(f"Could not remove intermediate SAM file {output_sam}: {e}")

        # Validate single BAM exists
        if not output_bam_path.exists():
            raise FileOperationError(
                f"pbsim3 simulation failed: Expected output BAM {output_bam} not found. "
                f"Also checked for diploid BAMs matching pattern: {multi_bam_pattern}. "
                f"Check pbsim3 logs for errors."
            )

        if output_bam_path.stat().st_size == 0:
            raise FileOperationError(
                f"pbsim3 simulation produced empty BAM file: {output_bam}. "
                f"This may indicate insufficient reference length or coverage settings."
            )

        logging.info(f"pbsim3 simulation complete: {output_bam}")
        logging.info(f"  Output size: {output_bam_path.stat().st_size / (1024 * 1024):.2f} MB")

        output_bams.append(output_bam)

    # Log total output statistics
    total_size = sum(Path(bam).stat().st_size for bam in output_bams)
    logging.info(f"pbsim3 simulation complete: {len(output_bams)} CLR BAM file(s)")
    logging.info(f"  Total size: {total_size / (1024 * 1024):.2f} MB")

    return output_bams


def validate_pbsim3_parameters(
    coverage: float,
    pass_num: int,
    accuracy_mean: float,
    accuracy_sd: float,
    accuracy_min: float,
    length_mean: int,
    length_min: int,
    length_max: int,
) -> None:
    """
    Validate pbsim3 simulation parameters before running.

    Args:
        coverage: Target sequencing coverage.
        pass_num: Number of passes per molecule.
        accuracy_mean: Mean per-base accuracy.
        accuracy_sd: Standard deviation of accuracy.
        accuracy_min: Minimum per-base accuracy.
        length_mean: Mean read length.
        length_min: Minimum read length.
        length_max: Maximum read length.

    Raises:
        ValueError: If any parameter is outside valid range.

    Example:
        Validate before simulation::

            from muc_one_up.read_simulator.wrappers.pbsim3_wrapper import validate_pbsim3_parameters

            validate_pbsim3_parameters(
                coverage=30,
                pass_num=3,
                accuracy_mean=0.85,
                accuracy_sd=0.05,
                accuracy_min=0.75,
                length_mean=15000,
                length_min=5000,
                length_max=30000
            )

    Notes:
        - Called automatically by run_pbsim3_simulation
        - Can be used for early validation in CLI/config parsing
        - Prevents cryptic pbsim3 errors from invalid parameters
    """
    from ..constants import MAX_COVERAGE, MAX_PASS_NUM, MIN_COVERAGE, MIN_PASS_NUM

    # Validate coverage range
    if coverage < MIN_COVERAGE or coverage > MAX_COVERAGE:
        raise ValueError(
            f"Invalid coverage: {coverage}. Must be between {MIN_COVERAGE} and {MAX_COVERAGE}"
        )

    # Validate pass_num range
    if pass_num < MIN_PASS_NUM or pass_num > MAX_PASS_NUM:
        raise ValueError(
            f"Invalid pass_num: {pass_num}. Must be between {MIN_PASS_NUM} and {MAX_PASS_NUM}"
        )

    # Validate accuracy parameters (0.0-1.0)
    if not 0.0 <= accuracy_mean <= 1.0:
        raise ValueError(f"Invalid accuracy_mean: {accuracy_mean}. Must be between 0.0 and 1.0")

    if not 0.0 <= accuracy_sd <= 1.0:
        raise ValueError(f"Invalid accuracy_sd: {accuracy_sd}. Must be between 0.0 and 1.0")

    if not 0.0 <= accuracy_min <= 1.0:
        raise ValueError(f"Invalid accuracy_min: {accuracy_min}. Must be between 0.0 and 1.0")

    # Validate accuracy_min <= accuracy_mean
    if accuracy_min > accuracy_mean:
        raise ValueError(
            f"Invalid accuracy parameters: accuracy_min ({accuracy_min}) > "
            f"accuracy_mean ({accuracy_mean})"
        )

    # Validate length parameters
    if length_min < 1:
        raise ValueError(f"Invalid length_min: {length_min}. Must be ≥ 1")

    if length_max < length_min:
        raise ValueError(
            f"Invalid length parameters: length_max ({length_max}) < length_min ({length_min})"
        )

    if length_mean < length_min or length_mean > length_max:
        raise ValueError(
            f"Invalid length_mean: {length_mean}. Must be between length_min ({length_min}) "
            f"and length_max ({length_max})"
        )
