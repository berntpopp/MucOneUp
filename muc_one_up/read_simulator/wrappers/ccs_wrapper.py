#!/usr/bin/env python3
"""
Wrapper module for CCS (Circular Consensus Sequencing).

This module provides a wrapper for the CCS tool, which generates high-accuracy
HiFi (High-Fidelity) reads from multi-pass PacBio CLR (Continuous Long Read)
sequencing data.

Key Functions:
    run_ccs_consensus: Generate HiFi consensus reads from multi-pass CLR BAM

Key Features:
    - Quality filtering (min_rq threshold)
    - Pass number filtering (min_passes requirement)
    - Multi-threaded processing
    - Empty output detection
    - Validation of input and output files

Example:
    Generate HiFi reads from multi-pass CLR::

        from muc_one_up.read_simulator.wrappers.ccs_wrapper import run_ccs_consensus

        hifi_bam = run_ccs_consensus(
            ccs_cmd="ccs",
            input_bam="multi_pass_clr.bam",
            output_bam="hifi_reads.bam",
            min_passes=3,
            min_rq=0.99,
            threads=8
        )

Notes:
    - min_rq=0.99 corresponds to Q20 accuracy (standard HiFi threshold)
    - min_passes controls trade-off between yield and accuracy
    - Empty output indicates insufficient pass coverage or too strict filtering
    - Output BAM is sorted and can be directly aligned or converted to FASTQ

References:
    - CCS documentation: https://ccs.how/
    - CCS FAQ: https://ccs.how/faq/
    - HiFi sequencing: https://www.pacb.com/technology/hifi-sequencing/
"""

import logging
import math
from pathlib import Path

from ...exceptions import FileOperationError
from ..command_utils import build_tool_command
from ..constants import (
    DEFAULT_CCS_MIN_PASSES,
    DEFAULT_CCS_MIN_RQ,
    DEFAULT_CCS_THREADS,
    DEFAULT_CCS_TIMEOUT,
    MIN_VALID_CCS_OUTPUT_SIZE,
)
from ..utils import run_command


def run_ccs_consensus(
    ccs_cmd: str,
    input_bam: str,
    output_bam: str,
    min_passes: int = DEFAULT_CCS_MIN_PASSES,
    min_rq: float = DEFAULT_CCS_MIN_RQ,
    threads: int = DEFAULT_CCS_THREADS,
    seed: int | None = None,
    timeout: int = DEFAULT_CCS_TIMEOUT,
) -> str:
    """
    Generate HiFi consensus reads from multi-pass CLR BAM using CCS.

    Processes multi-pass PacBio CLR reads to generate high-accuracy HiFi reads
    via circular consensus sequencing. Filters reads based on pass number and
    predicted quality.

    Args:
        ccs_cmd: Path to the ccs executable.
        input_bam: Input BAM file containing multi-pass CLR reads from pbsim3.
        output_bam: Output BAM file for HiFi consensus reads.
        min_passes: Minimum number of passes required for consensus (default: 3).
                   Higher values increase accuracy but reduce yield.
        min_rq: Minimum predicted accuracy (0.0-1.0) for HiFi reads (default: 0.99).
               0.99 = Q20 (99% accuracy, standard HiFi threshold)
               0.999 = Q30 (99.9% accuracy, ultra-high accuracy)
        threads: Number of threads to use (default: 4).
        seed: Random seed (accepted for API compatibility but NOT used by CCS).
             CCS does not support --seed parameter. Reproducibility is achieved
             through consistent input BAM ordering from pbsim3 seed.
        timeout: Timeout in seconds (default: 7200 = 2 hours).

    Returns:
        Path to the output HiFi BAM file.

    Raises:
        FileOperationError: If input BAM doesn't exist or output creation fails.
        ExternalToolError: If ccs command fails (propagated from run_command).

    Example:
        Standard HiFi generation (Q20 threshold)::

            from muc_one_up.read_simulator.wrappers.ccs_wrapper import run_ccs_consensus

            hifi_bam = run_ccs_consensus(
                ccs_cmd="ccs",
                input_bam="sim_clr.bam",
                output_bam="sim_hifi.bam",
                min_passes=3,
                min_rq=0.99,
                threads=8,
                seed=42
            )

        Ultra-high accuracy HiFi (Q30 threshold)::

            hifi_bam = run_ccs_consensus(
                ccs_cmd="ccs",
                input_bam="sim_clr.bam",
                output_bam="sim_hifi_q30.bam",
                min_passes=5,
                min_rq=0.999,
                threads=16
            )

    Notes:
        - Input BAM must contain multi-pass CLR reads (pass_num ≥ 2 from pbsim3)
        - min_passes controls minimum ZMW coverage (≥3 recommended)
        - min_rq filters low-quality consensus reads
        - Empty output may indicate:
          * Insufficient pass coverage (increase pbsim3 pass_num)
          * Too strict filtering (decrease min_passes or min_rq)
          * Invalid input BAM format
        - Output BAM contains HiFi reads ready for alignment
        - CCS automatically sorts output by read name
    """
    # Validate input BAM exists
    input_bam_path = Path(input_bam)
    if not input_bam_path.exists():
        raise FileOperationError(f"Input multi-pass CLR BAM not found: {input_bam}")

    # Validate input BAM is non-empty
    if input_bam_path.stat().st_size == 0:
        raise FileOperationError(f"Input multi-pass CLR BAM is empty: {input_bam}")

    # Build CCS command
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    cmd_args = [
        input_bam,  # Input multi-pass CLR BAM
        output_bam,  # Output HiFi BAM
        "--min-passes",
        min_passes,  # Minimum number of passes
        "--min-rq",
        min_rq,  # Minimum predicted accuracy (0.0-1.0)
        "--num-threads",
        threads,  # Number of threads (build_tool_command handles conversion)
    ]

    # NOTE: CCS does not support --seed parameter for reproducibility
    # The seed parameter is accepted for API compatibility but is not used
    # Reproducibility in CCS is achieved through consistent input BAM ordering

    cmd = build_tool_command(ccs_cmd, *cmd_args)

    # Log CCS parameters for debugging
    logging.info("Running CCS consensus generation:")
    logging.info(f"  Input CLR BAM: {input_bam}")
    logging.info(f"  Output HiFi BAM: {output_bam}")
    logging.info(f"  Min passes: {min_passes}")
    logging.info(f"  Min accuracy (RQ): {min_rq} (Q{int(-10 * math.log10(1 - min_rq))})")
    logging.info(f"  Threads: {threads}")
    if seed is not None:
        logging.info(f"  Note: Seed {seed} provided but CCS does not support --seed parameter")

    # Run CCS consensus generation
    # Let ExternalToolError propagate from run_command (don't catch/re-raise)
    run_command(cmd, timeout=timeout, stderr_prefix="[ccs] ", stderr_log_level=logging.INFO)

    # Validate output BAM exists
    output_bam_path = Path(output_bam)
    if not output_bam_path.exists():
        raise FileOperationError(
            f"CCS consensus generation failed: Expected output BAM {output_bam} not found. "
            f"Check CCS logs for errors."
        )

    # Detect empty or corrupted output (CCS may produce small BAM header even if no reads)
    output_size = output_bam_path.stat().st_size
    if output_size < MIN_VALID_CCS_OUTPUT_SIZE:
        raise FileOperationError(
            f"CCS produced empty or invalid output: {output_bam} ({output_size} bytes). "
            f"This may indicate:\n"
            f"  1. Insufficient pass coverage (increase pbsim3 pass_num)\n"
            f"  2. Too strict filtering (min_passes={min_passes}, min_rq={min_rq})\n"
            f"  3. Invalid input BAM format\n"
            f"Consider reducing min_passes or min_rq, or increasing pbsim3 pass_num."
        )

    logging.info(f"CCS consensus generation complete: {output_bam}")
    logging.info(f"  Output size: {output_size / (1024 * 1024):.2f} MB")

    return output_bam


def validate_ccs_parameters(
    min_passes: int,
    min_rq: float,
) -> None:
    """
    Validate CCS consensus parameters before running.

    Args:
        min_passes: Minimum number of passes required.
        min_rq: Minimum predicted accuracy (0.0-1.0).

    Raises:
        ValueError: If any parameter is outside valid range.

    Example:
        Validate before CCS run::

            from muc_one_up.read_simulator.wrappers.ccs_wrapper import validate_ccs_parameters

            validate_ccs_parameters(min_passes=3, min_rq=0.99)  # OK
            validate_ccs_parameters(min_passes=0, min_rq=0.99)  # Raises ValueError

    Notes:
        - min_passes must be ≥ 1 (typically ≥3 for HiFi)
        - min_rq must be in range [0.0, 1.0]
        - Called automatically by run_ccs_consensus
        - Can be used for early validation in CLI/config parsing
    """
    # Validate min_passes
    if min_passes < 1:
        raise ValueError(f"Invalid min_passes: {min_passes}. Must be ≥ 1")

    if min_passes > 50:
        logging.warning(
            f"Very high min_passes value: {min_passes}. This may result in very low yield."
        )

    # Validate min_rq range
    if not 0.0 <= min_rq <= 1.0:
        raise ValueError(f"Invalid min_rq: {min_rq}. Must be between 0.0 and 1.0")

    # Warn if min_rq is below standard HiFi threshold
    if min_rq < 0.99:
        logging.warning(
            f"min_rq={min_rq} is below standard HiFi threshold (0.99 = Q20). "
            f"Output reads may not meet HiFi quality standards."
        )

    # Warn if min_rq is very high (may reduce yield)
    if min_rq > 0.999:
        logging.warning(
            f"Very high min_rq value: {min_rq} (Q{int(-10 * math.log10(1 - min_rq))}). "
            f"This may result in very low yield."
        )


def estimate_hifi_yield(
    clr_bam_size: int,
    min_passes: int,
    min_rq: float,
) -> float:
    """
    Estimate HiFi yield from multi-pass CLR BAM parameters.

    Provides rough estimate of expected HiFi reads based on CLR BAM size
    and filtering parameters. Useful for validating parameter choices.

    Args:
        clr_bam_size: Size of input multi-pass CLR BAM in bytes.
        min_passes: Minimum number of passes required.
        min_rq: Minimum predicted accuracy threshold.

    Returns:
        Estimated yield as fraction of input (0.0-1.0).

    Example:
        Estimate yield before CCS run::

            from pathlib import Path
            from muc_one_up.read_simulator.wrappers.ccs_wrapper import estimate_hifi_yield

            clr_size = Path("sim_clr.bam").stat().st_size
            estimated_yield = estimate_hifi_yield(clr_size, min_passes=3, min_rq=0.99)
            print(f"Estimated HiFi yield: {estimated_yield * 100:.1f}% of CLR data")

    Notes:
        - Rough estimate only, actual yield varies by data quality
        - Higher min_passes and min_rq reduce yield
        - Typical HiFi yield: 50-80% of CLR data for standard parameters
        - Very strict filtering may yield <10% of input data
    """
    # Base yield estimate (typical CCS conversion rate)
    base_yield = 0.65  # 65% typical yield for standard parameters

    # Adjust for min_passes (higher = more stringent = lower yield)
    pass_penalty = max(0.0, (min_passes - 3) * 0.1)  # -10% per pass above 3

    # Adjust for min_rq (higher = more stringent = lower yield)
    rq_penalty = max(0.0, (min_rq - 0.99) * 5.0)  # -5% per 0.01 increase above 0.99

    # Calculate adjusted yield
    adjusted_yield = base_yield * (1.0 - pass_penalty) * (1.0 - rq_penalty)

    # Clamp to reasonable range [0.05, 0.95]
    adjusted_yield = max(0.05, min(0.95, adjusted_yield))

    return adjusted_yield
