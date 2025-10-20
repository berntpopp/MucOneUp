#!/usr/bin/env python3
"""
PacBio HiFi read simulation pipeline.

This module orchestrates the complete PacBio HiFi read simulation workflow,
integrating pbsim3 for CLR simulation, CCS for consensus generation, and
minimap2 for alignment.

Key Functions:
    simulate_pacbio_hifi_reads: Complete HiFi simulation pipeline

Pipeline Stages:
    1. Multi-pass CLR simulation (pbsim3)
    2. HiFi consensus generation (CCS)
    3. BAM→FASTQ conversion (samtools)
    4. Read alignment (minimap2 with map-hifi preset)
    5. Intermediate file cleanup

Example:
    Simulate PacBio HiFi reads::

        from muc_one_up.read_simulator.pacbio_pipeline import simulate_pacbio_hifi_reads

        final_bam = simulate_pacbio_hifi_reads(
            config={
                "tools": {
                    "pbsim3": "pbsim",
                    "ccs": "ccs",
                    "samtools": "samtools",
                    "minimap2": "minimap2"
                },
                "pacbio_params": {
                    "model_type": "qshmm",
                    "model_file": "/models/QSHMM-SEQUEL.model",
                    "coverage": 30,
                    "pass_num": 3,
                    "min_passes": 3,
                    "min_rq": 0.99,
                    "threads": 8
                }
            },
            input_fa="vntr_diploid.fa",
            output_dir="/output",
            output_base="sample",
            human_reference="/ref/hg38.fa"
        )

Notes:
    - Automatically manages intermediate files
    - Validates all inputs before processing
    - Provides detailed progress logging
    - Handles errors gracefully with cleanup

References:
    - pbsim3: https://github.com/yukiteruono/pbsim3
    - CCS: https://ccs.how/
    - minimap2: https://github.com/lh3/minimap2
"""

import logging
from pathlib import Path
from typing import Any

from ..exceptions import ExternalToolError, FileOperationError
from .constants import MINIMAP2_PRESET_PACBIO_HIFI
from .utils import cleanup_files
from .wrappers.ccs_wrapper import run_ccs_consensus, validate_ccs_parameters
from .wrappers.minimap2_wrapper import align_reads_with_minimap2
from .wrappers.pbsim3_wrapper import run_pbsim3_simulation, validate_pbsim3_parameters
from .wrappers.samtools_wrapper import convert_bam_to_fastq


def simulate_pacbio_hifi_reads(
    config: dict[str, Any],
    input_fa: str,
    human_reference: str | None = None,
) -> str:
    """
    Complete PacBio HiFi read simulation pipeline.

    Orchestrates multi-pass CLR simulation, HiFi consensus generation,
    FASTQ conversion, and alignment to produce publication-ready HiFi reads.

    Args:
        config: Configuration dictionary containing:
                - tools: {
                    "pbsim3": path to pbsim executable,
                    "ccs": path to ccs executable,
                    "samtools": path to samtools executable,
                    "minimap2": path to minimap2 executable
                  }
                - pacbio_params: {
                    "model_type": "qshmm" or "errhmm",
                    "model_file": path to .model file,
                    "coverage": target coverage (e.g., 30),
                    "pass_num": number of passes per molecule (≥2),
                    "min_passes": minimum passes for CCS (≥1),
                    "min_rq": minimum quality for HiFi (0.0-1.0),
                    "threads": number of threads,
                    "seed": random seed (optional),
                    # Optional length/accuracy parameters with defaults
                  }
        input_fa: Input reference FASTA file (e.g., diploid haplotypes).
                  Output files are generated in the same directory with
                  names derived from input filename.
        human_reference: Human reference genome for alignment (optional).
                        If None, alignment step is skipped.

    Returns:
        Path to final output file:
        - If human_reference provided: aligned, sorted, indexed BAM
        - If human_reference is None: HiFi FASTQ file

    Raises:
        RuntimeError: If any pipeline stage fails (wraps underlying exceptions).

    Example:
        Full pipeline with alignment::

            from muc_one_up.read_simulator.pacbio_pipeline import simulate_pacbio_hifi_reads

            final_bam = simulate_pacbio_hifi_reads(
                config={
                    "tools": {
                        "pbsim3": "pbsim",
                        "ccs": "ccs",
                        "samtools": "samtools",
                        "minimap2": "minimap2"
                    },
                    "pacbio_params": {
                        "model_type": "qshmm",
                        "model_file": "/models/QSHMM-SEQUEL.model",
                        "coverage": 30,
                        "pass_num": 3,
                        "min_passes": 3,
                        "min_rq": 0.99,
                        "threads": 8,
                        "seed": 42
                    }
                },
                input_fa="vntr_diploid.fa",
                human_reference="/ref/hg38.fa"
            )

        Simulation only (no alignment)::

            hifi_fastq = simulate_pacbio_hifi_reads(
                config=config,
                input_fa="vntr_diploid.fa",
                human_reference=None  # Skip alignment
            )

    Notes:
        - All intermediate files are automatically cleaned up
        - Progress is logged at INFO level
        - Errors include detailed diagnostic messages
        - Seed parameter enables reproducible simulations
        - If alignment is skipped, HiFi FASTQ is the final output
        - Pipeline is optimized for VNTR and complex genomic regions
    """
    # Extract configuration
    pacbio_params = config.get("pacbio_params", {})
    tools = config.get("tools", {})

    # Extract tool commands
    pbsim3_cmd = tools.get("pbsim3", "pbsim")
    ccs_cmd = tools.get("ccs", "ccs")
    samtools_cmd = tools.get("samtools", "samtools")
    minimap2_cmd = tools.get("minimap2", "minimap2")

    # Extract required parameters
    model_type = pacbio_params.get("model_type")
    model_file = pacbio_params.get("model_file")
    coverage = pacbio_params.get("coverage")
    pass_num = pacbio_params.get("pass_num")
    min_passes = pacbio_params.get("min_passes")
    min_rq = pacbio_params.get("min_rq")
    threads = pacbio_params.get("threads", 4)
    seed = pacbio_params.get("seed")

    # Extract optional length/accuracy parameters (use defaults from constants)
    from .constants import (
        DEFAULT_PBSIM3_ACCURACY_MEAN,
        DEFAULT_PBSIM3_ACCURACY_MIN,
        DEFAULT_PBSIM3_ACCURACY_SD,
        DEFAULT_PBSIM3_LENGTH_MAX,
        DEFAULT_PBSIM3_LENGTH_MEAN,
        DEFAULT_PBSIM3_LENGTH_MIN,
        DEFAULT_PBSIM3_LENGTH_SD,
    )

    accuracy_mean = pacbio_params.get("accuracy_mean", DEFAULT_PBSIM3_ACCURACY_MEAN)
    accuracy_sd = pacbio_params.get("accuracy_sd", DEFAULT_PBSIM3_ACCURACY_SD)
    accuracy_min = pacbio_params.get("accuracy_min", DEFAULT_PBSIM3_ACCURACY_MIN)
    length_mean = pacbio_params.get("length_mean", DEFAULT_PBSIM3_LENGTH_MEAN)
    length_sd = pacbio_params.get("length_sd", DEFAULT_PBSIM3_LENGTH_SD)
    length_min = pacbio_params.get("length_min", DEFAULT_PBSIM3_LENGTH_MIN)
    length_max = pacbio_params.get("length_max", DEFAULT_PBSIM3_LENGTH_MAX)

    # Validate required parameters are present
    required_params = {
        "model_type": model_type,
        "model_file": model_file,
        "coverage": coverage,
        "pass_num": pass_num,
        "min_passes": min_passes,
        "min_rq": min_rq,
    }

    missing_params = [key for key, value in required_params.items() if value is None]
    if missing_params:
        raise RuntimeError(
            f"Missing required PacBio parameters in config: {', '.join(missing_params)}"
        )

    # Validate parameter ranges
    try:
        validate_pbsim3_parameters(
            coverage=coverage,
            pass_num=pass_num,
            accuracy_mean=accuracy_mean,
            accuracy_sd=accuracy_sd,
            accuracy_min=accuracy_min,
            length_mean=length_mean,
            length_min=length_min,
            length_max=length_max,
        )

        validate_ccs_parameters(
            min_passes=min_passes,
            min_rq=min_rq,
        )
    except ValueError as e:
        raise RuntimeError(f"Invalid PacBio parameters: {e}") from e

    # Derive output directory and base name from input file path
    input_path = Path(input_fa)
    output_dir_path = input_path.parent
    output_base = input_path.stem  # e.g., "sample.001.simulated" from "sample.001.simulated.fa"

    # Create output directory (in case it doesn't exist)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    # Define intermediate and final file paths
    clr_prefix = str(output_dir_path / f"{output_base}_clr")
    clr_bam = f"{clr_prefix}.bam"

    hifi_bam = str(output_dir_path / f"{output_base}_hifi.bam")
    hifi_fastq = str(output_dir_path / f"{output_base}_hifi.fastq")
    aligned_bam = str(output_dir_path / f"{output_base}_aligned.bam")

    # Track intermediate files for cleanup
    intermediate_files = []

    try:
        # ==========================================================================
        # STAGE 1: Multi-pass CLR Simulation (pbsim3)
        # ==========================================================================
        logging.info("=" * 80)
        logging.info("STAGE 1: Multi-pass CLR simulation with pbsim3")
        logging.info("=" * 80)

        clr_bam = run_pbsim3_simulation(
            pbsim3_cmd=pbsim3_cmd,
            samtools_cmd=samtools_cmd,
            reference=input_fa,
            model_type=model_type,
            model_file=model_file,
            coverage=coverage,
            output_prefix=clr_prefix,
            pass_num=pass_num,
            accuracy_mean=accuracy_mean,
            accuracy_sd=accuracy_sd,
            accuracy_min=accuracy_min,
            length_mean=length_mean,
            length_sd=length_sd,
            length_min=length_min,
            length_max=length_max,
            seed=seed,
        )

        intermediate_files.append(clr_bam)
        intermediate_files.append(f"{clr_prefix}.maf")  # pbsim3 alignment file
        intermediate_files.append(f"{clr_prefix}.ref")  # pbsim3 reference file

        # ==========================================================================
        # STAGE 2: HiFi Consensus Generation (CCS)
        # ==========================================================================
        logging.info("=" * 80)
        logging.info("STAGE 2: HiFi consensus generation with CCS")
        logging.info("=" * 80)

        hifi_bam = run_ccs_consensus(
            ccs_cmd=ccs_cmd,
            input_bam=clr_bam,
            output_bam=hifi_bam,
            min_passes=min_passes,
            min_rq=min_rq,
            threads=threads,
            seed=seed + 1 if seed is not None else None,  # Offset seed for CCS
        )

        intermediate_files.append(hifi_bam)

        # ==========================================================================
        # STAGE 3: BAM → FASTQ Conversion (samtools)
        # ==========================================================================
        logging.info("=" * 80)
        logging.info("STAGE 3: Converting HiFi BAM to FASTQ")
        logging.info("=" * 80)

        hifi_fastq = convert_bam_to_fastq(
            samtools_cmd=samtools_cmd,
            input_bam=hifi_bam,
            output_fastq=hifi_fastq,
            threads=threads,
        )

        # If no human reference provided, stop here and return FASTQ
        if human_reference is None:
            logging.info("=" * 80)
            logging.info("PacBio HiFi simulation complete (no alignment requested)")
            logging.info(f"Final output: {hifi_fastq}")
            logging.info("=" * 80)

            # Clean up intermediate files (keep HiFi FASTQ as final output)
            cleanup_files(intermediate_files)

            return hifi_fastq

        intermediate_files.append(hifi_fastq)

        # ==========================================================================
        # STAGE 4: Read Alignment (minimap2 with map-hifi preset)
        # ==========================================================================
        logging.info("=" * 80)
        logging.info("STAGE 4: Aligning HiFi reads with minimap2 (map-hifi preset)")
        logging.info("=" * 80)

        aligned_bam = align_reads_with_minimap2(
            minimap2_cmd=minimap2_cmd,
            samtools_cmd=samtools_cmd,
            reference=human_reference,
            reads_fastq=hifi_fastq,
            output_bam=aligned_bam,
            preset=MINIMAP2_PRESET_PACBIO_HIFI,  # Use HiFi-specific preset
            threads=threads,
        )

        # ==========================================================================
        # STAGE 5: Cleanup
        # ==========================================================================
        logging.info("=" * 80)
        logging.info("STAGE 5: Cleaning up intermediate files")
        logging.info("=" * 80)

        cleanup_files(intermediate_files)

        logging.info("=" * 80)
        logging.info("PacBio HiFi simulation pipeline complete!")
        logging.info(f"Final output: {aligned_bam}")
        logging.info("=" * 80)

        return aligned_bam

    except (ExternalToolError, FileOperationError) as e:
        # Wrap tool-specific errors in RuntimeError with context
        logging.error(f"PacBio pipeline failed: {e}")
        raise RuntimeError(
            f"PacBio HiFi simulation pipeline failed during processing.\n"
            f"Error details: {e}\n\n"
            f"Troubleshooting:\n"
            f"  1. Check that all tools are installed: pbsim, ccs, samtools, minimap2\n"
            f"  2. Verify conda environment: conda activate env_pacbio\n"
            f"  3. Check model file exists: {model_file}\n"
            f"  4. Review parameter settings (coverage, pass_num, min_passes, min_rq)\n"
            f"  5. Check log output above for detailed error messages"
        ) from e
    except Exception as e:
        # Catch any unexpected errors
        logging.error(f"Unexpected error in PacBio pipeline: {e}")
        raise RuntimeError(f"PacBio HiFi simulation pipeline failed: {e}") from e
    finally:
        # Always attempt cleanup, even if pipeline fails
        # (Don't let cleanup errors mask the original error)
        try:
            if "intermediate_files" in locals():
                cleanup_files(intermediate_files)
        except Exception as cleanup_error:
            logging.warning(f"Cleanup failed (non-fatal): {cleanup_error}")
