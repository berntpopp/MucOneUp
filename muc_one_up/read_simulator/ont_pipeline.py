#!/usr/bin/env python3
"""
Oxford Nanopore read simulation pipeline.

This module implements the ONT read simulation pipeline using NanoSim,
with alignment using minimap2. Supports both standard and diploid split-simulation
modes to eliminate allelic bias in diploid references.

Key features:
- Automatic diploid detection and split-simulation
- Configurable coverage correction for training model efficiency
- Reliable timeout handling for all external tool calls
- Error handling with descriptive messages
- Input validation and appropriate parameter defaults
- Consistent output file naming
"""

import logging
from datetime import datetime
from pathlib import Path
from typing import Any

from ..bioinformatics.reference_validation import (
    get_reference_path_for_assembly,
    validate_reference_for_assembly,
)
from ..exceptions import FileOperationError, ValidationError
from .utils import is_diploid_reference, run_split_simulation, write_metadata_file
from .wrappers.nanosim_wrapper import (
    align_ont_reads_with_minimap2,
    run_nanosim_simulation,
)


def simulate_ont_reads_pipeline(
    config: dict[str, Any], input_fa: str, human_reference: str | None = None
) -> str:
    """
    Run the complete Oxford Nanopore read simulation pipeline.

    Supports two modes:
    1. **Diploid Split-Simulation** (default for diploid references):
       - Detects diploid references (2 sequences)
       - Simulates each haplotype independently
       - Merges reads to eliminate length-proportional bias
       - Reduces allelic imbalance from ~3.27:1 to ~1.40:1

    2. **Standard Simulation** (for haploid or when split mode disabled):
       - Simulates reads from combined reference
       - Applies coverage correction if enabled

    Pipeline steps:
    1. Validate input parameters from config
    2. Detect diploid reference (if enable_split_simulation=True)
    3. Run simulation (split or standard mode)
    4. Align reads to the reference using minimap2
    5. Create and index the final BAM file

    Args:
        config: Dictionary containing "tools" and "nanosim_params" sections.
               Required parameters:
               - tools: Dictionary containing "nanosim", "minimap2", "samtools"
               - nanosim_params: Dictionary containing:
                 - training_data_path: Path to NanoSim training model
                 - coverage: Target sequencing coverage
                 Optional parameters:
                 - num_threads: Number of threads to use (default: 4)
                 - min_read_length: Minimum read length
                 - max_read_length: Maximum read length
                 - other_options: Additional NanoSim options
                 - seed: Random seed for reproducibility
                 - correction_factor: Training model efficiency (default: 0.325)
                 - enable_split_simulation: Enable diploid split mode (default: True)
                 - enable_coverage_correction: Apply correction factor (default: True)
        input_fa: Input simulated FASTA file.
        human_reference: Path to human reference genome. If None, uses input_fa.

    Returns:
        Path to the final output BAM file (input_basename_ont.bam).

    Raises:
        ValueError: If required parameters are missing.
        RuntimeError: If any step in the pipeline fails.

    Example:
        >>> config = {
        ...     "tools": {"nanosim": "simulator.py", "minimap2": "minimap2", "samtools": "samtools"},
        ...     "nanosim_params": {
        ...         "training_data_path": "path/to/model",
        ...         "coverage": 200,
        ...         "correction_factor": 0.325,
        ...         "enable_split_simulation": True,
        ...         "seed": 42
        ...     }
        ... }
        >>> bam_file = simulate_ont_reads_pipeline(config, "diploid.fa")
        >>> print(f"Generated: {bam_file}")
    """
    # Start time
    start_time = datetime.now()
    logging.info(
        "Starting ONT read simulation pipeline at %s",
        start_time.strftime("%Y-%m-%d %H:%M:%S"),
    )

    # Extract tool commands and simulation parameters
    tools = config.get("tools", {})
    ns_params = config.get("nanosim_params", {})
    rs_config = config.get("read_simulation", {})

    # Validate required tools
    nanosim_cmd = tools.get("nanosim")
    minimap2_cmd = tools.get("minimap2")
    samtools_cmd = tools.get("samtools")

    if not nanosim_cmd:
        raise ValueError("Missing 'nanosim' command in tools configuration.")
    if not minimap2_cmd:
        raise ValueError("Missing 'minimap2' command in tools configuration.")
    if not samtools_cmd:
        raise ValueError("Missing 'samtools' command in tools configuration.")

    # Validate required simulation parameters
    training_model = ns_params.get("training_data_path")
    coverage = ns_params.get("coverage")

    if not training_model:
        raise ValueError("Missing 'training_data_path' in nanosim_params configuration.")
    if not coverage:
        raise ValueError("Missing 'coverage' in nanosim_params configuration.")

    # Optional parameters with defaults
    threads = ns_params.get("num_threads", rs_config.get("threads", 4))
    min_read_length = ns_params.get("min_read_length")
    max_read_length = ns_params.get("max_read_length")
    other_options = ns_params.get("other_options", "")
    seed = ns_params.get("seed")

    # New diploid split-simulation parameters
    correction_factor = ns_params.get("correction_factor", 0.325)
    enable_split_simulation = ns_params.get("enable_split_simulation", True)
    enable_coverage_correction = ns_params.get("enable_coverage_correction", True)

    # Setup output paths
    input_path = Path(input_fa)
    input_basename = input_path.stem
    output_dir = rs_config.get("output_dir", str(input_path.parent))
    output_prefix = str(Path(output_dir) / f"{input_basename}_ont")

    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Check if reference is diploid
    is_diploid = is_diploid_reference(input_fa)

    # Determine simulation mode
    use_split_simulation = is_diploid and enable_split_simulation

    if is_diploid:
        logging.info("=" * 70)
        logging.info("DIPLOID REFERENCE DETECTED")
        logging.info("=" * 70)
        if use_split_simulation:
            logging.info("Mode: Split-simulation (reduces allelic bias)")
            logging.info("Correction factor: %.3f", correction_factor)
        else:
            logging.info("Mode: Standard simulation (split-simulation disabled)")
        logging.info("=" * 70)
    else:
        logging.info("Haploid reference detected - using standard simulation")
        use_split_simulation = False  # Ensure split mode off for haploid

    # 1. Run simulation (split or standard)
    fastq_file = None

    if use_split_simulation:
        # === DIPLOID SPLIT-SIMULATION MODE ===
        logging.info("1. Starting diploid split-simulation")

        # Create simulation function wrapper for run_split_simulation
        def simulation_wrapper(reference_fasta, output_prefix, coverage, seed, **params):
            """Wrapper to call run_nanosim_simulation with expected signature."""
            return run_nanosim_simulation(
                nanosim_cmd=params["nanosim_cmd"],
                reference_fasta=reference_fasta,
                output_prefix=output_prefix,
                training_model=params["training_model"],
                coverage=coverage,
                threads=params["threads"],
                min_read_length=params.get("min_read_length"),
                max_read_length=params.get("max_read_length"),
                other_options=params.get("other_options", ""),
                seed=seed,
            )

        # Prepare simulation parameters
        sim_params = {
            "nanosim_cmd": nanosim_cmd,
            "training_model": training_model,
            "coverage": float(coverage),
            "threads": int(threads),
            "min_read_length": min_read_length,
            "max_read_length": max_read_length,
            "other_options": other_options,
        }

        # Run split-simulation
        try:
            result = run_split_simulation(
                diploid_fasta=input_fa,
                simulation_func=simulation_wrapper,
                simulation_params=sim_params,
                output_fastq=f"{output_prefix}_merged.fastq",
                correction_factor=correction_factor,
                keep_intermediate=False,
                seed=seed,
            )
            fastq_file = result.merged_fastq
            logging.info("Split-simulation completed successfully")
            logging.info("  Haplotype 1: %d reads", result.reads_hap1)
            logging.info("  Haplotype 2: %d reads", result.reads_hap2)
            logging.info("  Total: %d reads", result.reads_hap1 + result.reads_hap2)
        except Exception as e:
            logging.error("Split-simulation failed: %s", e)
            raise RuntimeError(f"ONT split-simulation failed: {e!s}") from e

    else:
        # === STANDARD SIMULATION MODE ===
        logging.info("1. Starting standard NanoSim simulation")

        # Apply coverage correction if enabled
        actual_coverage = float(coverage)
        if enable_coverage_correction:
            actual_coverage = float(coverage) / correction_factor
            logging.info(
                "Coverage correction enabled: %.1fx → %.1fx (factor=%.3f)",
                float(coverage),
                actual_coverage,
                correction_factor,
            )

        try:
            fastq_file = run_nanosim_simulation(
                nanosim_cmd=nanosim_cmd,
                reference_fasta=input_fa,
                output_prefix=output_prefix,
                training_model=training_model,
                coverage=actual_coverage,
                threads=int(threads),
                min_read_length=min_read_length,
                max_read_length=max_read_length,
                other_options=other_options,
                seed=seed,
            )
            logging.info("NanoSim simulation completed successfully")
        except Exception as e:
            logging.error("NanoSim simulation failed: %s", e)
            raise RuntimeError(f"ONT read simulation failed: {e!s}") from e

    # 2. Align reads with minimap2
    logging.info("2. Starting read alignment with minimap2")
    output_bam = f"{output_prefix}.bam"

    # Use human_reference if provided, otherwise get from config (Issue #28)
    if human_reference is None:
        try:
            # Get assembly from config (default: hg38)
            assembly = config.get("reference_assembly", "hg38")

            # Get reference path for assembly
            ref_path = get_reference_path_for_assembly(config, assembly)
            reference_for_alignment = str(ref_path)

            # Validate reference and indices for minimap2 (logs warnings automatically)
            warnings = validate_reference_for_assembly(config, assembly, aligner="minimap2")

            logging.info("Using reference from config: %s (%s)", reference_for_alignment, assembly)

            # Log index warnings if any
            if warnings:
                for warning in warnings:
                    logging.warning("%s", warning)
        except (ValidationError, FileOperationError) as e:
            logging.warning("Could not load reference from config: %s", e)
            logging.warning("Falling back to aligning against simulated reference")
            reference_for_alignment = input_fa
    else:
        reference_for_alignment = human_reference
        logging.info("Using user-provided reference: %s", reference_for_alignment)

    logging.info("Aligning ONT reads to reference: %s", reference_for_alignment)

    try:
        align_ont_reads_with_minimap2(
            minimap2_cmd=minimap2_cmd,
            samtools_cmd=samtools_cmd,
            human_reference=reference_for_alignment,
            reads_fastq=fastq_file,
            output_bam=output_bam,
            threads=int(threads),
        )
        logging.info("Read alignment completed successfully")
    except Exception as e:
        logging.error("Read alignment failed: %s", e)
        raise RuntimeError(f"ONT read alignment failed: {e!s}") from e

    # Calculate elapsed time
    end_time = datetime.now()
    duration = end_time - start_time
    logging.info(
        "ONT read simulation pipeline completed at %s (duration: %s)",
        end_time.strftime("%Y-%m-%d %H:%M:%S"),
        str(duration).split(".")[0],  # Remove microseconds for cleaner output
    )

    logging.info("Final outputs:")
    logging.info("  Aligned and indexed BAM: %s", output_bam)
    logging.info("  Reads FASTQ: %s", fastq_file)

    # Write metadata file with tool versions and provenance
    metadata_file = write_metadata_file(
        output_dir=output_dir,
        output_base=f"{input_basename}_ont",
        config=config,
        start_time=start_time,
        end_time=end_time,
        platform="ONT",
    )
    logging.info("  Metadata file: %s", metadata_file)

    return output_bam
