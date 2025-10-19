#!/usr/bin/env python3
"""
Oxford Nanopore read simulation pipeline.

This module implements the ONT read simulation pipeline using NanoSim,
with alignment using minimap2. It provides a standardized interface
for simulating Oxford Nanopore reads from a simulated MUC1 sequence.

Key features:
- Reliable timeout handling for all external tool calls
- Error handling with descriptive messages
- Input validation and appropriate parameter defaults
- Consistent output file naming
"""

import logging
from datetime import datetime
from pathlib import Path
from typing import Any

from .wrappers.nanosim_wrapper import (
    align_ont_reads_with_minimap2,
    run_nanosim_simulation,
)


def simulate_ont_reads_pipeline(
    config: dict[str, Any], input_fa: str, human_reference: str | None = None
) -> str:
    """
    Run the complete Oxford Nanopore read simulation pipeline.

    Pipeline steps:
    1. Validate input parameters from config
    2. Run NanoSim simulation to generate Oxford Nanopore reads
    3. Align reads to the reference using minimap2
    4. Create and index the final BAM file

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
        input_fa: Input simulated FASTA file.
        human_reference: Path to human reference genome. If None, uses input_fa.

    Returns:
        Path to the final output BAM file (input_basename_ont.bam).

    Raises:
        ValueError: If required parameters are missing.
        RuntimeError: If any step in the pipeline fails.
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

    # Setup output paths
    input_path = Path(input_fa)
    input_basename = input_path.stem
    output_dir = rs_config.get("output_dir", str(input_path.parent))
    output_prefix = str(Path(output_dir) / f"{input_basename}_ont")

    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # 1. Run NanoSim simulation
    logging.info("1. Starting NanoSim simulation")
    try:
        fastq_file = run_nanosim_simulation(
            nanosim_cmd=nanosim_cmd,
            reference_fasta=input_fa,
            output_prefix=output_prefix,
            training_model=training_model,
            coverage=float(coverage),
            threads=int(threads),
            min_read_length=min_read_length,
            max_read_length=max_read_length,
            other_options=other_options,
        )
        logging.info("NanoSim simulation completed successfully")
    except Exception as e:
        logging.error("NanoSim simulation failed: %s", e)
        raise RuntimeError(f"ONT read simulation failed: {e!s}") from e

    # 2. Align reads with minimap2
    logging.info("2. Starting read alignment with minimap2")
    output_bam = f"{output_prefix}.bam"

    # Use human_reference if provided, otherwise use input_fa
    reference_for_alignment = human_reference if human_reference else input_fa
    logging.info(f"Aligning ONT reads to reference: {reference_for_alignment}")

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

    return output_bam
