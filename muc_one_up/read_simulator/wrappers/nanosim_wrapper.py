#!/usr/bin/env python3
"""
NanoSim wrapper for Oxford Nanopore read simulation.

This module provides wrapper functions for NanoSim, a nanopore sequencing read
simulator. It enables integration with the MucOneUp pipeline and handles proper
command construction, execution, and error handling.
"""

import logging
from pathlib import Path

from ..command_utils import build_tool_command
from ..utils import run_command


def run_nanosim_simulation(
    nanosim_cmd: str,
    reference_fasta: str,
    output_prefix: str,
    training_model: str,
    coverage: float,
    threads: int = 4,
    min_read_length: int | None = None,
    max_read_length: int | None = None,
    other_options: str = "",
    timeout: int = 3600,
    seed: int | None = None,
) -> str:
    """
    Run NanoSim simulation to generate Oxford Nanopore reads.

    Args:
        nanosim_cmd: Path to the NanoSim simulator.py command
        reference_fasta: Path to reference FASTA file
        output_prefix: Prefix for output files
        training_model: Path to NanoSim training model
        coverage: Desired coverage
        threads: Number of threads to use (default: 4)
        min_read_length: Minimum read length (optional)
        max_read_length: Maximum read length (optional)
        other_options: Additional NanoSim options (optional)
        timeout: Timeout in seconds (default: 3600)
        seed: Random seed for reproducibility (optional)

    Returns:
        Path to the generated FASTQ file

    Raises:
        RuntimeError: If NanoSim execution fails
    """
    # Ensure output directory exists
    output_path = Path(output_prefix).resolve()
    output_dir = output_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    # Build the command with required parameters
    # SECURITY: Always use list form, never shell=True
    # Use centralized build_tool_command to safely handle multi-word commands (conda/mamba)
    cmd_list = build_tool_command(
        nanosim_cmd,
        "genome",
        "-rg",
        reference_fasta,
        "-c",
        training_model,
        "-o",
        output_prefix,
        "-t",
        threads,  # build_tool_command handles conversion
        "-x",
        coverage,  # build_tool_command handles conversion
    )

    # Add seed if provided
    if seed is not None:
        cmd_list.extend(["--seed", str(seed)])
        logging.info(f"[NanoSim] Using random seed: {seed}")

    # Add optional parameters
    if min_read_length:
        cmd_list.extend(["--min_len", str(min_read_length)])
    if max_read_length:
        cmd_list.extend(["--max_len", str(max_read_length)])
    if other_options:
        # Note: other_options may need splitting if it contains spaces
        # For now, assuming it's already a properly formatted string or list
        if isinstance(other_options, str) and " " in other_options:
            # Build a temporary command to split other_options safely
            other_parts = build_tool_command(other_options)
            cmd_list.extend(other_parts)
        else:
            cmd_list.append(other_options)
        if "--fastq" not in str(other_options):
            cmd_list.append("--fastq")
    else:
        cmd_list.append("--fastq")

    # Log the command
    logging.info("[NanoSim] Running command: %s", " ".join(cmd_list))

    try:
        # Run the command using the project's run_command utility
        # SECURITY: Never use shell=True
        run_command(
            cmd_list,
            timeout=timeout,
            stderr_log_level=logging.INFO,
            stderr_prefix="[NanoSim] ",
        )

        # NanoSim output is already logged by run_command
        # Standard output file pattern for NanoSim
        output_fastq = f"{output_prefix}_aligned_reads.fastq"

        if not Path(output_fastq).exists():
            # Check for alternative naming pattern
            alt_output_fastq = f"{output_prefix}.fastq"
            if Path(alt_output_fastq).exists():
                output_fastq = alt_output_fastq
            else:
                error_msg = f"Expected output file not found: {output_fastq}"
                logging.error("[NanoSim] %s", error_msg)
                raise RuntimeError(error_msg)

        logging.info("[NanoSim] Simulation completed successfully")
        logging.info("[NanoSim] Generated FASTQ: %s", output_fastq)

        return output_fastq

    except Exception as e:
        logging.error("[NanoSim] Unexpected error during execution: %s", e)
        msg = str(e)
        raise RuntimeError(f"NanoSim simulation failed: {msg}") from e

    finally:
        # No need to restore working directory since we didn't change it
        pass
