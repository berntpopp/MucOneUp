#!/usr/bin/env python3
"""
NanoSim wrapper for Oxford Nanopore read simulation.

This module provides wrapper functions for NanoSim, a nanopore sequencing read
simulator. It enables integration with the MucOneUp pipeline and handles proper
command construction, execution, and error handling.
"""

import logging
import subprocess
import tempfile
from pathlib import Path

from ...exceptions import ExternalToolError
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


def align_ont_reads_with_minimap2(
    minimap2_cmd: str,
    samtools_cmd: str,
    human_reference: str,
    reads_fastq: str,
    output_bam: str,
    threads: int = 4,
    timeout: int = 1800,
) -> str:
    """
    Align Oxford Nanopore reads to a reference using minimap2 and convert to BAM.

    Args:
        minimap2_cmd: Path to minimap2 command
        samtools_cmd: Path to samtools command
        human_reference: Path to human reference FASTA
        reads_fastq: Path to reads FASTQ
        output_bam: Path for output BAM
        threads: Number of threads to use (default: 4)
        timeout: Timeout in seconds (default: 1800)

    Returns:
        Path to the output BAM file

    Raises:
        RuntimeError: If alignment or conversion fails
    """
    # Ensure output directory exists
    Path(output_bam).resolve().parent.mkdir(parents=True, exist_ok=True)

    # Create an intermediate SAM file
    with tempfile.NamedTemporaryFile(suffix=".sam", delete=False) as temp_sam:
        sam_path = temp_sam.name

    try:
        # Construct alignment command
        # SECURITY: Always use list form with subprocess, never shell=True
        # Use centralized build_tool_command to safely handle multi-word commands (conda/mamba)
        align_cmd_list = build_tool_command(
            minimap2_cmd,
            "-t",
            threads,  # build_tool_command handles conversion
            "-ax",
            "map-ont",
            human_reference,
            reads_fastq,
        )

        logging.info("[minimap2] Running alignment: %s > %s", " ".join(align_cmd_list), sam_path)

        # Run with output redirection to SAM file using subprocess
        try:
            with Path(sam_path).open("w") as sam_file:
                result = subprocess.run(
                    align_cmd_list,
                    stdout=sam_file,
                    stderr=subprocess.PIPE,
                    timeout=timeout,
                    check=True,
                    text=True,
                )
            if result.stderr:
                logging.info("[minimap2] %s", result.stderr)
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr if e.stderr else "Unknown error"
            raise ExternalToolError(
                tool="minimap2",
                exit_code=e.returncode,
                stderr=error_msg,
                cmd=" ".join(align_cmd_list),
            ) from e
        except subprocess.TimeoutExpired as e:
            raise ExternalToolError(
                tool="minimap2",
                exit_code=-1,
                stderr=f"Timed out after {timeout}s",
                cmd=" ".join(align_cmd_list),
            ) from e

        # Convert SAM to BAM
        # SECURITY: Always use list form, never shell=True
        # Use centralized build_tool_command to safely handle multi-word commands (conda/mamba)
        sam_to_bam_cmd = build_tool_command(
            samtools_cmd,
            "view",
            "-@",
            threads,  # build_tool_command handles conversion
            "-b",
            "-h",
            "-F",
            "4",  # Filter unmapped reads
            "-o",
            output_bam + ".unsorted",
            sam_path,
        )

        logging.info("[samtools] Converting SAM to BAM: %s", " ".join(sam_to_bam_cmd))

        # Run SAM to BAM conversion
        run_command(
            sam_to_bam_cmd,
            timeout=timeout,
            stderr_log_level=logging.INFO,
            stderr_prefix="[samtools] ",
        )

        # Sort BAM
        # SECURITY: Always use list form, never shell=True
        # Use centralized build_tool_command to safely handle multi-word commands (conda/mamba)
        sort_cmd_list = build_tool_command(
            samtools_cmd,
            "sort",
            "-@",
            threads,  # build_tool_command handles conversion
            "-o",
            output_bam,
            output_bam + ".unsorted",
        )

        logging.info("[samtools] Sorting BAM: %s", " ".join(sort_cmd_list))

        # Run BAM sorting
        run_command(
            sort_cmd_list,
            timeout=timeout,
            stderr_log_level=logging.INFO,
            stderr_prefix="[samtools] ",
        )

        # Index BAM
        # SECURITY: Always use list form, never shell=True
        # Use centralized build_tool_command to safely handle multi-word commands (conda/mamba)
        index_cmd_list = build_tool_command(samtools_cmd, "index", output_bam)

        logging.info("[samtools] Indexing BAM: %s", " ".join(index_cmd_list))

        # Run BAM indexing
        run_command(
            index_cmd_list,
            timeout=timeout,
            stderr_log_level=logging.INFO,
            stderr_prefix="[samtools] ",
        )

        logging.info("[minimap2] ONT read alignment completed successfully")
        logging.info("[samtools] Generated BAM: %s", output_bam)

        return output_bam

    except Exception as e:
        logging.error("[minimap2] Unexpected error during ONT alignment: %s", e)
        error_msg = f"ONT read alignment failed: {str(e)[:60]}..."
        raise RuntimeError(error_msg) from e

    finally:
        # Clean up temporary files
        sam_path_obj = Path(sam_path)
        if sam_path_obj.exists():
            sam_path_obj.unlink()
        unsorted_path = Path(output_bam + ".unsorted")
        if unsorted_path.exists():
            unsorted_path.unlink()
