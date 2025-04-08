#!/usr/bin/env python3
"""
NanoSim wrapper for Oxford Nanopore read simulation.

This module provides wrapper functions for NanoSim, a nanopore sequencing read
simulator. It enables integration with the MucOneUp pipeline and handles proper
command construction, execution, and error handling.
"""

import logging
import os
import subprocess
import tempfile
from typing import Optional

from ..utils import run_command


def run_nanosim_simulation(
    nanosim_cmd: str,
    reference_fasta: str,
    output_prefix: str,
    training_model: str,
    coverage: float,
    threads: int = 4,
    min_read_length: Optional[int] = None,
    max_read_length: Optional[int] = None,
    other_options: str = "",
    timeout: int = 3600,
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

    Returns:
        Path to the generated FASTQ file

    Raises:
        RuntimeError: If NanoSim execution fails
    """
    # Ensure output directory exists
    output_dir = os.path.dirname(os.path.abspath(output_prefix))
    os.makedirs(output_dir, exist_ok=True)

    # Build the command with required parameters
    if isinstance(nanosim_cmd, str) and (" " in nanosim_cmd):
        # Command contains spaces, it's likely a conda/mamba run command
        # Let run_command handle it as a string
        cmd = f"{nanosim_cmd} genome -rg {reference_fasta}"
        cmd += f" -c {training_model} -o {output_prefix}"
        cmd += f" -t {threads} -x {coverage}"
    else:
        # Build as list for direct command
        cmd = [
            nanosim_cmd,
            "genome",
            "-rg",
            reference_fasta,
            "-c",
            training_model,
            "-o",
            output_prefix,
            "-t",
            str(threads),
            "-x",
            str(coverage),
        ]

    # Add optional parameters if provided
    if isinstance(cmd, list):
        if min_read_length:
            cmd.extend(["--min_len", str(min_read_length)])
        if max_read_length:
            cmd.extend(["--max_len", str(max_read_length)])
        if other_options:
            cmd.extend(other_options.split())
        if "--fastq" not in other_options:
            cmd.append("--fastq")
    else:  # String command
        if min_read_length:
            cmd += f" --min_len {min_read_length}"
        if max_read_length:
            cmd += f" --max_len {max_read_length}"
        if other_options:
            cmd += f" {other_options}"
        if "--fastq" not in other_options:
            cmd += " --fastq"

    # Log the command
    logging.info("[NanoSim] Running command: %s", cmd)

    try:
        # Run the command using the project's run_command utility
        run_command(
            cmd,
            shell=isinstance(cmd, str),
            timeout=timeout,
            stderr_log_level=logging.INFO,
            stderr_prefix="[NanoSim] ",
        )

        # NanoSim output is already logged by run_command
        # Standard output file pattern for NanoSim
        output_fastq = f"{output_prefix}_aligned_reads.fastq"

        if not os.path.exists(output_fastq):
            # Check for alternative naming pattern
            alt_output_fastq = f"{output_prefix}.fastq"
            if os.path.exists(alt_output_fastq):
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
        raise RuntimeError(f"NanoSim simulation failed: {msg}")

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
    os.makedirs(os.path.dirname(os.path.abspath(output_bam)), exist_ok=True)

    # Create an intermediate SAM file
    with tempfile.NamedTemporaryFile(suffix=".sam", delete=False) as temp_sam:
        sam_path = temp_sam.name

    try:
        # Construct alignment command
        if isinstance(minimap2_cmd, str) and (" " in minimap2_cmd):
            # Command has spaces, handle as a string command
            align_cmd = f"{minimap2_cmd} -t {threads} -ax map-ont "
            align_cmd += f"{human_reference} {reads_fastq} > {sam_path}"
            logging.info("[minimap2] Running alignment: %s", align_cmd)
            # Run alignment with shell=True since we're using redirection
            run_command(
                align_cmd,
                shell=True,
                timeout=timeout,
                stderr_log_level=logging.INFO,
                stderr_prefix="[minimap2] ",
            )
        else:
            # Use list form for direct command
            align_cmd = [
                minimap2_cmd,
                "-t",
                str(threads),
                "-ax",
                "map-ont",  # Preset for Oxford Nanopore reads
                human_reference,
                reads_fastq,
            ]
            align_cmd_str = " ".join(str(c) for c in align_cmd)
            logging.info("Running minimap2 alignment: %s", align_cmd_str)

            # Run with output redirection to SAM file
            with open(sam_path, "w") as sam_file:
                # Use run_command when possible, but for file redirection
                # we need to use subprocess.Popen directly
                process = subprocess.Popen(
                    align_cmd, stdout=sam_file, stderr=subprocess.PIPE, text=True
                )
                # Get stderr for logging
                _, stderr = process.communicate(timeout=timeout)
                if process.returncode != 0:
                    logging.error("[minimap2] Failed with code %d", process.returncode)
                    logging.error("[minimap2] stderr: %s", stderr)
                    return None

        # Convert SAM to BAM
        if isinstance(samtools_cmd, str) and (" " in samtools_cmd):
            # Execute as string command
            view_cmd = (
                f"{samtools_cmd} view -@ {threads} -b -h -F 4"
                f" -o {output_bam}.unsorted {sam_path}"
            )
            logging.info("[samtools] Converting SAM to BAM: %s", view_cmd)
            run_command(
                view_cmd,
                shell=True,
                timeout=timeout,
                stderr_log_level=logging.INFO,
                stderr_prefix="[samtools] ",
            )
        else:
            # List form for direct command
            sam_to_bam_cmd = [
                samtools_cmd,
                "view",
                "-@",
                str(threads),
                "-b",
                "-h",
                "-F",
                "4",  # Filter unmapped reads
                "-o",
                output_bam + ".unsorted",
                sam_path,
            ]

            cmd_str = " ".join(str(c) for c in sam_to_bam_cmd)
            logging.info("[samtools] Converting SAM to BAM: %s", cmd_str)

            # Run SAM to BAM conversion
            run_command(
                sam_to_bam_cmd,
                timeout=timeout,
                stderr_log_level=logging.INFO,
                stderr_prefix="[samtools] ",
            )

        # Sort BAM
        if isinstance(samtools_cmd, str) and (" " in samtools_cmd):
            # Execute as string command
            sort_cmd = (
                f"{samtools_cmd} sort -@ {threads} "
                f"-o {output_bam} {output_bam}.unsorted"
            )
            logging.info("[samtools] Sorting BAM: %s", sort_cmd)
            run_command(
                sort_cmd,
                shell=True,
                timeout=timeout,
                stderr_log_level=logging.INFO,
                stderr_prefix="[samtools] ",
            )
        else:
            # List form for direct command
            sort_cmd = [
                samtools_cmd,
                "sort",
                "-@",
                str(threads),
                "-o",
                output_bam,
                output_bam + ".unsorted",
            ]

            cmd_str = " ".join(str(c) for c in sort_cmd)
            logging.info("[samtools] Sorting BAM: %s", cmd_str)

            # Run BAM sorting
            run_command(
                sort_cmd,
                timeout=timeout,
                stderr_log_level=logging.INFO,
                stderr_prefix="[samtools] ",
            )

        # Index BAM
        if isinstance(samtools_cmd, str) and (" " in samtools_cmd):
            # Execute as string command
            index_cmd = f"{samtools_cmd} index {output_bam}"
            logging.info("[samtools] Indexing BAM: %s", index_cmd)
            run_command(
                index_cmd,
                shell=True,
                timeout=timeout,
                stderr_log_level=logging.INFO,
                stderr_prefix="[samtools] ",
            )
        else:
            # List form for direct command
            index_cmd = [samtools_cmd, "index", output_bam]

            cmd_str = " ".join(str(c) for c in index_cmd)
            logging.info("[samtools] Indexing BAM: %s", cmd_str)

            # Run BAM indexing
            run_command(
                index_cmd,
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
        raise RuntimeError(error_msg)

    finally:
        # Clean up temporary files
        if os.path.exists(sam_path):
            os.unlink(sam_path)
        if os.path.exists(output_bam + ".unsorted"):
            os.unlink(output_bam + ".unsorted")
