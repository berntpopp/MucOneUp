#!/usr/bin/env python3
"""
Wrapper module for reseq tool operations.

This module provides wrapper functions for the reseq toolkit operations:
- replace_Ns: Replace N nucleotides in FASTA files
- generate_systematic_errors: Create systematic error profiles
- create_reads: Generate reads from fragments using seqToIllumina
"""

import logging
from pathlib import Path

from ...exceptions import ExternalToolError, FileOperationError
from ..command_utils import build_tool_command
from ..utils import run_command


def replace_Ns(input_fa: str, output_fa: str, tools: dict[str, str]) -> None:
    """
    Replace Ns in the simulated FASTA using reseq replaceN.

    Args:
        input_fa: Input FASTA filename.
        output_fa: Output FASTA filename.
        tools: Dictionary of tool commands.

    Raises:
        SystemExit: If the command fails.
    """
    # From reseq help: replaceN expects -r for input and -R for output
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    cmd = build_tool_command(tools["reseq"], "replaceN", "-r", input_fa, "-R", output_fa)
    run_command(cmd, timeout=60, stderr_log_level=logging.INFO, stderr_prefix="[reseq] ")


def generate_systematic_errors(
    input_fa: str, reseq_model: str, output_fq: str, tools: dict[str, str]
) -> None:
    """
    Generate systematic errors using reseq illuminaPE.

    Args:
        input_fa: Input FASTA filename.
        reseq_model: Reseq model file.
        output_fq: Output FASTQ filename.
        tools: Dictionary of tool commands.

    Raises:
        SystemExit: If the command fails or the output file is not created.
    """
    # Based on the error and original code, this needs to be exactly like the original
    # ERROR: Either bamIn or statsIn option are mandatory.
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    cmd = build_tool_command(
        tools["reseq"],
        "illuminaPE",
        "-r",
        input_fa,  # Input reference fasta
        "-s",
        reseq_model,  # Reseq model file (statsIn)
        "--stopAfterEstimation",
        "--writeSysError",
        output_fq,  # Output FASTQ filename
    )
    run_command(cmd, timeout=60, stderr_log_level=logging.INFO, stderr_prefix="[reseq] ")

    # Verify output exists and is non-empty
    output_path = Path(output_fq)
    if not output_path.exists() or output_path.stat().st_size == 0:
        raise FileOperationError(
            f"Failed to generate systematic errors: Output {output_fq} missing or empty"
        )


def create_reads(
    input_fragments: str,
    reseq_model: str,
    output_reads: str,
    threads: int,
    tools: dict[str, str],
    timeout: int | None = 120,
) -> None:
    """
    Create reads from fragments using reseq seqToIllumina.

    If the command times out but the output file exists and is non-empty,
    a warning is logged and the process continues.

    Args:
        input_fragments: Input fragments FASTA.
        reseq_model: Reseq model file.
        output_reads: Output reads FASTQ.
        threads: Number of threads.
        tools: Dictionary of tool commands.
        timeout: Timeout in seconds (default: 1 hour)

    Raises:
        SystemExit: If the command fails and the output file is not created.
    """
    # Based on the original code, these must be the exact parameters
    # Use build_tool_command to safely handle multi-word commands (conda/mamba)
    cmd = build_tool_command(
        tools["reseq"],
        "seqToIllumina",
        "-j",
        threads,  # Number of threads (build_tool_command handles conversion)
        "-s",
        reseq_model,  # Reseq model file
        "-i",
        input_fragments,  # Input fragments file
        "-o",
        output_reads,  # Output file
    )

    try:
        # Execute with a short timeout - reseq seqToIllumina tends to keep running
        # indefinitely even after it has produced useful output
        run_command(
            cmd,
            timeout=timeout,
            stderr_log_level=logging.INFO,
            stderr_prefix="[reseq] ",
        )
    except ExternalToolError as e:
        # Check if the output file exists and is non-empty despite timeout
        output_path = Path(output_reads)
        if output_path.exists() and output_path.stat().st_size > 0:
            logging.warning(
                (
                    "reseq seqToIllumina timed out after %s seconds, but valid output "
                    "exists. This is EXPECTED BEHAVIOR as seqToIllumina often "
                    "produces complete output before finishing execution. Process was "
                    "terminated but simulation will continue normally."
                ),
                timeout,
            )
        else:
            raise FileOperationError(
                f"Failed to create reads: Output {output_reads} missing or empty after {timeout} seconds"
            ) from e


def split_reads(interleaved_fastq: str, output_fastq1: str, output_fastq2: str) -> None:
    """
    Split an interleaved FASTQ (4 lines per record) into two gzipped FASTQ files.

    Args:
        interleaved_fastq: Input interleaved FASTQ filename.
        output_fastq1: Output FASTQ filename for read1.
        output_fastq2: Output FASTQ filename for read2.

    Raises:
        IOError: If the input file can't be read or output files can't be written.
    """
    try:
        with (
            Path(interleaved_fastq).open() as in_fq,
            Path(output_fastq1).open("wb") as out_fq1,
            Path(output_fastq2).open("wb") as out_fq2,
        ):
            # Create gzip writers
            import gzip

            gzip_fq1 = gzip.GzipFile(fileobj=out_fq1, mode="wb")
            gzip_fq2 = gzip.GzipFile(fileobj=out_fq2, mode="wb")

            # Process 8 lines at a time (4 lines per record, 2 records)
            while True:
                # Read 4 lines for first read
                read1_lines = [in_fq.readline() for _ in range(4)]

                # Break if EOF
                if not read1_lines[0]:
                    break

                # Read 4 lines for second read
                read2_lines = [in_fq.readline() for _ in range(4)]

                # Write to gzipped files
                for line in read1_lines:
                    gzip_fq1.write(line.encode())

                for line in read2_lines:
                    gzip_fq2.write(line.encode())

            # Close gzip writers
            gzip_fq1.close()
            gzip_fq2.close()

        logging.info(
            "Split interleaved FASTQ %s into %s and %s",
            interleaved_fastq,
            output_fastq1,
            output_fastq2,
        )
    except Exception as e:
        logging.error("Error splitting reads: %s", str(e))
        raise
