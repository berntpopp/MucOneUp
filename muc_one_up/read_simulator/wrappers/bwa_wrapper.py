#!/usr/bin/env python3
"""
Wrapper module for BWA operations.

This module provides wrapper functions for BWA operations:
- align_reads: Align paired-end reads using bwa mem
"""

import logging
import subprocess
from pathlib import Path

from ...exceptions import ExternalToolError, FileOperationError


def align_reads(
    read1: str,
    read2: str,
    human_reference: str,
    output_bam: str,
    tools: dict[str, str],
    threads: int = 4,
    timeout: int | None = 300,
) -> None:
    """
    Align paired-end reads with bwa mem, sort with samtools, and index the BAM file.

    Args:
        read1: FASTQ filename for read1.
        read2: FASTQ filename for read2.
        human_reference: Human reference FASTA.
        output_bam: Output BAM filename.
        tools: Dictionary of tool commands.
        threads: Number of threads.

    Raises:
        FileOperationError: If input files not found or output files missing
    """
    # Check input files exist
    for file in [read1, read2, human_reference]:
        if not Path(file).exists():
            raise FileOperationError(f"Input file not found for BWA alignment: {file}")

    # Get output directory and filename base for intermediate files
    output_bam_path = Path(output_bam)
    output_dir = output_bam_path.parent
    output_base = output_bam_path.stem

    # Create intermediate filenames for unsorted BAM
    unsorted_bam = str(output_dir / f"{output_base}_unsorted.bam")

    # Run alignment with BWA mem and pipe to samtools to create BAM
    # SECURITY: Use subprocess.Popen for piping, never shell=True
    try:
        # Start BWA mem process
        bwa_cmd = [tools["bwa"], "mem", "-t", str(threads), human_reference, read1, read2]

        logging.info(f"Running BWA mem: {' '.join(bwa_cmd)}")
        bwa_process = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Start samtools view process (reads from BWA stdout)
        samtools_cmd = [tools["samtools"], "view", "-bS", "-"]

        with Path(unsorted_bam).open("wb") as bam_out:
            samtools_process = subprocess.Popen(
                samtools_cmd, stdin=bwa_process.stdout, stdout=bam_out, stderr=subprocess.PIPE
            )

        # Allow BWA to receive SIGPIPE if samtools exits
        if bwa_process.stdout:
            bwa_process.stdout.close()

        # Wait for both processes to complete
        samtools_stderr = samtools_process.communicate(timeout=timeout or 300)[1]
        bwa_stderr = bwa_process.communicate(timeout=5)[1]  # BWA should already be done

        # Check return codes
        if bwa_process.returncode != 0:
            error_msg = (
                bwa_stderr.decode("utf-8", errors="ignore") if bwa_stderr else "Unknown error"
            )
            raise ExternalToolError(
                tool="bwa mem",
                exit_code=bwa_process.returncode,
                stderr=error_msg,
                cmd=" ".join(bwa_cmd),
            )

        if samtools_process.returncode != 0:
            error_msg = (
                samtools_stderr.decode("utf-8", errors="ignore")
                if samtools_stderr
                else "Unknown error"
            )
            raise ExternalToolError(
                tool="samtools view",
                exit_code=samtools_process.returncode,
                stderr=error_msg,
                cmd=" ".join(samtools_cmd),
            )

        logging.info("[bwa+samtools] Alignment and BAM conversion completed")

    except subprocess.TimeoutExpired as e:
        # Kill processes if timeout
        if bwa_process:
            bwa_process.kill()
        if samtools_process:
            samtools_process.kill()
        raise ExternalToolError(
            tool="bwa/samtools",
            exit_code=-1,
            stderr=f"Pipeline timed out after {timeout}s",
        ) from e
    except Exception as e:
        raise ExternalToolError(
            tool="bwa/samtools",
            exit_code=-1,
            stderr=str(e),
        ) from e

    # Sort the BAM file using subprocess.run (secure - no shell)
    sort_cmd = [
        tools["samtools"],
        "sort",
        "-@",
        str(threads),
        "-o",
        output_bam,
        unsorted_bam,
    ]

    try:
        logging.info(f"[samtools] Sorting BAM: {' '.join(sort_cmd)}")
        result = subprocess.run(sort_cmd, capture_output=True, timeout=60, check=True)
        if result.stderr:
            logging.info(f"[samtools] {result.stderr.decode('utf-8', errors='ignore')}")
    except subprocess.CalledProcessError as e:
        raise ExternalToolError(
            tool="samtools sort",
            exit_code=e.returncode,
            stderr=e.stderr.decode("utf-8", errors="ignore") if e.stderr else "",
            cmd=" ".join(sort_cmd),
        ) from e
    except subprocess.TimeoutExpired as e:
        raise ExternalToolError(
            tool="samtools sort",
            exit_code=-1,
            stderr="Timed out after 60s",
            cmd=" ".join(sort_cmd),
        ) from e

    # Index the BAM file using subprocess.run (secure - no shell)
    index_cmd = [tools["samtools"], "index", output_bam]

    try:
        logging.info(f"[samtools] Indexing BAM: {' '.join(index_cmd)}")
        result = subprocess.run(index_cmd, capture_output=True, timeout=60, check=True)
        if result.stderr:
            logging.info(f"[samtools] {result.stderr.decode('utf-8', errors='ignore')}")
    except subprocess.CalledProcessError as e:
        raise ExternalToolError(
            tool="samtools index",
            exit_code=e.returncode,
            stderr=e.stderr.decode("utf-8", errors="ignore") if e.stderr else "",
            cmd=" ".join(index_cmd),
        ) from e
    except subprocess.TimeoutExpired as e:
        raise ExternalToolError(
            tool="samtools index",
            exit_code=-1,
            stderr="Timed out after 60s",
            cmd=" ".join(index_cmd),
        ) from e

    # Check output files exist
    for file in [output_bam, f"{output_bam}.bai"]:
        file_path = Path(file)
        if not file_path.exists() or file_path.stat().st_size == 0:
            raise FileOperationError(f"BWA alignment output file missing or empty: {file}")

    # Clean up intermediate files
    unsorted_bam_path = Path(unsorted_bam)
    if unsorted_bam_path.exists():
        try:
            unsorted_bam_path.unlink()
            logging.info(f"Removed intermediate file: {unsorted_bam}")
        except Exception as e:
            logging.warning(f"Could not remove intermediate file {unsorted_bam}: {e!s}")

    logging.info(f"Successfully aligned reads to {output_bam}")
