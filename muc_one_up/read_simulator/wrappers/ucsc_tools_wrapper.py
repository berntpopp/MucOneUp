#!/usr/bin/env python3
"""
Wrapper module for UCSC tools operations.

This module provides wrapper functions for UCSC tools operations:
- fa_to_twobit: Convert FASTA to 2bit format using faToTwoBit
- run_pblat: Align sequences using pblat (parallel BLAT)
"""

import logging
import sys
from pathlib import Path

from ..utils import run_command


def fa_to_twobit(
    input_fa: str, output_2bit: str, tools: dict[str, str], timeout: int | None = 60
) -> None:
    """
    Convert a FASTA file to 2bit format using faToTwoBit.

    Args:
        input_fa: Input FASTA filename.
        output_2bit: Output 2bit filename.
        tools: Dictionary of tool commands.

    Raises:
        SystemExit: If the conversion fails.
    """
    cmd = [tools["faToTwoBit"], input_fa, output_2bit]
    run_command(cmd, timeout=60, stderr_prefix="[ucsc] ", stderr_log_level=logging.INFO)

    # Check output exists and is non-empty
    output_path = Path(output_2bit)
    if not output_path.exists() or output_path.stat().st_size == 0:
        logging.error(
            "Failed to convert FASTA to 2bit: Output file %s missing or empty.",
            output_2bit,
        )
        sys.exit(1)


def run_pblat(
    twobit_file: str,
    subset_reference: str,
    output_psl: str,
    tools: dict[str, str],
    threads: int = 24,
    min_score: int = 95,
    min_identity: int = 95,
    timeout: int | None = 120,
) -> None:
    """
    Run pblat to align a 2bit file to a subset reference.

    Args:
        twobit_file: Input 2bit filename.
        subset_reference: Subset reference FASTA.
        output_psl: Output PSL filename.
        tools: Dictionary of tool commands.
        threads: Number of threads.
        min_score: Minimal score.
        min_identity: Minimal identity.
        timeout: Timeout in seconds (default: 1 hour)

    Raises:
        SystemExit: If the pblat alignment fails.
    """
    # Check input files exist
    for file in [twobit_file, subset_reference]:
        file_path = Path(file)
        if not file_path.exists() or file_path.stat().st_size == 0:
            logging.error(f"Input file missing or empty: {file}")
            sys.exit(1)

    cmd = [
        tools["pblat"],
        f"-threads={threads}",
        f"-minScore={min_score}",
        f"-minIdentity={min_identity}",
        twobit_file,
        subset_reference,
        output_psl,
    ]

    try:
        run_command(
            cmd,
            timeout=timeout,
            stderr_prefix="[pblat] ",
            stderr_log_level=logging.INFO,
        )
    except SystemExit:
        # Check if output was produced despite error
        output_path = Path(output_psl)
        if output_path.exists() and output_path.stat().st_size > 0:
            logging.warning(
                "pblat command failed, but output file exists. Continuing with caution."
            )
        else:
            logging.error("pblat alignment failed and no output was produced.")
            sys.exit(1)

    # Verify output exists and is non-empty
    output_path = Path(output_psl)
    if not output_path.exists() or output_path.stat().st_size == 0:
        logging.error(f"pblat output file missing or empty: {output_psl}")
        sys.exit(1)

    logging.info(f"Successfully created PSL file: {output_psl}")
