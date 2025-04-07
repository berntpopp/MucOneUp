#!/usr/bin/env python3
"""
Wrapper module for UCSC tools operations.

This module provides wrapper functions for UCSC tools operations:
- fa_to_twobit: Convert FASTA to 2bit format using faToTwoBit
- run_pblat: Align sequences using pblat (parallel BLAT)
"""

import logging
import os
import sys
from typing import Dict, Optional

from ..utils import run_command


def fa_to_twobit(
    input_fa: str, output_2bit: str, tools: Dict[str, str], timeout: Optional[int] = 60
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
    run_command(cmd, timeout=60)

    # Check output exists and is non-empty
    if not os.path.exists(output_2bit) or os.path.getsize(output_2bit) == 0:
        logging.error(
            "Failed to convert FASTA to 2bit: Output file %s missing or empty.",
            output_2bit,
        )
        sys.exit(1)


def run_pblat(
    twobit_file: str,
    subset_reference: str,
    output_psl: str,
    tools: Dict[str, str],
    threads: int = 24,
    minScore: int = 95,
    minIdentity: int = 95,
    timeout: Optional[int] = 120,
) -> None:
    """
    Run pblat to align a 2bit file to a subset reference.

    Args:
        twobit_file: Input 2bit filename.
        subset_reference: Subset reference FASTA.
        output_psl: Output PSL filename.
        tools: Dictionary of tool commands.
        threads: Number of threads.
        minScore: Minimal score.
        minIdentity: Minimal identity.
        timeout: Timeout in seconds (default: 1 hour)

    Raises:
        SystemExit: If the pblat alignment fails.
    """
    # Check input files exist
    for file in [twobit_file, subset_reference]:
        if not os.path.exists(file) or os.path.getsize(file) == 0:
            logging.error(f"Input file missing or empty: {file}")
            sys.exit(1)

    cmd = [
        tools["pblat"],
        "-threads={}".format(threads),
        "-minScore={}".format(minScore),
        "-minIdentity={}".format(minIdentity),
        twobit_file,
        subset_reference,
        output_psl,
    ]

    try:
        run_command(cmd, timeout=timeout)
    except SystemExit as e:
        # Check if output was produced despite error
        if os.path.exists(output_psl) and os.path.getsize(output_psl) > 0:
            logging.warning(
                "pblat command failed, but output file exists. Continuing with caution."
            )
        else:
            logging.error("pblat alignment failed and no output was produced.")
            sys.exit(1)

    # Verify output exists and is non-empty
    if not os.path.exists(output_psl) or os.path.getsize(output_psl) == 0:
        logging.error(f"pblat output file missing or empty: {output_psl}")
        sys.exit(1)

    logging.info(f"Successfully created PSL file: {output_psl}")
