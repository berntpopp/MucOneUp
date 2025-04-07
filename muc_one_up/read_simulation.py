#!/usr/bin/env python3
"""
read_simulation.py

Command-line entry point for the MUC1 read simulation pipeline.

This module provides a standardized interface to the read simulation pipeline,
which creates realistic Illumina-like reads from a simulated MUC1 haplotype FASTA.
It serves as a compatibility layer for the refactored modular implementation
in the muc_one_up.read_simulator package.

Key features:
- Reliable timeout handling for all external tool calls (60-300s depending on complexity)
- Consistent naming convention for input and output files
- Comprehensive simulation pipeline from FASTA sequence to aligned BAM output

The pipeline uses several external bioinformatics tools:
- reseq: For sequence manipulation and read simulation
- faToTwoBit/pblat: For sequence format conversion and alignment
- samtools: For BAM manipulation and coverage analysis
- bwa: For read alignment to reference genome

Output files:
- {input_basename}.bam: Aligned and indexed BAM file
- {input_basename}_R1.fastq.gz: Forward reads (gzipped FASTQ)
- {input_basename}_R2.fastq.gz: Reverse reads (gzipped FASTQ)

Usage:
    python read_simulation.py <config.json> <input_fasta>
where <input_fasta> is typically the output from muconeup
(e.g., muc1_simulated.fa).
"""

import logging
import sys
from typing import Dict, Any

# Import the refactored pipeline
from muc_one_up.read_simulator.pipeline import simulate_reads_pipeline

# Configure logging
logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)


def simulate_reads(config: Dict[str, Any], input_fa: str) -> str:
    """
    Run the complete read simulation pipeline.

    This function serves as the main API entry point for the read simulation pipeline.
    It delegates to the modular implementation in the read_simulator package.

    Args:
        config: Dictionary containing "tools" and "read_simulation" sections.
               See the pipeline documentation for detailed parameter information.
               Must include paths to required external tools and simulation parameters.
        input_fa: Input simulated FASTA file (e.g., muc1_simulated.fa).

    Returns:
        Path to the final output BAM file with name {input_basename}.bam.
    """
    return simulate_reads_pipeline(config, input_fa)


if __name__ == "__main__":
    import json

    if len(sys.argv) != 3:
        print("Usage: python read_simulation.py <config.json> <input_fasta>")
        sys.exit(1)
    config_file = sys.argv[1]
    input_fa = sys.argv[2]
    with open(config_file, "r") as fh:
        config = json.load(fh)
    simulate_reads(config, input_fa)
