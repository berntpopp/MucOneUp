#!/usr/bin/env python3
"""
read_simulation.py

Command-line entry point for the MUC1 read simulation pipeline.

This module provides a standardized interface to the read simulation pipeline,
which creates realistic sequencing reads from a simulated MUC1 haplotype FASTA.
It supports both Illumina short reads (via reseq/WeSSim) and Oxford Nanopore long
reads (via NanoSim).

It serves as a compatibility layer for the refactored modular implementations
in the muc_one_up.read_simulator package.

Key features:
- Simulator selection (Illumina or Oxford Nanopore)
- Reliable timeout handling for all external tool calls
- Consistent naming convention for input and output files
- Comprehensive simulation pipeline from FASTA sequence to aligned BAM output

Supported external bioinformatics tools:
- Illumina pipeline:
  - reseq: For sequence manipulation and read simulation
  - faToTwoBit/pblat: For sequence format conversion and alignment
  - samtools: For BAM manipulation and coverage analysis
  - bwa: For read alignment to reference genome
- Oxford Nanopore pipeline:
  - NanoSim: For ONT read simulation
  - minimap2: For ONT read alignment
  - samtools: For BAM manipulation

Output files:
- {input_basename}.bam: Aligned and indexed BAM file (Illumina)
- {input_basename}_R1.fastq.gz: Forward reads (gzipped FASTQ) (Illumina)
- {input_basename}_R2.fastq.gz: Reverse reads (gzipped FASTQ) (Illumina)
- {input_basename}_ont.bam: Aligned and indexed BAM file (ONT)
- {input_basename}_ont_aligned_reads.fastq: ONT reads (FASTQ)

Usage:
    python read_simulation.py <config.json> <input_fasta>
where <input_fasta> is typically the output from muconeup
(e.g., muc1_simulated.fa).
"""

import logging
import sys
from typing import Dict, Any, Optional

# Import the refactored pipelines
from muc_one_up.read_simulator.pipeline import (
    simulate_reads_pipeline as simulate_illumina_reads,
)
from muc_one_up.read_simulator.ont_pipeline import simulate_ont_reads_pipeline

# Configure logging
logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)


def simulate_reads(config: Dict[str, Any], input_fa: str) -> str:
    """
    Run the complete read simulation pipeline.

    This function serves as the main API entry point for the read simulation pipeline.
    It selects the appropriate pipeline based on the specified simulator in the config,
    and delegates to the modular implementation in the read_simulator package.

    Args:
        config: Dictionary containing configuration sections:
               - 'read_simulation': Contains 'simulator' ('illumina' or 'ont') and other parameters
               - 'tools': Contains paths to required external tools
               - 'nanosim_params': Required for ONT simulation
               See the respective pipeline documentation for detailed parameter information.
        input_fa: Input simulated FASTA file (e.g., muc1_simulated.fa).

    Returns:
        Path to the final output BAM file.
    """
    # Determine which simulator to use
    simulator = config.get("read_simulation", {}).get("simulator", "illumina")

    if simulator.lower() == "ont":
        logging.info(
            "Using Oxford Nanopore (ONT) read simulation pipeline with NanoSim"
        )
        # Extract human reference from config, just like the Illumina pipeline does
        human_reference = config.get("read_simulation", {}).get("human_reference")
        if not human_reference:
            logging.warning(
                "No human_reference specified in config for ONT alignment. "
                "Will use the simulated reference instead."
            )
        return simulate_ont_reads_pipeline(config, input_fa, human_reference)
    else:
        logging.info("Using Illumina read simulation pipeline with reseq/WeSSim")
        return simulate_illumina_reads(config, input_fa)


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
