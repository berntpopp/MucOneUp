#!/usr/bin/env python3
"""
read_simulation.py

Command-line entry point for the MUC1 read simulation pipeline.

This module provides a standardized interface to the read simulation pipeline,
which creates realistic sequencing reads from a simulated MUC1 haplotype FASTA.
It supports Illumina short reads (via reseq/WeSSim), Oxford Nanopore long reads
(via NanoSim), and PacBio HiFi reads (via pbsim3/CCS).

It serves as a compatibility layer for the refactored modular implementations
in the muc_one_up.read_simulator package.

Key features:
- Simulator selection (Illumina, Oxford Nanopore, or PacBio HiFi)
- Strategy Pattern for extensible simulator support
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
- PacBio HiFi pipeline:
  - pbsim3: For multi-pass CLR read simulation
  - CCS: For HiFi consensus generation
  - minimap2: For HiFi read alignment (map-hifi preset)
  - samtools: For BAM/FASTQ conversion

Output files:
- {input_basename}.bam: Aligned and indexed BAM file (Illumina)
- {input_basename}_R1.fastq.gz: Forward reads (gzipped FASTQ) (Illumina)
- {input_basename}_R2.fastq.gz: Reverse reads (gzipped FASTQ) (Illumina)
- {input_basename}_ont.bam: Aligned and indexed BAM file (ONT)
- {input_basename}_ont_aligned_reads.fastq: ONT reads (FASTQ)
- {input_basename}_aligned.bam: Aligned and indexed BAM file (PacBio)
- {input_basename}_hifi.fastq: PacBio HiFi reads (FASTQ)

Usage:
    python read_simulation.py <config.json> <input_fasta>
where <input_fasta> is typically the output from muconeup
(e.g., muc1_simulated.fa).
"""

import logging
from collections.abc import Callable
from pathlib import Path
from typing import Any

from muc_one_up.read_simulator.ont_pipeline import simulate_ont_reads_pipeline
from muc_one_up.read_simulator.pacbio_pipeline import simulate_pacbio_hifi_reads

# Import the refactored pipelines
from muc_one_up.read_simulator.pipeline import (
    simulate_reads_pipeline as simulate_illumina_reads,
)

# Configure logging
logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s")


# =============================================================================
# Strategy Pattern: Simulator Selection
# =============================================================================
#
# Maps simulator names to pipeline functions for extensible simulator support.
# This design pattern (Strategy Pattern) follows the Open/Closed Principle:
# - Open for extension: New simulators can be added without modifying existing code
# - Closed for modification: Adding a simulator doesn't change the router logic
#
# To add a new simulator:
# 1. Create a new pipeline module (e.g., nanopore_promethion_pipeline.py)
# 2. Import the pipeline function
# 3. Add an entry to SIMULATOR_MAP
# 4. No changes to simulate_reads() required!
#
# Example:
#     SIMULATOR_MAP["promethion"] = simulate_promethion_reads_pipeline
#
SIMULATOR_MAP: dict[str, Callable[[dict[str, Any], str, str | None], str]] = {
    "illumina": lambda config, input_fa, _: simulate_illumina_reads(config, input_fa),
    "ont": simulate_ont_reads_pipeline,
    "pacbio": simulate_pacbio_hifi_reads,
}


def simulate_reads(config: dict[str, Any], input_fa: str) -> str:
    """
    Run the complete read simulation pipeline using Strategy Pattern.

    This function serves as the main API entry point for the read simulation pipeline.
    It selects the appropriate pipeline based on the specified simulator in the config,
    and delegates to the modular implementation in the read_simulator package.

    The Strategy Pattern enables extensible simulator support without modifying this
    function - new simulators can be added by registering them in SIMULATOR_MAP.

    Args:
        config: Dictionary containing configuration sections:
               - 'read_simulation': Contains 'simulator' ('illumina', 'ont', or 'pacbio')
                 and other parameters
               - 'tools': Contains paths to required external tools
               - 'nanosim_params': Required for ONT simulation
               - 'pacbio_params': Required for PacBio simulation
               See the respective pipeline documentation for detailed parameter information.
        input_fa: Input simulated FASTA file (e.g., muc1_simulated.fa).

    Returns:
        Path to the final output BAM file (or FASTQ if alignment skipped).

    Raises:
        ValueError: If simulator is unknown or not registered in SIMULATOR_MAP.

    Example:
        Illumina simulation::

            config = {
                "read_simulation": {"simulator": "illumina", ...},
                "tools": {...}
            }
            bam = simulate_reads(config, "sample.fa")

        ONT simulation::

            config = {
                "read_simulation": {"simulator": "ont", ...},
                "nanosim_params": {...}
            }
            bam = simulate_reads(config, "sample.fa")

        PacBio HiFi simulation::

            config = {
                "read_simulation": {"simulator": "pacbio", ...},
                "pacbio_params": {...}
            }
            bam = simulate_reads(config, "sample.fa")

    Notes:
        - Simulator names are case-insensitive
        - Default simulator is "illumina" if not specified
        - Human reference is extracted from config for alignment
        - Strategy Pattern makes adding new simulators trivial
    """
    # Determine which simulator to use (case-insensitive, default to illumina)
    simulator = config.get("read_simulation", {}).get("simulator", "illumina").lower()

    # Extract human reference from config (needed for ONT and PacBio alignment)
    human_reference = config.get("read_simulation", {}).get("human_reference")

    # Validate simulator is registered
    if simulator not in SIMULATOR_MAP:
        valid_simulators = ", ".join(sorted(SIMULATOR_MAP.keys()))
        raise ValueError(
            f"Unknown simulator: '{simulator}'. "
            f"Valid options: {valid_simulators}. "
            f"To add a new simulator, register it in SIMULATOR_MAP."
        )

    # Log which simulator is being used
    simulator_names = {
        "illumina": "Illumina read simulation pipeline with reseq/WeSSim",
        "ont": "Oxford Nanopore (ONT) read simulation pipeline with NanoSim",
        "pacbio": "PacBio HiFi read simulation pipeline with pbsim3/CCS",
    }
    logging.info(f"Using {simulator_names.get(simulator, f'{simulator} simulator')}")

    # Warn if human reference is missing (for alignment-based pipelines)
    if simulator in ["ont", "pacbio"] and not human_reference:
        logging.warning(
            f"No human_reference specified in config for {simulator.upper()} alignment. "
            "Alignment step will be skipped (FASTQ output only)."
        )

    # Dispatch to appropriate pipeline using Strategy Pattern
    simulator_func = SIMULATOR_MAP[simulator]
    return simulator_func(config, input_fa, human_reference)


if __name__ == "__main__":  # OK: top-level entry point
    import json
    import sys

    if len(sys.argv) != 3:
        print("Usage: python read_simulation.py <config.json> <input_fasta>")
        sys.exit(1)
    config_file = sys.argv[1]
    input_fa = sys.argv[2]
    with Path(config_file).open() as fh:
        config = json.load(fh)
    simulate_reads(config, input_fa)
