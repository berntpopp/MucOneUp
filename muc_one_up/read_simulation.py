#!/usr/bin/env python3
"""
read_simulation.py

Command-line entry point for the MUC1 read simulation pipeline.

This module provides a standardized interface to the read simulation pipeline,
which creates realistic sequencing reads from a simulated MUC1 haplotype FASTA.
It supports Illumina short reads (via a hybrid pipeline combining Wessim2-style
fragment simulation with ReSeq error modeling and BWA alignment), Oxford Nanopore
long reads (via NanoSim), and PacBio HiFi reads (via pbsim3/CCS).

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

from __future__ import annotations

import logging
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from muc_one_up.read_simulator.output_config import OutputConfig


def _get_simulator(simulator_type: str) -> Callable[..., str]:
    """Return the pipeline function for the given simulator type.

    Imports only the selected backend to avoid loading heavy dependencies
    (e.g., Biopython for amplicon) when they are not needed.
    """
    if simulator_type == "illumina":
        from muc_one_up.read_simulator.pipeline import (
            simulate_reads_pipeline as simulate_illumina_reads,
        )

        return lambda config, input_fa, _, **kw: simulate_illumina_reads(config, input_fa, **kw)
    elif simulator_type == "ont":
        from muc_one_up.read_simulator.ont_pipeline import simulate_ont_reads_pipeline

        return lambda config, input_fa, human_reference, **kw: simulate_ont_reads_pipeline(
            config, input_fa, human_reference=human_reference, **kw
        )
    elif simulator_type == "pacbio":
        from muc_one_up.read_simulator.pacbio_pipeline import simulate_pacbio_hifi_reads

        return lambda config, input_fa, human_reference, **kw: simulate_pacbio_hifi_reads(
            config, input_fa, human_reference=human_reference, **kw
        )
    elif simulator_type == "amplicon":
        from muc_one_up.read_simulator.amplicon_pipeline import (
            simulate_amplicon_reads_pipeline,
        )

        return lambda config, input_fa, human_reference, **kw: simulate_amplicon_reads_pipeline(
            config, input_fa, human_reference=human_reference, **kw
        )
    else:
        valid = "amplicon, illumina, ont, pacbio"
        raise ValueError(f"Unknown simulator: '{simulator_type}'. Valid options: {valid}. ")


def simulate_reads(
    config: dict[str, Any],
    input_fa: str,
    source_tracker: Any | None = None,
    output_config: OutputConfig | None = None,
) -> str:
    """
    Run the complete read simulation pipeline using Strategy Pattern.

    This function serves as the main API entry point for the read simulation pipeline.
    It selects the appropriate pipeline based on the specified simulator in the config,
    and delegates to the modular implementation in the read_simulator package.

    The Strategy Pattern enables extensible simulator support without modifying this
    function - new simulators can be added by registering them in _get_simulator().

    Args:
        config: Dictionary containing configuration sections:
               - 'read_simulation': Contains 'simulator' ('illumina', 'ont', or 'pacbio')
                 and other parameters
               - 'tools': Contains paths to required external tools
               - 'nanosim_params': Required for ONT simulation
               - 'pacbio_params': Required for PacBio simulation
               See the respective pipeline documentation for detailed parameter information.
        input_fa: Input simulated FASTA file (e.g., muc1_simulated.fa).
        source_tracker: Optional read source tracker for provenance.
        output_config: Optional OutputConfig controlling output directory and
               base name. When provided, overrides default output placement
               (which derives paths from the input file).

    Returns:
        Path to the final output BAM file (or FASTQ if alignment skipped).

    Raises:
        ValueError: If simulator is unknown or not registered in _get_simulator().

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

    # Get the pipeline function for the selected simulator (lazy import)
    simulator_func = _get_simulator(simulator)

    # Log which simulator is being used
    simulator_names = {
        "illumina": "Illumina read simulation pipeline (Wessim2-style fragments + ReSeq error modeling)",
        "ont": "Oxford Nanopore (ONT) read simulation pipeline with NanoSim",
        "pacbio": "PacBio HiFi read simulation pipeline with pbsim3/CCS",
        "amplicon": "PacBio amplicon read simulation pipeline with pbsim3/CCS (template mode)",
    }
    logging.info(f"Using {simulator_names.get(simulator, f'{simulator} simulator')}")

    # Warn if human reference is missing (for alignment-based pipelines)
    if simulator in ["ont", "pacbio", "amplicon"] and not human_reference:
        logging.warning(
            f"No human_reference specified in config for {simulator.upper()} alignment. "
            "Alignment step will be skipped (FASTQ output only)."
        )

    # Dispatch to appropriate pipeline using Strategy Pattern
    return simulator_func(
        config,
        input_fa,
        human_reference,
        source_tracker=source_tracker,
        output_config=output_config,
    )


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
