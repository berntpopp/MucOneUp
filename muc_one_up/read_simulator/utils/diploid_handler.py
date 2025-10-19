#!/usr/bin/env python3
"""
Diploid split-simulation handler for NanoSim.

This module implements split-simulation to eliminate length-proportional
sampling bias in diploid references. Instead of simulating reads from a
combined diploid FASTA, it:
1. Extracts each haplotype to separate files
2. Simulates each haplotype independently at half coverage
3. Merges the resulting FASTQ files

This reduces allelic imbalance from ~3.27:1 to ~1.40:1 for divergent haplotypes.
"""

import logging
import tempfile
from collections.abc import Callable
from pathlib import Path
from typing import NamedTuple

from .fastq_utils import merge_fastq_files
from .reference_utils import extract_haplotypes, get_reference_info, is_diploid_reference


class DiploidSimulationResult(NamedTuple):
    """Result from diploid split-simulation.

    Attributes:
        merged_fastq: Path to merged FASTQ file containing reads from both haplotypes.
        hap1_fastq: Path to haplotype 1 FASTQ file (before merging).
        hap2_fastq: Path to haplotype 2 FASTQ file (before merging).
        hap1_reference: Path to haplotype 1 reference FASTA.
        hap2_reference: Path to haplotype 2 reference FASTA.
        reads_hap1: Number of reads generated for haplotype 1.
        reads_hap2: Number of reads generated for haplotype 2.
    """

    merged_fastq: str
    hap1_fastq: str
    hap2_fastq: str
    hap1_reference: str
    hap2_reference: str
    reads_hap1: int
    reads_hap2: int


def calculate_corrected_coverage(
    desired_coverage: float,
    correction_factor: float,
    is_split_simulation: bool = True,
) -> float:
    """
    Calculate NanoSim coverage parameter accounting for training model efficiency.

    NanoSim's `-x` parameter doesn't directly produce the requested coverage due to
    empirical kernel density estimation in the training model. The correction factor
    represents the ratio: actual_coverage / requested_coverage.

    For diploid split-simulation, we target half the desired coverage per haplotype
    (since reads will be merged), then apply the correction factor.

    Args:
        desired_coverage: Target coverage after merging (diploid) or total (haploid).
        correction_factor: Training model efficiency (e.g., 0.325 for human_giab_hg002).
        is_split_simulation: If True, calculate for per-haplotype simulation (half coverage).

    Returns:
        Coverage value to pass to NanoSim's `-x` parameter.

    Example:
        >>> # For 200x diploid coverage with 0.325 factor:
        >>> corrected = calculate_corrected_coverage(200, 0.325, is_split_simulation=True)
        >>> print(f"{corrected:.1f}x")
        307.7x
        >>> # This means: (200 / 2) / 0.325 = 307.7x per haplotype
        >>> # After merging: ~100x per haplotype = ~200x total
    """
    if correction_factor <= 0 or correction_factor > 1:
        raise ValueError(f"Correction factor must be between 0 and 1, got {correction_factor}")

    if desired_coverage <= 0:
        raise ValueError(f"Desired coverage must be positive, got {desired_coverage}")

    # For diploid: simulate each haplotype at half coverage, then merge
    # For haploid or standard simulation, use full coverage
    target_per_haplotype = desired_coverage / 2.0 if is_split_simulation else desired_coverage

    # Apply correction factor: actual = requested * factor
    # So: requested = actual / factor
    corrected = target_per_haplotype / correction_factor

    logging.debug(
        "Coverage calculation: desired=%.1fx, factor=%.3f, split=%s → corrected=%.1fx",
        desired_coverage,
        correction_factor,
        is_split_simulation,
        corrected,
    )

    return corrected


def prepare_diploid_simulation(
    diploid_fasta: str | Path,
    output_dir: str | Path,
    base_name: str | None = None,
) -> tuple[str, str, str]:
    """
    Prepare diploid reference for split-simulation by extracting haplotypes.

    Args:
        diploid_fasta: Path to diploid reference (must contain exactly 2 sequences).
        output_dir: Directory for temporary haplotype files.
        base_name: Base name for output files (default: use input filename).

    Returns:
        Tuple of (hap1_fasta, hap2_fasta, temp_dir) paths.

    Raises:
        ValidationError: If input is not a valid diploid reference.

    Example:
        >>> hap1, hap2, tmpdir = prepare_diploid_simulation("diploid.fa", "output")
        >>> print(hap1, hap2)
        output/diploid_hap1.fa output/diploid_hap2.fa
    """
    diploid_fasta = Path(diploid_fasta)
    output_dir = Path(output_dir)

    # Validate diploid reference
    if not is_diploid_reference(diploid_fasta):
        info = get_reference_info(diploid_fasta)
        raise ValueError(
            f"Expected diploid reference (2 sequences), found {info.num_sequences} "
            f"in {diploid_fasta}"
        )

    # Extract haplotypes
    hap1_path, hap2_path = extract_haplotypes(
        diploid_fasta,
        output_dir,
        base_name=base_name,
    )

    logging.info(
        "Prepared diploid simulation: %s → hap1=%s, hap2=%s",
        diploid_fasta.name,
        Path(hap1_path).name,
        Path(hap2_path).name,
    )

    return hap1_path, hap2_path, str(output_dir)


def run_split_simulation(
    diploid_fasta: str | Path,
    simulation_func: Callable,
    simulation_params: dict,
    output_fastq: str | Path,
    correction_factor: float = 0.325,
    keep_intermediate: bool = False,
    seed: int | None = None,
) -> DiploidSimulationResult:
    """
    Run split-simulation workflow for diploid references.

    This function orchestrates the complete split-simulation process:
    1. Extract haplotypes from diploid FASTA
    2. Calculate corrected coverage for each haplotype
    3. Run simulation on each haplotype independently
    4. Merge FASTQ files
    5. Clean up intermediate files (optional)

    Args:
        diploid_fasta: Path to diploid reference FASTA (2 sequences).
        simulation_func: Function to call for each haplotype simulation.
                        Must accept: (reference_fasta, output_prefix, coverage, seed, **params)
                        Must return: path to generated FASTQ file.
        simulation_params: Dictionary of additional parameters to pass to simulation_func.
        output_fastq: Path for merged output FASTQ file.
        correction_factor: NanoSim training model correction factor (default: 0.325).
        keep_intermediate: If True, keep individual haplotype FASTQs (default: False).
        seed: Random seed for reproducibility. If provided, hap2 will use seed+1.

    Returns:
        DiploidSimulationResult with paths and read counts.

    Raises:
        ValueError: If reference is not diploid.
        ValidationError: If simulation fails.

    Example:
        >>> from muc_one_up.read_simulator.wrappers.nanosim_wrapper import run_nanosim_simulation
        >>>
        >>> def sim_func(ref, prefix, cov, seed, **params):
        ...     return run_nanosim_simulation(
        ...         params['nanosim_cmd'], ref, prefix, params['model'],
        ...         cov, seed=seed, threads=params['threads']
        ...     )
        >>>
        >>> result = run_split_simulation(
        ...     "diploid.fa",
        ...     sim_func,
        ...     {'nanosim_cmd': 'simulator.py', 'model': 'path/to/model', 'threads': 8},
        ...     "output.fastq.gz",
        ...     seed=42
        ... )
        >>> print(f"Merged: {result.merged_fastq}")
        Merged: output.fastq.gz
    """
    diploid_fasta = Path(diploid_fasta)
    output_fastq = Path(output_fastq)

    # Get desired coverage from params
    desired_coverage = simulation_params.get("coverage", 30)

    # Create clean params dict without coverage (will be passed separately with correction)
    clean_params = {k: v for k, v in simulation_params.items() if k != "coverage"}

    # Create temporary directory for intermediate files
    with tempfile.TemporaryDirectory(prefix="diploid_sim_") as temp_dir:
        temp_path = Path(temp_dir)

        logging.info("=" * 70)
        logging.info("DIPLOID SPLIT-SIMULATION WORKFLOW")
        logging.info("=" * 70)
        logging.info("Reference: %s", diploid_fasta)
        logging.info(
            "Desired coverage: %.1fx (per haplotype after merge: %.1fx)",
            desired_coverage,
            desired_coverage / 2,
        )
        logging.info("Correction factor: %.3f", correction_factor)
        logging.info("Seed: %s", seed if seed is not None else "random")
        logging.info("=" * 70)

        # Step 1: Extract haplotypes
        logging.info("Step 1/4: Extracting haplotypes...")
        hap1_ref, hap2_ref, _ = prepare_diploid_simulation(
            diploid_fasta,
            temp_path,
            base_name="sim",
        )

        # Step 2: Calculate corrected coverage
        logging.info("Step 2/4: Calculating corrected coverage...")
        corrected_coverage = calculate_corrected_coverage(
            desired_coverage,
            correction_factor,
            is_split_simulation=True,
        )
        logging.info(
            "Corrected coverage per haplotype: %.1fx (NanoSim -x parameter)",
            corrected_coverage,
        )

        # Step 3: Simulate each haplotype
        logging.info("Step 3/4: Simulating haplotypes...")

        # Haplotype 1
        logging.info("  Simulating haplotype 1...")
        hap1_prefix = str(temp_path / "hap1_sim")
        hap1_seed = seed if seed is not None else None
        hap1_fastq = simulation_func(
            reference_fasta=hap1_ref,
            output_prefix=hap1_prefix,
            coverage=corrected_coverage,
            seed=hap1_seed,
            **clean_params,
        )
        logging.info("    Haplotype 1 FASTQ: %s", hap1_fastq)

        # Haplotype 2 (use seed+1 for different but reproducible sequence)
        logging.info("  Simulating haplotype 2...")
        hap2_prefix = str(temp_path / "hap2_sim")
        hap2_seed = (seed + 1) if seed is not None else None
        hap2_fastq = simulation_func(
            reference_fasta=hap2_ref,
            output_prefix=hap2_prefix,
            coverage=corrected_coverage,
            seed=hap2_seed,
            **clean_params,
        )
        logging.info("    Haplotype 2 FASTQ: %s", hap2_fastq)

        # Step 4: Merge FASTQs
        logging.info("Step 4/4: Merging FASTQ files...")
        merged_fastq = merge_fastq_files(
            [hap1_fastq, hap2_fastq],
            output_fastq,
            validate_inputs=False,  # Already validated by simulation
        )
        logging.info("  Merged FASTQ: %s", merged_fastq)

        # Get read counts (import here to avoid circular dependency)
        from .fastq_utils import count_fastq_reads

        reads_hap1 = count_fastq_reads(hap1_fastq)
        reads_hap2 = count_fastq_reads(hap2_fastq)

        logging.info("=" * 70)
        logging.info("SPLIT-SIMULATION COMPLETE")
        logging.info("  Haplotype 1: %d reads", reads_hap1)
        logging.info("  Haplotype 2: %d reads", reads_hap2)
        logging.info("  Total: %d reads", reads_hap1 + reads_hap2)
        logging.info("  Merged FASTQ: %s", merged_fastq)
        logging.info("=" * 70)

        # Optionally copy intermediate files before temp dir is cleaned up
        final_hap1_fastq = hap1_fastq
        final_hap2_fastq = hap2_fastq

        if keep_intermediate:
            import shutil

            output_dir = output_fastq.parent
            final_hap1_fastq = str(output_dir / "hap1_reads.fastq.gz")
            final_hap2_fastq = str(output_dir / "hap2_reads.fastq.gz")
            shutil.copy(hap1_fastq, final_hap1_fastq)
            shutil.copy(hap2_fastq, final_hap2_fastq)
            logging.info("Kept intermediate files: %s, %s", final_hap1_fastq, final_hap2_fastq)

        return DiploidSimulationResult(
            merged_fastq=merged_fastq,
            hap1_fastq=final_hap1_fastq,
            hap2_fastq=final_hap2_fastq,
            hap1_reference=hap1_ref,
            hap2_reference=hap2_ref,
            reads_hap1=reads_hap1,
            reads_hap2=reads_hap2,
        )
