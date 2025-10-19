"""
Orchestration functions for MucOneUp CLI.

Single Responsibility: Orchestrate complete simulation iterations.
"""

import logging
import time

from .analysis import (
    run_orf_prediction,
    run_read_simulation,
    write_simulation_statistics,
)
from .haplotypes import generate_haplotypes
from .mutations import apply_mutation_pipeline
from .outputs import (
    write_fasta_outputs,
    write_mutated_units,
    write_structure_files,
)


def run_single_simulation_iteration(
    args,
    config,
    out_dir,
    out_base,
    sim_index,
    fixed_conf,
    predefined_chains,
    dual_mutation_mode,
    mutation_pair,
    structure_mutation_info,
):
    """Run complete simulation iteration with all processing steps.

    Orchestrates a single simulation iteration from haplotype generation
    through mutation application, FASTA output, ORF prediction, read
    simulation, and statistics generation. Handles both normal and dual
    mutation modes.

    Args:
        args: Namespace containing all CLI arguments from argparse
        config: Configuration dictionary from load_config()
        out_dir: Output directory path for all generated files
        out_base: Base filename for outputs (without iteration suffix)
        sim_index: Current iteration index (for numbered outputs)
        fixed_conf: Fixed length configuration or "from_structure" flag
        predefined_chains: Predefined repeat chains from structure file (or None)
        dual_mutation_mode: Boolean indicating dual simulation mode (normal + mutated)
        mutation_pair: Tuple of (normal_name, mutated_name) for dual mode
        structure_mutation_info: Mutation metadata from structure file comments

    Returns:
        None (all outputs written to files)

    Raises:
        SimulationError: If haplotype generation fails
        ValidationError: If configuration or parameters are invalid
        FileOperationError: If file writing fails

    Side Effects:
        Creates multiple output files:
        - FASTA files (*.fa or *.normal.fa + *.mut.fa)
        - Structure files (*.structure.txt)
        - Mutated units files (*.mutated_units.txt)
        - Statistics JSON (*.simulation_stats.json)
        - Optional: ORF peptides (*.pep.fa)
        - Optional: Read simulation BAM files

    Note:
        Single Responsibility: Orchestrate one complete simulation iteration.
        This function delegates to specialized modules (haplotypes, mutations,
        outputs, analysis) following SOLID principles.
    """
    iteration_start = time.time()

    # Generate haplotypes
    results = generate_haplotypes(args, config, fixed_conf, predefined_chains)
    logging.info("Haplotype simulation completed successfully for iteration %d.", sim_index)

    # Apply mutations
    results, mutated_results, mutated_units, mutation_positions = apply_mutation_pipeline(
        args, config, results, args.mutation_name, dual_mutation_mode, mutation_pair
    )

    # Write FASTA outputs (includes SNP integration)
    results, mutated_results, applied_snp_info_normal, applied_snp_info_mut = write_fasta_outputs(
        args,
        config,
        out_dir,
        out_base,
        sim_index,
        results,
        mutated_results,
        dual_mutation_mode,
        mutation_pair,
        mutation_positions,
        structure_mutation_info,
    )

    # Write mutated units
    write_mutated_units(args, out_dir, out_base, sim_index, mutated_units, dual_mutation_mode)

    # Write structure files
    write_structure_files(
        args,
        out_dir,
        out_base,
        sim_index,
        results,
        mutated_results,
        dual_mutation_mode,
        mutation_pair,
        mutation_positions,
        structure_mutation_info,
    )

    # ORF prediction and toxic detection
    run_orf_prediction(
        args, config, out_dir, out_base, sim_index, results, mutated_results, dual_mutation_mode
    )

    # Read simulation
    run_read_simulation(args, config, out_dir, out_base, sim_index, dual_mutation_mode)

    # Statistics
    iteration_end = time.time()
    write_simulation_statistics(
        args,
        config,
        out_dir,
        out_base,
        sim_index,
        iteration_start,
        iteration_end,
        results,
        mutated_results,
        dual_mutation_mode,
        mutation_pair,
        applied_snp_info_normal,
        applied_snp_info_mut,
    )
