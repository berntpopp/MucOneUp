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
    """
    Run a single simulation iteration with all processing steps.

    Single Responsibility: Orchestrate one complete simulation iteration.
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
