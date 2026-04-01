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
        args: SimulationOptions instance (or compatible namespace) with CLI arguments
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

    # Build read source tracker(s) if requested
    source_tracker = None
    source_tracker_mut = None
    if getattr(args, "track_read_source", False) and getattr(args, "simulate_reads", None):
        from pathlib import Path

        from ..read_simulator.source_tracking import ReadSourceTracker

        # Build repeat chains dict (1-based haplotype keys)
        # results is list[tuple[str, list[str]]] — (sequence, chain)
        repeat_chains = {}
        for i, result in enumerate(results):
            _sequence, chain = result
            repeat_chains[i + 1] = chain

        # Get left constant length from config
        ref_assembly = getattr(args, "reference_assembly", None) or config.get(
            "reference_assembly", "hg38"
        )
        left_const = config.get("constants", {}).get(ref_assembly, {}).get("left", "")
        left_const_len = len(left_const)

        # Get repeats dict
        repeats_dict = config.get("repeats", {})

        # Get mutation info
        mutation_name_str = getattr(args, "mutation_name", None)
        mut_positions = []
        mut_name = None
        if mutation_positions:
            mut_positions = mutation_positions
            if mutation_name_str and mutation_name_str != "normal":
                parts = mutation_name_str.split(",")
                for p in parts:
                    if p != "normal":
                        mut_name = p
                        break

        # Get SNP info (0-based haplotype indices)
        snp_info_dict = {}
        if applied_snp_info_normal:
            for i, snp_list in enumerate(applied_snp_info_normal):
                if snp_list:
                    snp_info_dict[i] = snp_list

        if dual_mutation_mode:
            # Build separate trackers: normal (no mutation) and mutated (with mutation)
            source_tracker = ReadSourceTracker(
                repeat_chains=repeat_chains,
                repeats_dict=repeats_dict,
                left_const_len=left_const_len,
                mutation_positions=[],
                mutation_name=None,
                snp_info=snp_info_dict if snp_info_dict else None,
            )

            # Build mutated chains from mutated_results
            mut_chains = {}
            if mutated_results:
                for i, result in enumerate(mutated_results):
                    _sequence, chain = result
                    mut_chains[i + 1] = chain

            # Get mutated SNP info
            mut_snp_info_dict: dict[int, list[dict[str, object]]] = {}
            if applied_snp_info_mut:
                for i, snp_list in enumerate(applied_snp_info_mut):
                    if snp_list:
                        mut_snp_info_dict[i] = snp_list

            source_tracker_mut = ReadSourceTracker(
                repeat_chains=mut_chains if mut_chains else repeat_chains,
                repeats_dict=repeats_dict,
                left_const_len=left_const_len,
                mutation_positions=mut_positions,
                mutation_name=mut_name,
                snp_info=mut_snp_info_dict if mut_snp_info_dict else None,
            )

            # Write coordinate maps for both variants
            normal_coord_path = str(
                Path(out_dir) / f"{out_base}.{sim_index:03d}.normal.repeat_coordinates.tsv"
            )
            source_tracker.write_coordinate_map(normal_coord_path)
            logging.info("Normal repeat coordinate map written: %s", normal_coord_path)

            mut_coord_path = str(
                Path(out_dir) / f"{out_base}.{sim_index:03d}.mut.repeat_coordinates.tsv"
            )
            source_tracker_mut.write_coordinate_map(mut_coord_path)
            logging.info("Mutated repeat coordinate map written: %s", mut_coord_path)
        else:
            source_tracker = ReadSourceTracker(
                repeat_chains=repeat_chains,
                repeats_dict=repeats_dict,
                left_const_len=left_const_len,
                mutation_positions=mut_positions,
                mutation_name=mut_name,
                snp_info=snp_info_dict if snp_info_dict else None,
            )

            coord_map_path = str(
                Path(out_dir) / f"{out_base}.{sim_index:03d}.repeat_coordinates.tsv"
            )
            source_tracker.write_coordinate_map(coord_map_path)
            logging.info("Repeat coordinate map written: %s", coord_map_path)

    # Read simulation
    run_read_simulation(
        args,
        config,
        out_dir,
        out_base,
        sim_index,
        dual_mutation_mode,
        source_tracker=source_tracker,
        source_tracker_mut=source_tracker_mut,
    )

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
