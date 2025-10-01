"""
Analysis functions for MucOneUp CLI.

Single Responsibility: Run optional analyses (ORF prediction, read simulation, statistics).
"""

import json
import logging
from pathlib import Path

from ..exceptions import FileOperationError, ReadSimulationError
from ..read_simulation import simulate_reads as simulate_reads_pipeline
from ..simulation_statistics import (
    generate_simulation_statistics,
    write_statistics_report,
)
from ..translate import run_orf_finder_in_memory
from .config import numbered_filename


def run_orf_prediction(
    args, config, out_dir, out_base, sim_index, results, mutated_results, dual_mutation_mode
):
    """Run ORF prediction and toxic protein detection."""
    if not args.output_orfs:
        return

    if dual_mutation_mode:
        normal_orf_out = numbered_filename(out_dir, out_base, sim_index, "orfs.fa", variant="normal")
        mut_orf_out = numbered_filename(out_dir, out_base, sim_index, "orfs.fa", variant="mut")

        try:
            run_orf_finder_in_memory(
                results,
                output_pep=normal_orf_out,
                orf_min_aa=args.orf_min_aa,
                required_prefix=args.orf_aa_prefix,
            )
            run_orf_finder_in_memory(
                mutated_results,
                output_pep=mut_orf_out,
                orf_min_aa=args.orf_min_aa,
                required_prefix=args.orf_aa_prefix,
            )
            logging.info(
                "ORF finding completed; peptide FASTA outputs written: %s and %s",
                normal_orf_out,
                mut_orf_out,
            )
        except Exception as e:
            raise FileOperationError(f"ORF finding failed (dual mode): {e}") from e

        # Toxic protein detection
        try:
            from ..toxic_protein_detector import scan_orf_fasta

            left_const_val = config.get("constants", {}).get("left")
            right_const_val = config.get("constants", {}).get("right")

            normal_stats = scan_orf_fasta(normal_orf_out, left_const=left_const_val, right_const=right_const_val)
            mut_stats = scan_orf_fasta(mut_orf_out, left_const=left_const_val, right_const=right_const_val)

            stats_file_normal = numbered_filename(out_dir, out_base, sim_index, "orf_stats.txt", variant="normal")
            stats_file_mut = numbered_filename(out_dir, out_base, sim_index, "orf_stats.txt", variant="mut")

            with Path(stats_file_normal).open("w") as nf:
                json.dump(normal_stats, nf, indent=4)
            with Path(stats_file_mut).open("w") as mf:
                json.dump(mut_stats, mf, indent=4)

            logging.info("Toxic protein detection stats written: %s and %s", stats_file_normal, stats_file_mut)
        except ImportError as e:
            raise FileOperationError(f"Failed to import toxic_protein_detector module: {e}") from e
    else:
        orf_out = numbered_filename(out_dir, out_base, sim_index, "orfs.fa")

        try:
            run_orf_finder_in_memory(
                results,
                output_pep=orf_out,
                orf_min_aa=args.orf_min_aa,
                required_prefix=args.orf_aa_prefix,
            )
            logging.info("ORF finding completed; peptide FASTA written to %s", orf_out)
        except Exception as e:
            raise FileOperationError(f"ORF finding failed: {e}") from e

        # Toxic protein detection
        try:
            from ..toxic_protein_detector import scan_orf_fasta

            left_const_val = config.get("constants", {}).get("left")
            right_const_val = config.get("constants", {}).get("right")

            stats = scan_orf_fasta(orf_out, left_const=left_const_val, right_const=right_const_val)
            stats_file = numbered_filename(out_dir, out_base, sim_index, "orf_stats.txt")

            with Path(stats_file).open("w") as sf:
                json.dump(stats, sf, indent=4)

            logging.info("Toxic protein detection stats written: %s", stats_file)
        except ImportError as e:
            raise FileOperationError(f"Failed to import toxic_protein_detector module: {e}") from e


def run_read_simulation(
    args, config, out_dir, out_base, sim_index, dual_mutation_mode
):
    """Run read simulation pipeline if requested."""
    if not args.simulate_reads:
        return

    if "read_simulation" not in config:
        config["read_simulation"] = {}
    config["read_simulation"]["simulator"] = args.simulate_reads

    simulator_type = args.simulate_reads
    simulator_name = "Oxford Nanopore" if simulator_type == "ont" else "Illumina"

    if dual_mutation_mode:
        normal_fa = numbered_filename(out_dir, out_base, sim_index, "simulated.fa", variant="normal")
        mut_fa = numbered_filename(out_dir, out_base, sim_index, "simulated.fa", variant="mut")

        try:
            logging.info(
                "Starting %s read simulation pipeline for iteration %d (normal variant).",
                simulator_name,
                sim_index,
            )
            simulate_reads_pipeline(config, normal_fa)
            logging.info(
                "%s read simulation pipeline completed for iteration %d (normal variant).",
                simulator_name,
                sim_index,
            )
        except Exception as e:
            raise ReadSimulationError(f"{simulator_name} read simulation pipeline (normal variant) failed: {e}") from e

        try:
            logging.info(
                "Starting %s read simulation pipeline for iteration %d (mutated variant).",
                simulator_name,
                sim_index,
            )
            simulate_reads_pipeline(config, mut_fa)
            logging.info(
                "%s read simulation pipeline completed for iteration %d (mutated variant).",
                simulator_name,
                sim_index,
            )
        except Exception as e:
            raise ReadSimulationError(f"{simulator_name} read simulation pipeline (mutated variant) failed: {e}") from e
    else:
        sim_fa = numbered_filename(out_dir, out_base, sim_index, "simulated.fa")

        try:
            logging.info(
                "Starting %s read simulation pipeline for iteration %d.", simulator_name, sim_index
            )
            simulate_reads_pipeline(config, sim_fa)
            logging.info("%s read simulation pipeline completed for iteration %d.", simulator_name, sim_index)
        except Exception as e:
            raise ReadSimulationError(f"{simulator_name} read simulation pipeline failed: {e}") from e


def write_simulation_statistics(
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
):
    """Generate and write simulation statistics."""
    vntr_coverage_stats = {}

    if dual_mutation_mode:
        normal_stats_report = generate_simulation_statistics(
            start_time=iteration_start,
            end_time=iteration_end,
            simulation_results=results,
            config=config,
            mutation_info={"mutation_name": "normal"},
            vntr_coverage=vntr_coverage_stats,
            applied_snp_info=applied_snp_info_normal,
        )
        mutated_stats_report = generate_simulation_statistics(
            start_time=iteration_start,
            end_time=iteration_end,
            simulation_results=mutated_results,
            config=config,
            mutation_info={"mutation_name": mutation_pair[1]},
            vntr_coverage=vntr_coverage_stats,
            applied_snp_info=applied_snp_info_mut,
        )

        stats_file_normal = numbered_filename(
            out_dir, out_base, sim_index, "simulation_stats.json", variant="normal"
        )
        stats_file_mut = numbered_filename(
            out_dir, out_base, sim_index, "simulation_stats.json", variant="mut"
        )

        write_statistics_report(normal_stats_report, stats_file_normal)
        write_statistics_report(mutated_stats_report, stats_file_mut)
    else:
        mutation_info = {}
        if args.mutation_name:
            mutation_info = {
                "mutation_name": args.mutation_name,
                "mutation_targets": args.mutation_targets,
            }

        stats_report = generate_simulation_statistics(
            start_time=iteration_start,
            end_time=iteration_end,
            simulation_results=results,
            config=config,
            mutation_info=mutation_info,
            vntr_coverage=vntr_coverage_stats,
            applied_snp_info=applied_snp_info_normal,
        )

        stats_output_file = numbered_filename(out_dir, out_base, sim_index, "simulation_stats.json")
        write_statistics_report(stats_report, stats_output_file)
