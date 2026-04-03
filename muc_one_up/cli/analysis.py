"""
Analysis functions for MucOneUp CLI.

Single Responsibility: Run optional analyses (ORF prediction, read simulation, statistics).
"""

import json
import logging
from pathlib import Path

from ..exceptions import FileOperationError, ReadSimulationError
from ..provenance import collect_provenance_metadata
from ..simulation_statistics import (
    generate_simulation_statistics,
    write_statistics_report,
)
from .config import numbered_filename


def _run_toxic_detection(orf_path: str, config: dict) -> dict:
    """Run toxic protein detection on an ORF FASTA file.

    Extracts detection parameters from config and delegates to scan_orf_fasta.
    Returns the stats dict from the scanner.
    """
    from ..toxic_protein_detector import scan_orf_fasta

    ref_assembly = config.get("reference_assembly", "hg38")
    assembly_constants = config.get("constants", {}).get(ref_assembly, {})
    left_const_val = assembly_constants.get("left")
    right_const_val = assembly_constants.get("right")

    toxic_config = config.get("toxic_protein_detection", {})
    detection_kwargs = {
        "consensus": toxic_config.get("consensus_motif", "RCHLGPGHQAGPGLHR"),
        "identity_threshold": toxic_config.get("identity_threshold", 0.8),
        "key_residues": toxic_config.get("key_residues", ["R", "C", "H"]),
        "expected_repeat_count": toxic_config.get("expected_repeat_count", 10),
        "w_repeat": toxic_config.get("weights", {}).get("repeat", 0.6),
        "w_composition": toxic_config.get("weights", {}).get("composition", 0.4),
        "toxic_detection_cutoff": toxic_config.get("toxic_cutoff", 0.5),
    }

    return scan_orf_fasta(
        orf_path,
        left_const=left_const_val,
        right_const=right_const_val,
        **detection_kwargs,
    )


def run_orf_prediction(
    args, config, out_dir, out_base, sim_index, results, mutated_results, dual_mutation_mode
):
    """Run ORF prediction and toxic protein detection."""
    if not args.output_orfs:
        return

    from ..translate import run_orf_finder_in_memory

    if dual_mutation_mode:
        normal_orf_out = numbered_filename(
            out_dir, out_base, sim_index, "orfs.fa", variant="normal"
        )
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
            normal_stats = _run_toxic_detection(normal_orf_out, config)
            mut_stats = _run_toxic_detection(mut_orf_out, config)

            stats_file_normal = numbered_filename(
                out_dir, out_base, sim_index, "orf_stats.txt", variant="normal"
            )
            stats_file_mut = numbered_filename(
                out_dir, out_base, sim_index, "orf_stats.txt", variant="mut"
            )

            with Path(stats_file_normal).open("w") as nf:
                json.dump(normal_stats, nf, indent=4)
            with Path(stats_file_mut).open("w") as mf:
                json.dump(mut_stats, mf, indent=4)

            logging.info(
                "Toxic protein detection stats written: %s and %s",
                stats_file_normal,
                stats_file_mut,
            )
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
            stats = _run_toxic_detection(orf_out, config)
            stats_file = numbered_filename(out_dir, out_base, sim_index, "orf_stats.txt")

            with Path(stats_file).open("w") as sf:
                json.dump(stats, sf, indent=4)

            logging.info("Toxic protein detection stats written: %s", stats_file)
        except ImportError as e:
            raise FileOperationError(f"Failed to import toxic_protein_detector module: {e}") from e


def run_read_simulation(
    args,
    config,
    out_dir,
    out_base,
    sim_index,
    dual_mutation_mode,
    source_tracker=None,
    source_tracker_mut=None,
):
    """Run read simulation pipeline if requested."""
    if not args.simulate_reads:
        return

    from ..read_simulation import simulate_reads as simulate_reads_pipeline

    if "read_simulation" not in config:
        config["read_simulation"] = {}
    config["read_simulation"]["simulator"] = args.simulate_reads

    simulator_type = args.simulate_reads
    simulator_name = "Oxford Nanopore" if simulator_type == "ont" else "Illumina"

    if dual_mutation_mode:
        normal_fa = numbered_filename(
            out_dir, out_base, sim_index, "simulated.fa", variant="normal"
        )
        mut_fa = numbered_filename(out_dir, out_base, sim_index, "simulated.fa", variant="mut")

        try:
            logging.info(
                "Starting %s read simulation pipeline for iteration %d (normal variant).",
                simulator_name,
                sim_index,
            )
            simulate_reads_pipeline(config, normal_fa, source_tracker=source_tracker)
            logging.info(
                "%s read simulation pipeline completed for iteration %d (normal variant).",
                simulator_name,
                sim_index,
            )
        except Exception as e:
            raise ReadSimulationError(
                f"{simulator_name} read simulation pipeline (normal variant) failed: {e}"
            ) from e

        try:
            logging.info(
                "Starting %s read simulation pipeline for iteration %d (mutated variant).",
                simulator_name,
                sim_index,
            )
            simulate_reads_pipeline(
                config, mut_fa, source_tracker=source_tracker_mut or source_tracker
            )
            logging.info(
                "%s read simulation pipeline completed for iteration %d (mutated variant).",
                simulator_name,
                sim_index,
            )
        except Exception as e:
            raise ReadSimulationError(
                f"{simulator_name} read simulation pipeline (mutated variant) failed: {e}"
            ) from e
    else:
        sim_fa = numbered_filename(out_dir, out_base, sim_index, "simulated.fa")

        try:
            logging.info(
                "Starting %s read simulation pipeline for iteration %d.", simulator_name, sim_index
            )
            simulate_reads_pipeline(config, sim_fa, source_tracker=source_tracker)
            logging.info(
                "%s read simulation pipeline completed for iteration %d.", simulator_name, sim_index
            )
        except Exception as e:
            raise ReadSimulationError(
                f"{simulator_name} read simulation pipeline failed: {e}"
            ) from e


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
    """Generate and write simulation statistics with provenance metadata."""
    from typing import Any

    vntr_coverage_stats: dict[str, Any] = {}

    # Collect provenance metadata for reproducibility
    # Defense-in-depth: collect_provenance_metadata() handles errors internally,
    # but we add outer protection in case of unexpected issues
    try:
        provenance_info = collect_provenance_metadata(
            config=config,
            start_time=iteration_start,
            end_time=iteration_end,
        )
    except Exception as e:
        # Graceful degradation: simulation continues even if provenance fails
        logging.error(f"Unexpected error collecting provenance metadata: {e}", exc_info=True)
        provenance_info = None

    if dual_mutation_mode:
        normal_stats_report = generate_simulation_statistics(
            start_time=iteration_start,
            end_time=iteration_end,
            simulation_results=results,
            config=config,
            mutation_info={"mutation_name": "normal"},
            vntr_coverage=vntr_coverage_stats,
            applied_snp_info=applied_snp_info_normal,
            provenance_info=provenance_info,
        )
        mutated_stats_report = generate_simulation_statistics(
            start_time=iteration_start,
            end_time=iteration_end,
            simulation_results=mutated_results,
            config=config,
            mutation_info={"mutation_name": mutation_pair[1]},
            vntr_coverage=vntr_coverage_stats,
            applied_snp_info=applied_snp_info_mut,
            provenance_info=provenance_info,
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
            provenance_info=provenance_info,
        )

        stats_output_file = numbered_filename(out_dir, out_base, sim_index, "simulation_stats.json")
        write_statistics_report(stats_report, stats_output_file)


def run_orf_analysis_standalone(
    input_fasta: str,
    out_dir: str,
    out_base: str,
    orf_min_aa: int,
    orf_aa_prefix: str | None,
    left_const: str | None,
    right_const: str | None,
) -> None:
    """Run ORF prediction and toxic protein detection on a FASTA file.

    Single implementation for standalone ORF analysis, used by the
    ``analyze orfs`` CLI command.
    """
    from ..read_simulator.utils.common_utils import run_command

    orf_output = Path(out_dir) / f"{out_base}.orfs.fa"

    orf_filename = f"{out_base}.orfs.fa"
    cmd = [
        "orfipy",
        input_fasta,
        "--outdir",
        str(out_dir),
        "--pep",
        orf_filename,
        "--min",
        str(orf_min_aa * 3),
        "--start",
        "ATG",
    ]

    try:
        run_command(cmd, capture=True)
    except Exception:
        logging.error("orfipy failed for %s", input_fasta)
        return

    logging.info("ORF prediction completed: %s", orf_output)

    # Filter by amino acid prefix if specified
    if orf_aa_prefix and orf_output.exists():
        from Bio import SeqIO

        filtered_orfs = []
        total_orfs = 0

        for record in SeqIO.parse(str(orf_output), "fasta"):
            total_orfs += 1
            if str(record.seq).startswith(orf_aa_prefix):
                filtered_orfs.append(record)

        SeqIO.write(filtered_orfs, str(orf_output), "fasta")
        logging.info(
            "Filtered ORFs by prefix '%s': %d/%d ORFs retained",
            orf_aa_prefix,
            len(filtered_orfs),
            total_orfs,
        )

    # Toxic protein detection
    if orf_output.exists():
        from ..toxic_protein_detector import scan_orf_fasta

        stats = scan_orf_fasta(str(orf_output), left_const=left_const, right_const=right_const)

        stats_file = Path(out_dir) / f"{out_base}.orf_stats.json"
        with stats_file.open("w") as f:
            json.dump(stats, f, indent=4)

        logging.info("Toxic protein stats written: %s", stats_file)
