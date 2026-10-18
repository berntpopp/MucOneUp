"""
Analysis functions for MucOneUp CLI.

Single Responsibility: Run optional analyses (standalone ORF analysis, statistics).
"""

import json
import logging
import re
from pathlib import Path
from typing import Any

from ..provenance import collect_provenance_metadata
from ..simulation_statistics import (
    generate_simulation_statistics,
    write_statistics_report,
)
from .config import numbered_filename


def _add_companion_fields(
    report, args, config, results, mutation_positions, mutation_name, mutated_units=None
):
    """Add read-source companion fields to a stats report if tracking was requested.

    These fields let ``ReadSourceTracker.from_companion_files`` rebuild the
    tracker when reads are simulated later with ``reads ... --track-read-source``.
    """
    if getattr(args, "track_read_source", False) is not True:
        return
    from ..read_simulator.source_tracking import companion_stats_fields

    fields = companion_stats_fields(
        results, config, mutation_positions, mutation_name, mutated_units=mutated_units
    )
    report.update(fields)


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
    mutation_positions=None,
    mutated_units=None,
):
    """Generate and write simulation statistics with provenance metadata.

    When ``args.track_read_source`` is set, the report also carries the
    companion fields needed to reconstruct read source tracking.
    """
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
            mutated_units=mutated_units,
        )
        _add_companion_fields(
            normal_stats_report, args, config, results, None, None, mutated_units=None
        )
        _add_companion_fields(
            mutated_stats_report,
            args,
            config,
            mutated_results,
            mutation_positions,
            mutation_pair[1],
            mutated_units=mutated_units,
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
            mutated_units=mutated_units,
        )
        _add_companion_fields(
            stats_report,
            args,
            config,
            results,
            mutation_positions,
            args.mutation_name,
            mutated_units=mutated_units,
        )

        stats_output_file = numbered_filename(out_dir, out_base, sim_index, "simulation_stats.json")
        write_statistics_report(stats_report, stats_output_file)


def _is_dna(seq: str | None) -> bool:
    """Check if sequence is a DNA string."""
    if not seq or len(seq) < 3:
        return False
    return set(seq.upper()) <= set("ACGTN")


def _translate_dna_flanks(
    left_const: str | None,
    right_const: str | None,
    orf_fasta_path: Path | None = None,
) -> tuple[str | list[str] | None, str | list[str] | None]:
    """Translate DNA flank constants to in-frame peptide sequences for ORF scanning.

    If left_const/right_const are already peptide strings (or None), they are returned as-is.
    If DNA, candidate reading frame translations are generated and matched against predicted
    peptides in orf_fasta_path (or orfipy header coordinates).
    """
    from Bio import SeqIO
    from Bio.Seq import Seq

    left_pep: str | list[str] | None = left_const
    right_pep: str | list[str] | None = right_const

    candidate_left_prefixes: list[str] = []
    if left_const and _is_dna(left_const):
        left_dna = left_const.upper()
        left_len = len(left_dna)
        f0 = str(Seq(left_dna[: (left_len // 3) * 3]).translate(to_stop=True))
        if f0:
            candidate_left_prefixes.append(f0)
        for f in (0, 1, 2):
            sub = left_dna[f : f + ((left_len - f) // 3) * 3]
            if len(sub) >= 3:
                full_t = str(Seq(sub).translate())
                to_stop = str(Seq(sub).translate(to_stop=True))
                if to_stop:
                    candidate_left_prefixes.append(to_stop)
                if "*" in full_t:
                    last_stop = full_t.rfind("*")
                    after_stop = full_t[last_stop + 1 :]
                    m_idx = after_stop.find("M")
                    if m_idx != -1:
                        candidate_left_prefixes.append(after_stop[m_idx:])
                else:
                    m_idx = full_t.find("M")
                    if m_idx != -1:
                        candidate_left_prefixes.append(full_t[m_idx:])
        candidate_left_prefixes = sorted(
            [p for p in dict.fromkeys(candidate_left_prefixes) if len(p) >= 1],
            key=len,
            reverse=True,
        )
        left_pep = f0 if f0 else (candidate_left_prefixes[0] if candidate_left_prefixes else None)

    candidate_right_suffixes: list[str] = []
    if right_const and _is_dna(right_const):
        right_dna = right_const.upper()
        right_len = len(right_dna)
        f0_r = str(Seq(right_dna[: (right_len // 3) * 3]).translate(to_stop=True))
        if f0_r:
            candidate_right_suffixes.append(f0_r)
        for f in (0, 1, 2):
            sub = right_dna[f : f + ((right_len - f) // 3) * 3]
            if len(sub) >= 3:
                to_stop = str(Seq(sub).translate(to_stop=True))
                if to_stop:
                    candidate_right_suffixes.append(to_stop)
        candidate_right_suffixes = sorted(
            [s for s in dict.fromkeys(candidate_right_suffixes) if len(s) >= 1],
            key=len,
            reverse=True,
        )
        right_pep = (
            f0_r if f0_r else (candidate_right_suffixes[0] if candidate_right_suffixes else None)
        )

    # If we have an ORF output file, check actual peptide records to match the active frame
    if orf_fasta_path and orf_fasta_path.exists():
        try:
            records = list(SeqIO.parse(str(orf_fasta_path), "fasta"))
            if left_const and _is_dna(left_const):
                for rec in records:
                    m = re.search(r"\[(\d+)-(\d+)\]\(([+-])\)", rec.description)
                    if m and m.group(3) == "+":
                        start_nt = int(m.group(1))
                        if start_nt < len(left_const):
                            aa_count = (len(left_const) - start_nt) // 3
                            seq_str = str(rec.seq)
                            if 0 < aa_count <= len(seq_str):
                                candidate_left_prefixes.append(seq_str[:aa_count])
                candidate_left_prefixes = sorted(
                    [p for p in dict.fromkeys(candidate_left_prefixes) if len(p) >= 1],
                    key=len,
                    reverse=True,
                )

            matched_left: list[str] = []
            matched_right: list[str] = []
            for rec in records:
                seq_str = str(rec.seq)
                if left_const and _is_dna(left_const):
                    for pfx in candidate_left_prefixes:
                        if seq_str.startswith(pfx):
                            if pfx not in matched_left:
                                matched_left.append(pfx)
                            break
                elif left_const and left_const not in matched_left:
                    matched_left.append(left_const)

                if right_const and _is_dna(right_const):
                    for sfx in candidate_right_suffixes:
                        if seq_str.endswith(sfx):
                            if sfx not in matched_right:
                                matched_right.append(sfx)
                            break
                elif right_const and right_const not in matched_right:
                    matched_right.append(right_const)

            if matched_left:
                matched_left = sorted(matched_left, key=len, reverse=True)
                left_pep = matched_left[0] if len(matched_left) == 1 else matched_left
            if matched_right:
                matched_right = sorted(matched_right, key=len, reverse=True)
                right_pep = matched_right[0] if len(matched_right) == 1 else matched_right
        except Exception as e:
            logging.debug("Could not inspect ORF records for flank matching: %s", e)

    return left_pep, right_pep


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
    from ..exceptions import ExternalToolError
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
    except ExternalToolError as e:
        logging.error("orfipy failed for %s: %s", input_fasta, e)
        return
    except FileNotFoundError:
        logging.error("orfipy not found — is it installed?")
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

        left_pep, right_pep = _translate_dna_flanks(left_const, right_const, orf_output)
        stats = scan_orf_fasta(str(orf_output), left_const=left_pep, right_const=right_pep)

        stats_file = Path(out_dir) / f"{out_base}.orf_stats.json"
        with stats_file.open("w") as f:
            json.dump(stats, f, indent=4)

        logging.info("Toxic protein stats written: %s", stats_file)
