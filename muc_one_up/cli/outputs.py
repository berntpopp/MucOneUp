"""
Output writing functions for MucOneUp CLI.

Single Responsibility: Write simulation outputs (FASTA, structure files, mutated units).
"""

import logging
from pathlib import Path
from typing import Any

from ..exceptions import FileOperationError
from ..fasta_writer import write_fasta
from .config import numbered_filename
from .snps import integrate_snps_unified


def _format_targets(targets: Any) -> list[tuple[int, int]]:
    """Format mutation targets to list of (haplotype, repeat) tuples."""
    return [
        t.as_tuple()
        if hasattr(t, "as_tuple")
        else (int(t.split(",")[0]), int(t.split(",")[1]))
        if isinstance(t, str)
        else tuple(t)
        for t in (targets or [])
    ]


def generate_output_base(input_path: Path, suffix: str) -> str:
    """Generate output base name from input file path.

    Args:
        input_path: Path to input FASTA file.
        suffix: Suffix to append (e.g., '_reads', '_orfs').

    Returns:
        Output base name (e.g., 'sample.001.simulated_reads').
    """
    return input_path.stem + suffix


def write_fasta_outputs(
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
    rng=None,
):
    """
    Write FASTA output files.

    Single Responsibility: FASTA file writing.
    """
    if dual_mutation_mode:
        normal_out = numbered_filename(
            out_dir, out_base, sim_index, "simulated.fa", variant="normal"
        )
        mut_out = numbered_filename(out_dir, out_base, sim_index, "simulated.fa", variant="mut")

        try:
            # Keep original sequences to compute indel shifts
            orig_normal_seqs = [hr.sequence for hr in results]
            orig_mut_seqs = [hr.sequence for hr in mutated_results] if mutated_results else []

            # Apply SNPs to normal variant first
            results, applied_snp_info_normal = integrate_snps_unified(
                args, config, results, rng=rng
            )

            # Translate SNPs for mutated variant
            translated_snps: list[dict[str, Any]] = []
            ref_assembly = getattr(args, "reference_assembly", None)
            if not isinstance(ref_assembly, str):
                ref_assembly = config.get("reference_assembly", "hg38")
            left_len = len(config.get("constants", {}).get(ref_assembly, {}).get("left", ""))

            for hap_idx, snps in applied_snp_info_normal.items():
                if not snps or mutated_results is None or hap_idx >= len(mutated_results):
                    continue
                normal_seq = orig_normal_seqs[hap_idx]
                mut_seq = orig_mut_seqs[hap_idx]
                delta = len(mut_seq) - len(normal_seq)

                norm_lengths = getattr(results[hap_idx], "_repeat_lengths", None) or [
                    len(config["repeats"][ru.symbol]) for ru in results[hap_idx].chain
                ]
                mut_lengths = getattr(mutated_results[hap_idx], "_repeat_lengths", None)
                if mut_lengths is not None and len(mut_lengths) == len(
                    mutated_results[hap_idx].chain
                ):
                    mut_lengths = list(mut_lengths)
                else:
                    mut_lengths = [
                        len(config["repeats"][ru.symbol]) for ru in mutated_results[hap_idx].chain
                    ]
                    if delta != sum(mut_lengths) - sum(norm_lengths):
                        mut_indices = [
                            i
                            for i, ru in enumerate(mutated_results[hap_idx].chain)
                            if getattr(ru, "mutated", False)
                        ]
                        if len(mut_indices) == 1:
                            mut_lengths[mut_indices[0]] += delta - (
                                sum(mut_lengths) - sum(norm_lengths)
                            )
                        elif not mut_indices:
                            hap_targets = [
                                t
                                for t in _format_targets(
                                    mutation_positions or getattr(args, "mutation_targets", None)
                                )
                                if t[0] == hap_idx + 1
                            ]
                            if len(hap_targets) == 1 and 0 <= hap_targets[0][1] - 1 < len(
                                mut_lengths
                            ):
                                mut_lengths[hap_targets[0][1] - 1] += delta - (
                                    sum(mut_lengths) - sum(norm_lengths)
                                )

                repeat_slices: list[tuple[int, int, int, int]] = []
                norm_cur = left_len
                mut_cur = left_len
                for r_i in range(min(len(norm_lengths), len(mut_lengths))):
                    n_len = norm_lengths[r_i]
                    m_len = mut_lengths[r_i]
                    repeat_slices.append((norm_cur, norm_cur + n_len, mut_cur, mut_cur + m_len))
                    norm_cur += n_len
                    mut_cur += m_len

                for snp in snps:
                    pos = snp["position"]
                    new_snp = dict(snp)
                    new_pos = pos
                    if pos < left_len:
                        new_pos = pos
                    elif pos >= left_len + sum(norm_lengths):
                        new_pos = pos + delta
                    else:
                        deleted = False
                        for norm_r_start, norm_r_end, mut_r_start, mut_r_end in repeat_slices:
                            if norm_r_start <= pos < norm_r_end:
                                rel_pos = pos - norm_r_start
                                norm_rep = normal_seq[norm_r_start:norm_r_end]
                                mut_rep = mut_seq[mut_r_start:mut_r_end]
                                r_delta = len(mut_rep) - len(norm_rep)
                                if r_delta == 0:
                                    new_pos = mut_r_start + rel_pos
                                else:
                                    diff_idx = 0
                                    min_len = min(len(norm_rep), len(mut_rep))
                                    while (
                                        diff_idx < min_len
                                        and norm_rep[diff_idx] == mut_rep[diff_idx]
                                    ):
                                        diff_idx += 1
                                    if rel_pos < diff_idx:
                                        new_pos = mut_r_start + rel_pos
                                    elif r_delta < 0 and diff_idx <= rel_pos < diff_idx - r_delta:
                                        deleted = True
                                    else:
                                        new_pos = mut_r_start + rel_pos + r_delta
                                break
                        if deleted:
                            continue

                    new_snp["position"] = new_pos
                    if 0 <= new_pos < len(mut_seq):
                        new_snp["ref_base"] = mut_seq[new_pos].upper()
                    translated_snps.append(new_snp)

            mutated_results, applied_snp_info_mut = integrate_snps_unified(
                args,
                config,
                mutated_results,
                skip_reference_check=True,
                rng=rng,
                predefined_snps=translated_snps,
            )

            write_fasta(
                [hr.sequence for hr in results],
                normal_out,
                comment="Normal sequence (no mutations applied)",
            )
            target_tuples = _format_targets(mutation_positions or args.mutation_targets)
            mutation_comment = f"Mutation Applied: {mutation_pair[1]} (Targets: {target_tuples})"
            write_fasta(
                [hr.sequence for hr in mutated_results],
                mut_out,
                comment=mutation_comment,
            )
            logging.info("Dual FASTA outputs written: %s and %s", normal_out, mut_out)

            return results, mutated_results, applied_snp_info_normal, applied_snp_info_mut
        except Exception as e:
            raise FileOperationError(f"Writing dual FASTA outputs failed: {e}") from e
    else:
        out_file = numbered_filename(out_dir, out_base, sim_index, "simulated.fa")

        try:
            # Create haplotype-specific comments
            haplotype_comments: list[str | None] = [None] * len(results)

            if structure_mutation_info or args.mutation_name:
                mutation_name = (
                    structure_mutation_info["name"]
                    if structure_mutation_info
                    else args.mutation_name
                )
                if structure_mutation_info:
                    mutation_targets = structure_mutation_info["targets"]
                else:
                    mutation_targets = _format_targets(args.mutation_targets or mutation_positions)

                for hap_idx, _rep_idx in mutation_targets:
                    if 1 <= hap_idx <= len(results):
                        haplotype_comments[hap_idx - 1] = (
                            f"Mutation Applied: {mutation_name} (Targets: {mutation_targets})"
                        )

            # Apply SNPs
            results, applied_snp_info = integrate_snps_unified(args, config, results, rng=rng)

            write_fasta(
                [hr.sequence for hr in results],
                out_file,
                prefix="haplotype",
                comments=haplotype_comments,
            )
            logging.info("FASTA output written to %s", out_file)

            return results, None, applied_snp_info, {}
        except Exception as e:
            raise FileOperationError(f"Writing FASTA output failed: {e}") from e


def write_mutated_units(args, out_dir, out_base, sim_index, mutated_units, dual_mutation_mode):
    """Write mutated VNTR unit FASTA if mutation was applied."""
    if args.mutation_name and mutated_units is not None:
        variant_suffix = "mut" if dual_mutation_mode else ""
        mutated_unit_out = numbered_filename(
            out_dir, out_base, sim_index, "mutated_unit.fa", variant=variant_suffix
        )
        try:
            with Path(mutated_unit_out).open("w") as muf:
                for hap_idx, muts in mutated_units.items():
                    for rep_idx, unit_seq in muts:
                        muf.write(f">haplotype_{hap_idx}_repeat_{rep_idx}\n{unit_seq}\n")
            logging.info("Mutated VNTR unit FASTA output written: %s", mutated_unit_out)
        except Exception as e:
            raise FileOperationError(f"Writing mutated VNTR unit FASTA failed: {e}") from e


def write_structure_files(
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
):
    """Write VNTR structure files."""
    if not args.output_structure:
        return

    if dual_mutation_mode:
        normal_struct_out = numbered_filename(
            out_dir, out_base, sim_index, "vntr_structure.txt", variant="normal"
        )
        mut_struct_out = numbered_filename(
            out_dir, out_base, sim_index, "vntr_structure.txt", variant="mut"
        )
        try:
            with Path(normal_struct_out).open("w") as nf:
                nf.write("# Normal sequence (no mutations applied)\n")
                for i, hr in enumerate(results, start=1):
                    chain_str = "-".join(str(ru) for ru in hr.chain)
                    nf.write(f"haplotype_{i}\t{chain_str}\n")

            with Path(mut_struct_out).open("w") as mf:
                target_tuples = _format_targets(mutation_positions or args.mutation_targets)
                mutation_comment = (
                    f"# Mutation Applied: {mutation_pair[1]} (Targets: {target_tuples})\n"
                )
                mf.write(mutation_comment)
                for i, hr in enumerate(mutated_results, start=1):
                    chain_str = "-".join(str(ru) for ru in hr.chain)
                    mf.write(f"haplotype_{i}\t{chain_str}\n")

            logging.info("Structure files written: %s and %s", normal_struct_out, mut_struct_out)
        except Exception as e:
            raise FileOperationError(f"Writing dual structure files failed: {e}") from e
    else:
        struct_out = numbered_filename(out_dir, out_base, sim_index, "vntr_structure.txt")
        try:
            with Path(struct_out).open("w") as struct_fh:
                if structure_mutation_info:
                    struct_fh.write(
                        f"# Mutation Applied: {structure_mutation_info['name']} (Targets: {structure_mutation_info['targets']})\n"
                    )
                elif args.mutation_name:
                    target_tuples = _format_targets(mutation_positions or args.mutation_targets)
                    if target_tuples:
                        struct_fh.write(
                            f"# Mutation Applied: {args.mutation_name} (Targets: {target_tuples})\n"
                        )
                    else:
                        struct_fh.write(
                            f"# Mutation Applied: {args.mutation_name} (Target: random)\n"
                        )

                for i, hr in enumerate(results, start=1):
                    chain_str = "-".join(str(ru) for ru in hr.chain)
                    struct_fh.write(f"haplotype_{i}\t{chain_str}\n")

            logging.info("Structure file written to %s", struct_out)
        except Exception as e:
            raise FileOperationError(f"Writing structure file failed: {e}") from e
