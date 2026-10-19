"""
Output writing functions for MucOneUp CLI.

Single Responsibility: Write simulation outputs (FASTA, structure files, mutated units).
"""

import logging
from pathlib import Path
from typing import Any

from ..exceptions import FileOperationError
from ..fasta_writer import write_fasta
from ..type_defs import HaplotypeResult
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


def _indel_match_score(norm_rep: str, cand_rep: str) -> tuple[int, int]:
    """Score how well cand_rep matches norm_rep as a single-indel/mutation derivative."""
    lcp = 0
    min_len = min(len(norm_rep), len(cand_rep))
    while lcp < min_len and norm_rep[lcp] == cand_rep[lcp]:
        lcp += 1
    lcs = 0
    while (
        lcs < min_len - lcp
        and norm_rep[len(norm_rep) - 1 - lcs] == cand_rep[len(cand_rep) - 1 - lcs]
    ):
        lcs += 1
    covered = lcp + lcs
    return (1 if covered == min_len else 0, covered)


def _infer_target_repeat_lengths(
    normal_seq: str,
    mut_seq: str,
    left_len: int,
    norm_lengths: list[int],
    mut_lengths: list[int],
    targets: list[int],
    delta: int,
) -> None:
    """Infer mutated repeat lengths in-place when _repeat_lengths and mutated_units are absent."""
    sorted_targets = sorted(set(targets))
    if len(sorted_targets) == 1:
        mut_lengths[sorted_targets[0]] += delta - (sum(mut_lengths) - sum(norm_lengths))
        return

    for idx, t_i in enumerate(sorted_targets):
        if idx == len(sorted_targets) - 1:
            rem_disc = delta - (sum(mut_lengths) - sum(norm_lengths))
            mut_lengths[t_i] += rem_disc
            break

        t_next = sorted_targets[idx + 1]
        norm_start = left_len + sum(norm_lengths[:t_i])
        norm_rep = normal_seq[norm_start : norm_start + norm_lengths[t_i]]
        unmut_start = norm_start + norm_lengths[t_i]
        unmut_end = left_len + sum(norm_lengths[:t_next])
        unmut_seq = normal_seq[unmut_start:unmut_end]
        next_norm_rep = normal_seq[unmut_end : unmut_end + norm_lengths[t_next]]

        mut_start = left_len + sum(mut_lengths[:t_i])
        max_cand_len = len(mut_seq) - mut_start - len(unmut_seq)
        best_len: int | None = None
        best_score = (-1, -1, -1, -1)

        for cand_len in range(1, max_cand_len + 1):
            if mut_seq[mut_start + cand_len : mut_start + cand_len + len(unmut_seq)] != unmut_seq:
                continue
            cand_rep = mut_seq[mut_start : mut_start + cand_len]
            s1_full, s1_cov = _indel_match_score(norm_rep, cand_rep)
            next_mut_start = mut_start + cand_len + len(unmut_seq)
            if idx + 1 == len(sorted_targets) - 1:
                rem_after_next = sum(norm_lengths[t_next + 1 :]) + (
                    len(normal_seq) - left_len - sum(norm_lengths)
                )
                next_cand_end = len(mut_seq) - rem_after_next
                next_cand_rep = mut_seq[next_mut_start : max(next_mut_start, next_cand_end)]
                s2_full, s2_cov = _indel_match_score(next_norm_rep, next_cand_rep)
            else:
                next_prefix = mut_seq[next_mut_start : next_mut_start + len(next_norm_rep)]
                s2_full, s2_cov = _indel_match_score(next_norm_rep, next_prefix)
            score = (s1_full + s2_full, s1_cov + s2_cov, s1_full, s1_cov)
            if score > best_score:
                best_score = score
                best_len = cand_len

        if best_len is not None:
            mut_lengths[t_i] = best_len
        else:
            rem_targets = len(sorted_targets) - idx
            rem_disc = delta - (sum(mut_lengths) - sum(norm_lengths))
            mut_lengths[t_i] += rem_disc // rem_targets


def _translate_snps_to_mutated(
    args: Any,
    config: dict[str, Any],
    normal_results: list[HaplotypeResult],
    mutated_results: list[HaplotypeResult],
    applied_snp_info_normal: dict[int, list[dict[str, Any]]],
    mutation_positions: Any = None,
    orig_normal_seqs: list[str] | None = None,
    orig_mut_seqs: list[str] | None = None,
    mutated_units: dict[int, list[tuple[int, str]]] | None = None,
) -> list[dict[str, Any]]:
    """Translate applied SNP coordinates from normal sequence to mutated sequence across indels.

    Shared helper used by both dual mutation mode and single mutation mode.
    """
    translated_snps: list[dict[str, Any]] = []
    ref_assembly = getattr(args, "reference_assembly", None)
    if not isinstance(ref_assembly, str):
        ref_assembly = config.get("reference_assembly", "hg38")
    left_len = len(config.get("constants", {}).get(ref_assembly, {}).get("left", ""))

    if orig_normal_seqs is None:
        orig_normal_seqs = [hr.sequence for hr in normal_results]
    if orig_mut_seqs is None:
        orig_mut_seqs = [hr.sequence for hr in mutated_results]

    for hap_idx, snps in applied_snp_info_normal.items():
        if not snps or hap_idx >= len(mutated_results):
            continue
        normal_seq = orig_normal_seqs[hap_idx]
        mut_seq = orig_mut_seqs[hap_idx]
        delta = len(mut_seq) - len(normal_seq)

        norm_lengths = getattr(normal_results[hap_idx], "_repeat_lengths", None) or [
            len(config["repeats"][ru.symbol]) for ru in normal_results[hap_idx].chain
        ]
        mut_lengths = getattr(mutated_results[hap_idx], "_repeat_lengths", None)
        if mut_lengths is not None and len(mut_lengths) == len(mutated_results[hap_idx].chain):
            mut_lengths = list(mut_lengths)
        else:
            mut_lengths = [
                len(config["repeats"][ru.symbol]) for ru in mutated_results[hap_idx].chain
            ]
            if mutated_units and (hap_idx + 1) in mutated_units:
                for rep_idx, unit_seq in mutated_units[hap_idx + 1]:
                    if 1 <= rep_idx <= len(mut_lengths):
                        mut_lengths[rep_idx - 1] = len(unit_seq)

            discrepancy = delta - (sum(mut_lengths) - sum(norm_lengths))
            if discrepancy != 0:
                mut_indices = [
                    i
                    for i, ru in enumerate(mutated_results[hap_idx].chain)
                    if getattr(ru, "mutated", False)
                ]
                targets = mut_indices
                if not targets:
                    hap_targets = [
                        t
                        for t in _format_targets(
                            mutation_positions or getattr(args, "mutation_targets", None)
                        )
                        if t[0] == hap_idx + 1
                    ]
                    targets = [t[1] - 1 for t in hap_targets if 0 <= t[1] - 1 < len(mut_lengths)]
                if targets:
                    _infer_target_repeat_lengths(
                        normal_seq, mut_seq, left_len, norm_lengths, mut_lengths, targets, delta
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
                            while diff_idx < min_len and norm_rep[diff_idx] == mut_rep[diff_idx]:
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

    return translated_snps


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
    normal_results=None,
    mutated_units=None,
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
            translated_snps = _translate_snps_to_mutated(
                args=args,
                config=config,
                normal_results=results,
                mutated_results=mutated_results or [],
                applied_snp_info_normal=applied_snp_info_normal,
                mutation_positions=mutation_positions,
                orig_normal_seqs=orig_normal_seqs,
                orig_mut_seqs=orig_mut_seqs,
                mutated_units=mutated_units,
            )

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
                mutation_name = args.mutation_name or structure_mutation_info["name"]
                if mutation_positions or args.mutation_targets:
                    mutation_targets = _format_targets(mutation_positions or args.mutation_targets)
                elif structure_mutation_info and (
                    not args.mutation_name or args.mutation_name == structure_mutation_info["name"]
                ):
                    mutation_targets = structure_mutation_info["targets"]
                else:
                    mutation_targets = []

                for hap_idx, _rep_idx in mutation_targets:
                    if 1 <= hap_idx <= len(results):
                        haplotype_comments[hap_idx - 1] = (
                            f"Mutation Applied: {mutation_name} (Targets: {mutation_targets})"
                        )

            # Apply SNPs
            if getattr(args, "snp_input_file", None) and normal_results is not None:
                orig_normal_seqs = [hr.sequence for hr in normal_results]
                _, applied_snp_info_normal = integrate_snps_unified(
                    args, config, list(normal_results), rng=rng
                )
                translated_snps = _translate_snps_to_mutated(
                    args=args,
                    config=config,
                    normal_results=normal_results,
                    mutated_results=results,
                    applied_snp_info_normal=applied_snp_info_normal,
                    mutation_positions=mutation_positions,
                    orig_normal_seqs=orig_normal_seqs,
                    mutated_units=mutated_units,
                )
                results, applied_snp_info = integrate_snps_unified(
                    args,
                    config,
                    results,
                    skip_reference_check=True,
                    rng=rng,
                    predefined_snps=translated_snps,
                )
            else:
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
                if structure_mutation_info or args.mutation_name:
                    mutation_name = args.mutation_name or structure_mutation_info["name"]
                    if mutation_positions or args.mutation_targets:
                        target_tuples = _format_targets(mutation_positions or args.mutation_targets)
                    elif structure_mutation_info and (
                        not args.mutation_name
                        or args.mutation_name == structure_mutation_info["name"]
                    ):
                        target_tuples = structure_mutation_info["targets"]
                    else:
                        target_tuples = []
                    if target_tuples:
                        struct_fh.write(
                            f"# Mutation Applied: {mutation_name} (Targets: {target_tuples})\n"
                        )
                    else:
                        struct_fh.write(f"# Mutation Applied: {mutation_name} (Target: random)\n")

                for i, hr in enumerate(results, start=1):
                    chain_str = "-".join(str(ru) for ru in hr.chain)
                    struct_fh.write(f"haplotype_{i}\t{chain_str}\n")

            logging.info("Structure file written to %s", struct_out)
        except Exception as e:
            raise FileOperationError(f"Writing structure file failed: {e}") from e
