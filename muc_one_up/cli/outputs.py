"""
Output writing functions for MucOneUp CLI.

Single Responsibility: Write simulation outputs (FASTA, structure files, mutated units).
"""

import logging
from pathlib import Path

from ..exceptions import FileOperationError
from ..fasta_writer import write_fasta
from .config import numbered_filename
from .snps import integrate_snps_unified


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
            # Apply SNPs to both variants
            results, applied_snp_info_normal = integrate_snps_unified(args, config, results)
            mutated_results, applied_snp_info_mut = integrate_snps_unified(
                args, config, mutated_results, skip_reference_check=True
            )

            write_fasta(
                [seq for seq, chain in results],
                normal_out,
                comment="Normal sequence (no mutations applied)",
            )
            mutation_comment = (
                f"Mutation Applied: {mutation_pair[1]} (Targets: {mutation_positions})"
            )
            write_fasta(
                [seq for seq, chain in mutated_results],
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
                mutation_targets = []

                if structure_mutation_info:
                    mutation_targets = structure_mutation_info["targets"]
                elif args.mutation_targets:
                    for t in args.mutation_targets:
                        if isinstance(t, str):
                            hap_str, rep_str = t.split(",")
                            mutation_targets.append((int(hap_str), int(rep_str)))
                        elif isinstance(t, tuple) and len(t) == 2:
                            mutation_targets.append(t)

                for hap_idx, _rep_idx in mutation_targets:
                    if 1 <= hap_idx <= len(results):
                        haplotype_comments[hap_idx - 1] = (
                            f"Mutation Applied: {mutation_name} (Targets: {mutation_targets})"
                        )

            # Apply SNPs
            results, applied_snp_info = integrate_snps_unified(args, config, results)

            write_fasta(
                [seq for seq, chain in results],
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
                for i, (_sequence, chain) in enumerate(results, start=1):
                    chain_str = "-".join(chain)
                    nf.write(f"haplotype_{i}\t{chain_str}\n")

            with Path(mut_struct_out).open("w") as mf:
                mutation_comment = (
                    f"# Mutation Applied: {mutation_pair[1]} (Targets: {mutation_positions})\n"
                )
                mf.write(mutation_comment)
                for i, (_sequence, chain) in enumerate(mutated_results, start=1):
                    chain_str = "-".join(chain)
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
                    if args.mutation_targets:
                        struct_fh.write(
                            f"# Mutation Applied: {args.mutation_name} (Targets: {args.mutation_targets})\n"
                        )
                    else:
                        struct_fh.write(
                            f"# Mutation Applied: {args.mutation_name} (Target: random)\n"
                        )

                for i, (_sequence, chain) in enumerate(results, start=1):
                    chain_str = "-".join(chain)
                    struct_fh.write(f"haplotype_{i}\t{chain_str}\n")

            logging.info("Structure file written to %s", struct_out)
        except Exception as e:
            raise FileOperationError(f"Writing structure file failed: {e}") from e
