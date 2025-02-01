# muc_one_up/cli.py
import argparse
import json
import logging
import random
import sys
import os
import itertools

from .mutate import apply_mutations
from .read_simulation import simulate_reads as simulate_reads_pipeline
from .simulate import simulate_diploid
from .translate import run_orf_finder_in_memory
from .fasta_writer import write_fasta  # Assumes this helper writes FASTA files


def build_parser():
    """Build and return the argument parser."""
    parser = argparse.ArgumentParser(
        prog="MucOneUp",
        description="Tool to simulate MUC1 VNTR diploid references."
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to JSON configuration file."
    )
    # New output flags replacing previous individual output flags.
    parser.add_argument(
        "--out-base",
        default="muc1_simulated",
        help=("Base name for all output files. All outputs (simulation FASTA, VNTR structure, ORF FASTA, "
              "and read simulation outputs) will be named based on this base name. Default is 'muc1_simulated'.")
    )
    parser.add_argument(
        "--out-dir",
        default=".",
        help=("Output folder where all files will be written. Default is the current directory.")
    )
    parser.add_argument(
        "--num-haplotypes",
        default=2,
        type=int,
        help="Number of haplotypes to simulate. Typically 2 for diploid."
    )
    # Allow fixed-lengths as strings so that ranges (e.g. "20-40") are allowed.
    parser.add_argument(
        "--fixed-lengths",
        nargs="+",
        type=str,
        default=None,
        help=("One or more fixed lengths (in VNTR repeats). If one value is provided, it will apply to all haplotypes. "
              "If multiple, must match --num-haplotypes. Values may be a single integer (e.g. '60') or a range (e.g. '20-40').")
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible simulations."
    )
    # Mutation arguments remain largely unchanged.
    parser.add_argument(
        "--mutation-name",
        default=None,
        help=("Name of the mutation to apply. To run dual simulations (normal and mutated), provide a comma-separated pair "
              "(e.g. 'normal,dupC').")
    )
    parser.add_argument(
        "--mutation-targets",
        nargs="+",
        default=None,
        help=("One or more 'haplotype_index,repeat_index' pairs, e.g. '1,5' '2,7'. If provided, each pair indicates the haplotype "
              "and repeat number at which to apply the named mutation once. If omitted, the mutation (if any) is applied once "
              "at a random allowed repeat.")
    )
    # New flags to control whether to output structure and ORF files.
    parser.add_argument(
        "--output-structure",
        action="store_true",
        help="If provided, output a VNTR structure file using the normalized naming scheme."
    )
    parser.add_argument(
        "--output-orfs",
        action="store_true",
        help="If provided, run ORF prediction and output an ORF FASTA file using the normalized naming scheme."
    )
    # New ORF parameters remain the same.
    parser.add_argument(
        "--orf-min-aa",
        type=int,
        default=100,
        help="Minimum peptide length (in amino acids) to report (default=100)."
    )
    parser.add_argument(
        "--orf-aa-prefix",
        nargs='?',
        const='MTSSV',
        default=None,
        help=("Filter resulting peptides to only those beginning with this prefix. "
              "If used without a value, defaults to 'MTSSV'. If omitted entirely, no prefix filter is applied.")
    )
    # Read simulation flag remains unchanged.
    parser.add_argument(
        "--simulate-reads",
        action="store_true",
        help=("If provided, run the read simulation pipeline on the simulated FASTA. "
              "This pipeline will produce an aligned and indexed BAM and gzipped paired FASTQ files.")
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL", "NONE"],
        help=("Set logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL, NONE). Default is INFO.")
    )
    return parser


def configure_logging(level_str):
    """
    Configure root logging based on the provided level string.

    If level_str is 'NONE', disable logging.
    """
    if level_str.upper() == "NONE":
        logging.disable(logging.CRITICAL)
    else:
        level = getattr(logging, level_str.upper(), logging.INFO)
        root_logger = logging.getLogger()
        if root_logger.handlers:
            for handler in root_logger.handlers[:]:
                root_logger.removeHandler(handler)
        logging.basicConfig(
            level=level,
            format="%(asctime)s - %(levelname)s - %(message)s"
        )
        logging.info("Logging configured at level: %s", level_str.upper())


def parse_fixed_lengths(fixed_lengths_args, num_haplotypes):
    """
    Parse the fixed-length values provided as strings. Each value may be a single integer
    or a range in the format "start-end". Returns a list of lists, one per haplotype.
    """
    if len(fixed_lengths_args) not in (1, num_haplotypes):
        logging.error("--fixed-lengths must have either 1 value or N=--num-haplotypes")
        sys.exit(1)
    result = []
    for val in fixed_lengths_args:
        if "-" in val:
            try:
                start_str, end_str = val.split("-")
                start = int(start_str)
                end = int(end_str)
                result.append(list(range(start, end + 1)))
            except ValueError:
                logging.error("Invalid fixed-length range format: '%s'", val)
                sys.exit(1)
        else:
            try:
                fixed = int(val)
                result.append([fixed])
            except ValueError:
                logging.error("Invalid fixed-length value: '%s'", val)
                sys.exit(1)
    if len(result) == 1 and num_haplotypes > 1:
        result = result * num_haplotypes
    return result


def build_cartesian_fixed_length_configs(fixed_matrix):
    """
    Given a list of lists (each corresponding to one haplotype's possible fixed lengths),
    return the Cartesian product as a list of fixed-length configurations (each a list).
    """
    return [list(prod) for prod in itertools.product(*fixed_matrix)]


def numbered_filename(out_dir: str, out_base: str, iteration: int, file_type: str, variant: str = "") -> str:
    """
    Build a filename by combining the output directory, base name, iteration number,
    variant suffix, and a file type suffix.

    :param out_dir: The output folder.
    :param out_base: The base name for outputs.
    :param iteration: Simulation iteration number.
    :param file_type: A string describing the file type (e.g., 'simulated.fa', 'vntr_structure.txt').
    :param variant: Optional variant suffix (e.g., 'mut' or 'normal').
    :return: A string filename.
    """
    iter_str = f".{iteration:03d}"
    variant_str = f".{variant}" if variant else ""
    return os.path.join(out_dir, f"{out_base}{iter_str}{variant_str}.{file_type}")


def main():
    """Main entry point for the CLI."""
    parser = build_parser()
    args = parser.parse_args()

    # Configure logging
    configure_logging(args.log_level)

    # Load config
    try:
        with open(args.config) as fh:
            config = json.load(fh)
        logging.info("Configuration loaded from %s", args.config)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        logging.error("Could not load config: %s", e)
        sys.exit(1)

    # Normalize output folder and base name.
    out_dir = args.out_dir
    out_base = args.out_base
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Process --fixed-lengths logic
    if args.fixed_lengths is not None:
        fixed_matrix = parse_fixed_lengths(args.fixed_lengths, args.num_haplotypes)
        sequence_mode = any(len(lst) > 1 for lst in fixed_matrix)
        if sequence_mode:
            simulation_configs = build_cartesian_fixed_length_configs(fixed_matrix)
        else:
            simulation_configs = [[lst[0] for lst in fixed_matrix]]
    else:
        simulation_configs = [None]  # Use random lengths

    # Process mutation-name: check for comma-separated pair.
    dual_mutation_mode = False
    mutation_pair = None
    if args.mutation_name:
        if ',' in args.mutation_name:
            mutation_pair = [s.strip() for s in args.mutation_name.split(',')]
            if mutation_pair[0].lower() != "normal":
                logging.error("In dual simulation mode, the first mutation-name must be 'normal'.")
                sys.exit(1)
            dual_mutation_mode = True

    sim_index = 1
    for fixed_conf in simulation_configs:
        try:
            results = simulate_diploid(
                config=config,
                num_haplotypes=args.num_haplotypes,
                fixed_lengths=fixed_conf,
                seed=args.seed
            )
            logging.info("Haplotype simulation completed successfully for iteration %d.", sim_index)
        except Exception as e:
            logging.error("Simulation failed: %s", e)
            sys.exit(1)

        # MUTATION LOGIC
        if args.mutation_name:
            logging.info("Applying mutation: %s", args.mutation_name)
            if dual_mutation_mode:
                normal_results = results
                mutation_positions = None
                if args.mutation_targets:
                    mutation_positions = []
                    for t in args.mutation_targets:
                        try:
                            hap_str, rep_str = t.split(",")
                            hap_i = int(hap_str)
                            rep_i = int(rep_str)
                            mutation_positions.append((hap_i, rep_i))
                        except ValueError:
                            logging.error("Invalid --mutation-targets format: '%s'", t)
                            sys.exit(1)
                try:
                    mutated_results = apply_mutations(
                        config=config,
                        results=[(seq, chain.copy()) for seq, chain in results],
                        mutation_name=mutation_pair[1],
                        targets=mutation_positions if mutation_positions else None
                    )
                    logging.info("Dual mutation: mutation applied for mutated version.")
                except Exception as e:
                    logging.error("Mutation failed: %s", e)
                    sys.exit(1)
            else:
                if args.mutation_targets:
                    mutation_positions = []
                    for t in args.mutation_targets:
                        try:
                            hap_str, rep_str = t.split(",")
                            hap_i = int(hap_str)
                            rep_i = int(rep_str)
                            mutation_positions.append((hap_i, rep_i))
                        except ValueError:
                            logging.error("Invalid --mutation-targets format: '%s'", t)
                            sys.exit(1)
                    try:
                        results = apply_mutations(
                            config=config,
                            results=results,
                            mutation_name=args.mutation_name,
                            targets=mutation_positions
                        )
                        logging.info("Mutation applied at specified targets.")
                    except Exception as e:
                        logging.error("Mutation failed: %s", e)
                        sys.exit(1)
                else:
                    try:
                        if ("mutations" not in config or
                                args.mutation_name not in config["mutations"]):
                            raise ValueError(
                                f"Mutation '{args.mutation_name}' not in config['mutations']"
                            )
                        mut_def = config["mutations"][args.mutation_name]
                        allowed_repeats = set(mut_def["allowed_repeats"])
                        possible_targets = []
                        for hap_idx, (seq, chain) in enumerate(results, start=1):
                            for rep_idx, sym in enumerate(chain, start=1):
                                pure_sym = sym.replace("m", "")
                                if pure_sym in allowed_repeats:
                                    possible_targets.append((hap_idx, rep_idx))
                        if not possible_targets:
                            raise ValueError(
                                f"No repeats match 'allowed_repeats' for '{args.mutation_name}'"
                            )
                        rand_target = random.choice(possible_targets)
                        results = apply_mutations(
                            config=config,
                            results=results,
                            mutation_name=args.mutation_name,
                            targets=[rand_target]
                        )
                        logging.info("Mutation applied at random target: %s", str(rand_target))
                    except Exception as e:
                        logging.error("Random-target mutation failed: %s", e)
                        sys.exit(1)

        # Write final FASTA outputs using normalized naming.
        if dual_mutation_mode:
            normal_out = numbered_filename(out_dir, out_base, sim_index, "simulated.fa", variant="normal")
            mut_out = numbered_filename(out_dir, out_base, sim_index, "simulated.fa", variant="mut")
            try:
                # Extract only the assembled sequence (first element of each tuple)
                write_fasta([seq for seq, chain in normal_results], normal_out)
                write_fasta([seq for seq, chain in mutated_results], mut_out)
                logging.info("Dual FASTA outputs written: %s and %s", normal_out, mut_out)
            except Exception as e:
                logging.error("Writing FASTA failed: %s", e)
                sys.exit(1)
        else:
            out_file = numbered_filename(out_dir, out_base, sim_index, "simulated.fa")
            try:
                write_fasta([seq for seq, chain in results], out_file)
                logging.info("FASTA output written to %s", out_file)
            except Exception as e:
                logging.error("Writing FASTA failed: %s", e)
                sys.exit(1)

        # Optionally write VNTR structure file if flag is set.
        if args.output_structure:
            struct_out = numbered_filename(out_dir, out_base, sim_index, "vntr_structure.txt")
            try:
                with open(struct_out, "w") as struct_fh:
                    for i, (sequence, chain) in enumerate(results, start=1):
                        chain_str = "-".join(chain)
                        struct_fh.write(f"haplotype_{i}\t{chain_str}\n")
                logging.info("Structure file written to %s", struct_out)
            except Exception as e:
                logging.error("Writing structure file failed: %s", e)
                sys.exit(1)

        # ORF logic: if flag is set, run ORF prediction.
        if args.output_orfs:
            orf_out = numbered_filename(out_dir, out_base, sim_index, "orfs.fa")
            try:
                run_orf_finder_in_memory(
                    results,
                    output_pep=orf_out,
                    orf_min_aa=args.orf_min_aa,
                    required_prefix=args.orf_aa_prefix
                )
                logging.info("ORF finding completed; peptide FASTA written to %s", orf_out)
            except Exception as e:
                logging.error("ORF finding failed: %s", e)
                sys.exit(1)

        # Run read simulation pipeline if requested.
        if args.simulate_reads:
            sim_fa_for_reads = numbered_filename(out_dir, out_base, sim_index, "simulated.fa")
            try:
                logging.info("Starting read simulation pipeline for iteration %d.", sim_index)
                simulate_reads_pipeline(config, sim_fa_for_reads)
                logging.info("Read simulation pipeline completed for iteration %d.", sim_index)
            except Exception as e:
                logging.error("Read simulation pipeline failed: %s", e)
                sys.exit(1)

        sim_index += 1


if __name__ == "__main__":
    main()
