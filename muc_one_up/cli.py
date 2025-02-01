#!/usr/bin/env python3
"""
cli.py

This module implements the command-line interface for MucOneUp.
It supports:
  - Generating multiple simulation iterations using a fixed-length range.
    If the --simulate-series flag is provided, a simulation is run for each value in the range(s);
    an optional step size can be provided (default is 1, meaning every value).
    Without the flag, a fixed-length range produces a single simulation iteration by choosing one random value per range.
  - Optionally producing dual outputs (normal and mutated) if a comma-separated mutation pair is provided.
  - Printing the version of the tool.
  - Running ORF prediction and, when activated via --output-orfs, scanning the resulting ORF FASTA
    for toxic protein sequence features and generating a statistics output.
"""

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
from .fasta_writer import write_fasta  # Helper for writing FASTA files
from .version import __version__  # Import version from the single source

# (Other functions remain unchanged.)

def build_parser():
    """Build and return the argument parser."""
    parser = argparse.ArgumentParser(
        prog="MucOneUp",
        description="Tool to simulate MUC1 VNTR diploid references."
    )
    # Add a version flag.
    parser.add_argument(
        "-v", "--version",
        action="version",
        version="%(prog)s " + __version__,
        help="Print the version number and exit."
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to JSON configuration file."
    )
    # Base name and output directory for all output files.
    parser.add_argument(
        "--out-base",
        default="muc1_simulated",
        help=("Base name for all output files. All outputs (simulation FASTA, VNTR structure, ORF FASTA, "
              "read simulation outputs) will be named based on this base name. Default is 'muc1_simulated'.")
    )
    parser.add_argument(
        "--out-dir",
        default=".",
        help="Output folder where all files will be written. Default is the current directory."
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
        help=("One or more fixed lengths (in VNTR repeats). If one value is provided, it applies to all haplotypes. "
              "If multiple, the number must match --num-haplotypes. Values may be a single integer (e.g. '60') or a range (e.g. '20-40').")
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible simulations."
    )
    # New flag to control series simulation with an optional step size.
    parser.add_argument(
        "--simulate-series",
        nargs="?",
        const=1,
        type=int,
        help=("If provided, generate a simulation for each value in the fixed-length range(s). "
              "Optionally supply a step size (e.g. '--simulate-series 5' to use a step of 5). "
              "If omitted, the step defaults to 1 (i.e. every value). Without this flag, a fixed-length range results in a single simulation iteration with a random pick from each range.")
    )
    # Mutation arguments.
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
        help=("One or more 'haplotype_index,repeat_index' pairs (1-based). E.g., '1,5 2,7'. "
              "If provided, each pair indicates which haplotype and repeat to mutate. If omitted, the mutation is applied at a random allowed repeat.")
    )
    # Flags for additional outputs.
    parser.add_argument(
        "--output-structure",
        action="store_true",
        help="If provided, output a VNTR structure file using the normalized naming scheme."
    )
    parser.add_argument(
        "--output-orfs",
        action="store_true",
        help=("If provided, run ORF prediction and output an ORF FASTA file using the normalized naming scheme. "
              "Additionally, the resulting ORF file will be scanned for toxic protein sequence features and "
              "a statistics file will be generated.")
    )
    # ORF parameters.
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
              "If used without a value, defaults to 'MTSSV'. If omitted, no prefix filter is applied.")
    )
    # Read simulation flag.
    parser.add_argument(
        "--simulate-reads",
        action="store_true",
        help=("If provided, run the read simulation pipeline on the simulated FASTA. "
              "This pipeline produces an aligned and indexed BAM and gzipped paired FASTQ files.")
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL", "NONE"],
        help="Set logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL, NONE). Default is INFO."
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
    Parse the fixed-length values provided as strings.
    Each value may be a single integer or a range in the format "start-end".
    Returns a list of lists, one per haplotype.
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
    """
    iter_str = f".{iteration:03d}"
    variant_str = f".{variant}" if variant else ""
    return os.path.join(out_dir, f"{out_base}{iter_str}{variant_str}.{file_type}")


def main():
    """Main entry point for the CLI."""
    parser = build_parser()
    args = parser.parse_args()

    # Configure logging.
    configure_logging(args.log_level)

    # Load configuration.
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

    # Process --fixed-lengths.
    if args.fixed_lengths is not None:
        fixed_matrix = parse_fixed_lengths(args.fixed_lengths, args.num_haplotypes)
        if args.simulate_series is not None:
            step = args.simulate_series
            new_fixed_matrix = []
            for lst in fixed_matrix:
                if len(lst) > 1:
                    start = lst[0]
                    end = lst[-1]
                    new_lst = list(range(start, end + 1, step))
                    if new_lst[-1] != end:
                        new_lst.append(end)
                    new_fixed_matrix.append(new_lst)
                else:
                    new_fixed_matrix.append(lst)
            fixed_matrix = new_fixed_matrix
            simulation_configs = build_cartesian_fixed_length_configs(fixed_matrix)
            logging.info("Series mode enabled with step %d: %d simulation iterations generated from fixed-length ranges.",
                         step, len(simulation_configs))
        else:
            simulation_configs = [[random.choice(lst) for lst in fixed_matrix]]
            logging.info("Single simulation iteration generated using a random choice from each fixed-length range.")
    else:
        simulation_configs = [None]  # Use random lengths if not provided.

    # Process --mutation-name to support dual simulation.
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
                mutation_positions = []
                if args.mutation_targets:
                    for t in args.mutation_targets:
                        try:
                            hap_str, rep_str = t.split(",")
                            hap_i = int(hap_str)
                            rep_i = int(rep_str)
                            mutation_positions.append((hap_i, rep_i))
                        except ValueError:
                            logging.error("Invalid --mutation-targets format: '%s'", t)
                            sys.exit(1)
                else:
                    if "mutations" not in config or mutation_pair[1] not in config["mutations"]:
                        logging.error("Mutation '%s' not found in config['mutations'].", mutation_pair[1])
                        sys.exit(1)
                    mut_def = config["mutations"][mutation_pair[1]]
                    allowed_repeats = set(mut_def["allowed_repeats"])
                    possible_targets = []
                    for hap_idx, (seq, chain) in enumerate(results, start=1):
                        for rep_idx, sym in enumerate(chain, start=1):
                            pure_sym = sym.replace("m", "")
                            if pure_sym in allowed_repeats:
                                possible_targets.append((hap_idx, rep_idx))
                    if not possible_targets:
                        logging.error("No repeats match 'allowed_repeats' for mutation '%s'", mutation_pair[1])
                        sys.exit(1)
                    mutation_positions.append(random.choice(possible_targets))
                try:
                    mutated_results = apply_mutations(
                        config=config,
                        results=[(seq, chain.copy()) for seq, chain in results],
                        mutation_name=mutation_pair[1],
                        targets=mutation_positions
                    )
                    logging.info("Dual mutation applied for mutated version.")
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

        # Write FASTA outputs.
        if dual_mutation_mode:
            normal_out = numbered_filename(out_dir, out_base, sim_index, "simulated.fa", variant="normal")
            mut_out = numbered_filename(out_dir, out_base, sim_index, "simulated.fa", variant="mut")
            try:
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

        # Write VNTR structure file.
        if args.output_structure:
            if dual_mutation_mode:
                normal_struct_out = numbered_filename(out_dir, out_base, sim_index, "vntr_structure.txt", variant="normal")
                mut_struct_out = numbered_filename(out_dir, out_base, sim_index, "vntr_structure.txt", variant="mut")
                try:
                    with open(normal_struct_out, "w") as nf:
                        for i, (sequence, chain) in enumerate(normal_results, start=1):
                            chain_str = "-".join(chain)
                            nf.write(f"haplotype_{i}\t{chain_str}\n")
                    with open(mut_struct_out, "w") as mf:
                        for i, (sequence, chain) in enumerate(mutated_results, start=1):
                            chain_str = "-".join(chain)
                            mf.write(f"haplotype_{i}\t{chain_str}\n")
                    logging.info("Structure files written: %s and %s", normal_struct_out, mut_struct_out)
                except Exception as e:
                    logging.error("Writing structure file failed: %s", e)
                    sys.exit(1)
            else:
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

        # ORF logic.
        if args.output_orfs:
            if dual_mutation_mode:
                normal_orf_out = numbered_filename(out_dir, out_base, sim_index, "orfs.fa", variant="normal")
                mut_orf_out = numbered_filename(out_dir, out_base, sim_index, "orfs.fa", variant="mut")
                try:
                    run_orf_finder_in_memory(
                        normal_results,
                        output_pep=normal_orf_out,
                        orf_min_aa=args.orf_min_aa,
                        required_prefix=args.orf_aa_prefix
                    )
                    run_orf_finder_in_memory(
                        mutated_results,
                        output_pep=mut_orf_out,
                        orf_min_aa=args.orf_min_aa,
                        required_prefix=args.orf_aa_prefix
                    )
                    logging.info("ORF finding completed; peptide FASTA outputs written: %s and %s", normal_orf_out, mut_orf_out)
                except Exception as e:
                    logging.error("ORF finding failed: %s", e)
                    sys.exit(1)
            else:
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

            # ----------------------------------------------------------------
            # Toxic Protein Detection Integration:
            # After ORF prediction, run the toxic protein detector on the generated ORF FASTA file(s)
            # and write a JSON statistics file with detection metrics for each ORF.
            try:
                from .toxic_protein_detector import scan_orf_fasta
            except ImportError as e:
                logging.error("Failed to import toxic_protein_detector module: %s", e)
                sys.exit(1)
            # Optionally, use constant flanks from the config if available.
            left_const_val = config.get("constants", {}).get("left")
            right_const_val = config.get("constants", {}).get("right")
            if dual_mutation_mode:
                normal_stats = scan_orf_fasta(normal_orf_out, left_const=left_const_val, right_const=right_const_val)
                mut_stats = scan_orf_fasta(mut_orf_out, left_const=left_const_val, right_const=right_const_val)
                stats_file_normal = numbered_filename(out_dir, out_base, sim_index, "orf_stats.txt", variant="normal")
                stats_file_mut = numbered_filename(out_dir, out_base, sim_index, "orf_stats.txt", variant="mut")
                with open(stats_file_normal, "w") as nf:
                    json.dump(normal_stats, nf, indent=4)
                with open(stats_file_mut, "w") as mf:
                    json.dump(mut_stats, mf, indent=4)
                logging.info("Toxic protein detection stats written: %s and %s", stats_file_normal, stats_file_mut)
            else:
                stats_file = numbered_filename(out_dir, out_base, sim_index, "orf_stats.txt")
                stats = scan_orf_fasta(orf_out, left_const=left_const_val, right_const=right_const_val)
                with open(stats_file, "w") as sf:
                    json.dump(stats, sf, indent=4)
                logging.info("Toxic protein detection stats written: %s", stats_file)
            # ----------------------------------------------------------------

        # Run read simulation pipeline if requested.
        if args.simulate_reads:
            if dual_mutation_mode:
                normal_fa_for_reads = numbered_filename(out_dir, out_base, sim_index, "simulated.fa", variant="normal")
                mut_fa_for_reads = numbered_filename(out_dir, out_base, sim_index, "simulated.fa", variant="mut")
                try:
                    logging.info("Starting read simulation pipeline for iteration %d (normal variant).", sim_index)
                    simulate_reads_pipeline(config, normal_fa_for_reads)
                    logging.info("Read simulation pipeline completed for iteration %d (normal variant).", sim_index)
                except Exception as e:
                    logging.error("Read simulation pipeline (normal variant) failed: %s", e)
                    sys.exit(1)
                try:
                    logging.info("Starting read simulation pipeline for iteration %d (mutated variant).", sim_index)
                    simulate_reads_pipeline(config, mut_fa_for_reads)
                    logging.info("Read simulation pipeline completed for iteration %d (mutated variant).", sim_index)
                except Exception as e:
                    logging.error("Read simulation pipeline (mutated variant) failed: %s", e)
                    sys.exit(1)
            else:
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
