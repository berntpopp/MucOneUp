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
import time  # NEW: imported for simulation statistics timing

from .mutate import apply_mutations
from .read_simulation import simulate_reads as simulate_reads_pipeline
from .simulate import simulate_diploid, simulate_from_chains
from .translate import run_orf_finder_in_memory
from .fasta_writer import write_fasta  # Helper for writing FASTA files
from .version import __version__  # Import version from the single source
from .simulation_statistics import (
    generate_simulation_statistics,
    write_statistics_report,
)  # NEW: import simulation statistics module
from .io import parse_vntr_structure_file  # NEW: import structure file parser


def build_parser():
    """Build and return the argument parser."""
    parser = argparse.ArgumentParser(
        prog="MucOneUp", description="Tool to simulate MUC1 VNTR diploid references."
    )
    # Add a version flag.
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s " + __version__,
        help="Print the version number and exit.",
    )
    parser.add_argument(
        "--config", required=True, help="Path to JSON configuration file."
    )
    # Base name and output directory for all output files.
    parser.add_argument(
        "--out-base",
        default="muc1_simulated",
        help=(
            "Base name for all output files. All outputs (simulation FASTA, VNTR structure, ORF FASTA, "
            "read simulation outputs) will be named based on this base name. Default is 'muc1_simulated'."
        ),
    )
    parser.add_argument(
        "--out-dir",
        default=".",
        help="Output folder where all files will be written. Default is the current directory.",
    )
    parser.add_argument(
        "--num-haplotypes",
        default=2,
        type=int,
        help="Number of haplotypes to simulate. Typically 2 for diploid.",
    )

    # Create mutually exclusive group for simulation modes
    simulation_mode_group = parser.add_mutually_exclusive_group()

    # Move fixed-lengths to the exclusive group
    simulation_mode_group.add_argument(
        "--fixed-lengths",
        nargs="+",
        type=str,
        default=None,
        help=(
            "One or more fixed lengths (in VNTR repeats). If one value is provided, it applies to all haplotypes. "
            "If multiple, the number must match --num-haplotypes. Values may be a single integer (e.g. '60') or a range (e.g. '20-40')."
        ),
    )

    # New argument for providing exact VNTR structure from a file
    simulation_mode_group.add_argument(
        "--input-structure",
        type=str,
        default=None,
        help=(
            "Path to a file containing specific VNTR repeat compositions to use. "
            "The file should contain one line per haplotype in the format: haplotype_N<TAB>1-2-3-4-5-...-9."
        ),
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible simulations.",
    )
    # New flag to control series simulation with an optional step size.
    simulation_mode_group.add_argument(
        "--simulate-series",
        nargs="?",
        const=1,
        type=int,
        help=(
            "If provided, generate a simulation for each value in the fixed-length range(s). "
            "Optionally supply a step size (e.g. '--simulate-series 5' to use a step of 5). "
            "If omitted, the step defaults to 1 (i.e. every value). Without this flag, a fixed-length range results in a single simulation iteration with a random pick from each range."
        ),
    )
    # Mutation arguments.
    parser.add_argument(
        "--mutation-name",
        default=None,
        help=(
            "Name of the mutation to apply. To run dual simulations (normal and mutated), provide a comma-separated pair "
            "(e.g. 'normal,dupC')."
        ),
    )
    parser.add_argument(
        "--mutation-targets",
        nargs="+",
        default=None,
        help=(
            "One or more 'haplotype_index,repeat_index' pairs (1-based). E.g., '1,5 2,7'. "
            "If provided, each pair indicates which haplotype and repeat to mutate. If omitted, the mutation is applied at a random allowed repeat."
        ),
    )
    # Flags for additional outputs.
    parser.add_argument(
        "--output-structure",
        action="store_true",
        help="If provided, output a VNTR structure file using the normalized naming scheme.",
    )
    parser.add_argument(
        "--output-orfs",
        action="store_true",
        help=(
            "If provided, run ORF prediction and output an ORF FASTA file using the normalized naming scheme. "
            "Additionally, the resulting ORF file will be scanned for toxic protein sequence features and "
            "a statistics file will be generated."
        ),
    )
    # ORF parameters.
    parser.add_argument(
        "--orf-min-aa",
        type=int,
        default=100,
        help="Minimum peptide length (in amino acids) to report (default=100).",
    )
    parser.add_argument(
        "--orf-aa-prefix",
        nargs="?",
        const="MTSSV",
        default=None,
        help=(
            "Filter resulting peptides to only those beginning with this prefix. "
            "If used without a value, defaults to 'MTSSV'. If omitted, no prefix filter is applied."
        ),
    )
    # Read simulation flag.
    parser.add_argument(
        "--simulate-reads",
        action="store_true",
        help=(
            "If provided, run the read simulation pipeline on the simulated FASTA. "
            "This pipeline produces an aligned and indexed BAM and gzipped paired FASTQ files."
        ),
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL", "NONE"],
        help="Set logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL, NONE). Default is INFO.",
    )
    # Reference assembly parameter to override the config file setting
    parser.add_argument(
        "--reference-assembly",
        choices=["hg19", "hg38"],
        help="Specify reference assembly (hg19 or hg38). Overrides the setting in config file.",
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
            level=level, format="%(asctime)s - %(levelname)s - %(message)s"
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


def numbered_filename(
    out_dir: str, out_base: str, iteration: int, file_type: str, variant: str = ""
) -> str:
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

        # Override reference_assembly if specified in command line
        if args.reference_assembly:
            current_assembly = config.get("reference_assembly", "hg38")
            config["reference_assembly"] = args.reference_assembly
            logging.info(
                "Reference assembly overridden by command line: %s -> %s",
                current_assembly,
                args.reference_assembly,
            )
    except (FileNotFoundError, json.JSONDecodeError) as e:
        logging.error("Could not load config: %s", e)
        sys.exit(1)

    # Normalize output folder and base name.
    out_dir = args.out_dir
    out_base = args.out_base
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Process --input-structure if provided
    predefined_chains = None
    structure_mutation_info = None
    if args.input_structure is not None:
        try:
            logging.info("Using input structure file: %s", args.input_structure)
            predefined_chains, structure_mutation_info = parse_vntr_structure_file(
                args.input_structure, config
            )
            num_chains = len(predefined_chains)
            logging.info("Loaded %d haplotype chains from structure file", num_chains)

            if structure_mutation_info:
                logging.info(
                    "Found mutation information in structure file: %s",
                    structure_mutation_info["name"],
                )

            # Check if the number of chains matches the requested number of haplotypes
            if num_chains != args.num_haplotypes:
                logging.warning(
                    "--num-haplotypes=%d specified but structure file has %d chains. "
                    "Using the number of chains from the structure file.",
                    args.num_haplotypes,
                    num_chains,
                )

            # If --simulate-series is provided, warn that it's ignored
            if args.simulate_series is not None:
                logging.warning(
                    "--simulate-series is ignored when using --input-structure"
                )
        except Exception as e:
            logging.error("Error parsing structure file: %s", e)
            sys.exit(1)

    # Process --fixed-lengths if --input-structure is not provided
    if args.input_structure is None and args.fixed_lengths is not None:
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
            logging.info(
                "Series mode enabled with step %d: %d simulation iterations generated from fixed-length ranges.",
                step,
                len(simulation_configs),
            )
        else:
            simulation_configs = [[random.choice(lst) for lst in fixed_matrix]]
            logging.info(
                "Single simulation iteration generated using a random choice from each fixed-length range."
            )
    elif predefined_chains is not None:
        # For structure files, we create a single simulation configuration
        simulation_configs = ["from_structure"]
        logging.info("Using predefined VNTR chains from structure file.")
    else:
        simulation_configs = [None]  # Use random lengths if not provided.

    # Process mutation information - prioritize structure file over CLI args
    dual_mutation_mode = False
    mutation_pair = None

    # Convert structure file mutation info to CLI-style mutation info if present
    if structure_mutation_info:
        # Override the mutation_name with structure file data
        args.mutation_name = structure_mutation_info["name"]

        # Convert tuple format [(1, 25)] to string format ["1,25"] for CLI compatibility
        target_tuples = structure_mutation_info["targets"]
        args.mutation_targets = [
            f"{hap_idx},{rep_idx}" for hap_idx, rep_idx in target_tuples
        ]

        logging.info(
            "Using mutation information from structure file: %s at targets %s",
            args.mutation_name,
            args.mutation_targets,
        )
    # Process standard mutation args
    if args.mutation_name:
        if "," in args.mutation_name:
            mutation_pair = [s.strip() for s in args.mutation_name.split(",")]
            if mutation_pair[0].lower() != "normal":
                logging.error(
                    "In dual simulation mode, the first mutation-name must be 'normal'."
                )
                sys.exit(1)
            dual_mutation_mode = True
            logging.info("Using mutations from command line: %s", mutation_pair)
        else:
            logging.info("Using mutation from command line: %s", args.mutation_name)
    sim_index = 1
    for fixed_conf in simulation_configs:
        iteration_start = time.time()  # NEW: record iteration start time
        try:
            if fixed_conf == "from_structure":
                # Use predefined chains from structure file without applying mutations here
                # Mutations will be applied later using the standard pipeline
                results = simulate_from_chains(
                    predefined_chains=predefined_chains, config=config
                )
            else:
                # Use the standard simulation function
                results = simulate_diploid(
                    config=config,
                    num_haplotypes=args.num_haplotypes,
                    fixed_lengths=fixed_conf,
                    seed=args.seed,
                )
            logging.info(
                "Haplotype simulation completed successfully for iteration %d.",
                sim_index,
            )
        except Exception as e:
            logging.error("Simulation failed: %s", e)
            sys.exit(1)

        # Initialize mutated_units as None.
        mutated_units = None

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
                    if (
                        "mutations" not in config
                        or mutation_pair[1] not in config["mutations"]
                    ):
                        logging.error(
                            "Mutation '%s' not found in config['mutations'].",
                            mutation_pair[1],
                        )
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
                        logging.error(
                            "No repeats match 'allowed_repeats' for mutation '%s'",
                            mutation_pair[1],
                        )
                        sys.exit(1)
                    mutation_positions.append(random.choice(possible_targets))
                try:
                    mutated_results, mutated_units = apply_mutations(
                        config=config,
                        results=[(seq, chain.copy()) for seq, chain in results],
                        mutation_name=mutation_pair[1],
                        targets=mutation_positions,
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
                            # Handle either string format "1,25" or tuple format (1, 25)
                            if isinstance(t, str):
                                hap_str, rep_str = t.split(",")
                                hap_i = int(hap_str)
                                rep_i = int(rep_str)
                            elif isinstance(t, tuple) and len(t) == 2:
                                hap_i, rep_i = t
                            else:
                                raise ValueError(f"Unexpected target format: {t}")

                            mutation_positions.append((hap_i, rep_i))
                        except Exception as e:
                            logging.error(
                                "Invalid --mutation-targets format: '%s' (%s)", t, e
                            )
                            sys.exit(1)
                    try:
                        results, mutated_units = apply_mutations(
                            config=config,
                            results=results,
                            mutation_name=args.mutation_name,
                            targets=mutation_positions,
                        )
                        logging.info("Mutation applied at specified targets.")
                    except Exception as e:
                        logging.error("Mutation failed: %s", e)
                        sys.exit(1)
                else:
                    try:
                        if (
                            "mutations" not in config
                            or args.mutation_name not in config["mutations"]
                        ):
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
                        results, mutated_units = apply_mutations(
                            config=config,
                            results=results,
                            mutation_name=args.mutation_name,
                            targets=[rand_target],
                        )
                        logging.info(
                            "Mutation applied at random target: %s", str(rand_target)
                        )
                    except Exception as e:
                        logging.error("Random-target mutation failed: %s", e)
                        sys.exit(1)

        # Write FASTA outputs.
        if dual_mutation_mode:
            normal_out = numbered_filename(
                out_dir, out_base, sim_index, "simulated.fa", variant="normal"
            )
            mut_out = numbered_filename(
                out_dir, out_base, sim_index, "simulated.fa", variant="mut"
            )
            try:
                # Add no mutation information for normal results
                write_fasta(
                    [seq for seq, chain in normal_results],
                    normal_out,
                    comment="Normal sequence (no mutations applied)",
                )

                # Add mutation information for mutated results
                mutation_comment = f"Mutation Applied: {mutation_pair[1]} (Targets: {mutation_positions})"
                write_fasta(
                    [seq for seq, chain in mutated_results],
                    mut_out,
                    comment=mutation_comment,
                )

                logging.info(
                    "Dual FASTA outputs written: %s and %s", normal_out, mut_out
                )
            except Exception as e:
                logging.error("Writing FASTA failed: %s", e)
                sys.exit(1)
        else:
            out_file = numbered_filename(out_dir, out_base, sim_index, "simulated.fa")
            try:
                # Create haplotype-specific comments based on mutation information
                # Initialize with no comments
                haplotype_comments = [None] * len(results)

                # If we have mutation information, apply it only to the targeted haplotypes
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
                        # Process CLI mutation targets
                        for t in args.mutation_targets:
                            if isinstance(t, str):
                                hap_str, rep_str = t.split(",")
                                mutation_targets.append((int(hap_str), int(rep_str)))
                            elif isinstance(t, tuple) and len(t) == 2:
                                mutation_targets.append(t)

                    # Apply comments only to targeted haplotypes
                    for hap_idx, rep_idx in mutation_targets:
                        if 1 <= hap_idx <= len(results):  # 1-indexed haplotype numbers
                            haplotype_comments[hap_idx - 1] = (
                                f"Mutation Applied: {mutation_name} (Targets: {mutation_targets})"
                            )

                write_fasta(
                    [seq for seq, chain in results],
                    out_file,
                    prefix="haplotype",
                    comments=haplotype_comments,
                )
                logging.info("FASTA output written to %s", out_file)
            except Exception as e:
                logging.error("Writing FASTA failed: %s", e)
                sys.exit(1)

        # Write mutated VNTR unit FASTA if a mutation was applied.
        if args.mutation_name and mutated_units is not None:
            variant_suffix = "mut" if dual_mutation_mode else ""
            mutated_unit_out = numbered_filename(
                out_dir, out_base, sim_index, "mutated_unit.fa", variant=variant_suffix
            )
            try:
                with open(mutated_unit_out, "w") as muf:
                    for hap_idx, muts in mutated_units.items():
                        for rep_idx, unit_seq in muts:
                            muf.write(
                                f">haplotype_{hap_idx}_repeat_{rep_idx}\n{unit_seq}\n"
                            )
                logging.info(
                    "Mutated VNTR unit FASTA output written: %s", mutated_unit_out
                )
            except Exception as e:
                logging.error("Writing mutated VNTR unit FASTA failed: %s", e)
                sys.exit(1)

        # Write VNTR structure file.
        if args.output_structure:
            if dual_mutation_mode:
                normal_struct_out = numbered_filename(
                    out_dir, out_base, sim_index, "vntr_structure.txt", variant="normal"
                )
                mut_struct_out = numbered_filename(
                    out_dir, out_base, sim_index, "vntr_structure.txt", variant="mut"
                )
                try:
                    with open(normal_struct_out, "w") as nf:
                        # Add comment indicating this is a normal (non-mutated) structure
                        nf.write("# Normal sequence (no mutations applied)\n")
                        for i, (sequence, chain) in enumerate(normal_results, start=1):
                            chain_str = "-".join(chain)
                            nf.write(f"haplotype_{i}\t{chain_str}\n")
                    with open(mut_struct_out, "w") as mf:
                        # Add comment indicating the mutation that was applied
                        mutation_comment = f"# Mutation Applied: {mutation_pair[1]} (Targets: {mutation_positions})\n"
                        mf.write(mutation_comment)
                        for i, (sequence, chain) in enumerate(mutated_results, start=1):
                            chain_str = "-".join(chain)
                            mf.write(f"haplotype_{i}\t{chain_str}\n")
                    logging.info(
                        "Structure files written: %s and %s",
                        normal_struct_out,
                        mut_struct_out,
                    )
                except Exception as e:
                    logging.error("Writing structure file failed: %s", e)
                    sys.exit(1)
            else:
                struct_out = numbered_filename(
                    out_dir, out_base, sim_index, "vntr_structure.txt"
                )
                try:
                    with open(struct_out, "w") as struct_fh:
                        # Add mutation information - prioritize info from structure file
                        if structure_mutation_info:
                            # Preserve the original mutation information from input structure
                            struct_fh.write(
                                f"# Mutation Applied: {structure_mutation_info['name']} (Targets: {structure_mutation_info['targets']})\n"
                            )
                        elif args.mutation_name:
                            # Fall back to command-line mutation info
                            if args.mutation_targets:
                                struct_fh.write(
                                    f"# Mutation Applied: {args.mutation_name} (Targets: {args.mutation_targets})\n"
                                )
                            else:
                                struct_fh.write(
                                    f"# Mutation Applied: {args.mutation_name} (Target: random)\n"
                                )

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
                normal_orf_out = numbered_filename(
                    out_dir, out_base, sim_index, "orfs.fa", variant="normal"
                )
                mut_orf_out = numbered_filename(
                    out_dir, out_base, sim_index, "orfs.fa", variant="mut"
                )
                try:
                    run_orf_finder_in_memory(
                        normal_results,
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
                    logging.error("ORF finding failed: %s", e)
                    sys.exit(1)
            else:
                orf_out = numbered_filename(out_dir, out_base, sim_index, "orfs.fa")
                try:
                    run_orf_finder_in_memory(
                        results,
                        output_pep=orf_out,
                        orf_min_aa=args.orf_min_aa,
                        required_prefix=args.orf_aa_prefix,
                    )
                    logging.info(
                        "ORF finding completed; peptide FASTA written to %s", orf_out
                    )
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
            left_const_val = config.get("constants", {}).get("left")
            right_const_val = config.get("constants", {}).get("right")
            if dual_mutation_mode:
                normal_stats = scan_orf_fasta(
                    normal_orf_out,
                    left_const=left_const_val,
                    right_const=right_const_val,
                )
                mut_stats = scan_orf_fasta(
                    mut_orf_out, left_const=left_const_val, right_const=right_const_val
                )
                stats_file_normal = numbered_filename(
                    out_dir, out_base, sim_index, "orf_stats.txt", variant="normal"
                )
                stats_file_mut = numbered_filename(
                    out_dir, out_base, sim_index, "orf_stats.txt", variant="mut"
                )
                with open(stats_file_normal, "w") as nf:
                    json.dump(normal_stats, nf, indent=4)
                with open(stats_file_mut, "w") as mf:
                    json.dump(mut_stats, mf, indent=4)
                logging.info(
                    "Toxic protein detection stats written: %s and %s",
                    stats_file_normal,
                    stats_file_mut,
                )
            else:
                stats_file = numbered_filename(
                    out_dir, out_base, sim_index, "orf_stats.txt"
                )
                stats = scan_orf_fasta(
                    orf_out, left_const=left_const_val, right_const=right_const_val
                )
                with open(stats_file, "w") as sf:
                    json.dump(stats, sf, indent=4)
                logging.info("Toxic protein detection stats written: %s", stats_file)
            # ----------------------------------------------------------------

        # Run read simulation pipeline if requested.
        if args.simulate_reads:
            if dual_mutation_mode:
                normal_fa_for_reads = numbered_filename(
                    out_dir, out_base, sim_index, "simulated.fa", variant="normal"
                )
                mut_fa_for_reads = numbered_filename(
                    out_dir, out_base, sim_index, "simulated.fa", variant="mut"
                )
                try:
                    logging.info(
                        "Starting read simulation pipeline for iteration %d (normal variant).",
                        sim_index,
                    )
                    simulate_reads_pipeline(config, normal_fa_for_reads)
                    logging.info(
                        "Read simulation pipeline completed for iteration %d (normal variant).",
                        sim_index,
                    )
                except Exception as e:
                    logging.error(
                        "Read simulation pipeline (normal variant) failed: %s", e
                    )
                    sys.exit(1)
                try:
                    logging.info(
                        "Starting read simulation pipeline for iteration %d (mutated variant).",
                        sim_index,
                    )
                    simulate_reads_pipeline(config, mut_fa_for_reads)
                    logging.info(
                        "Read simulation pipeline completed for iteration %d (mutated variant).",
                        sim_index,
                    )
                except Exception as e:
                    logging.error(
                        "Read simulation pipeline (mutated variant) failed: %s", e
                    )
                    sys.exit(1)
            else:
                sim_fa_for_reads = numbered_filename(
                    out_dir, out_base, sim_index, "simulated.fa"
                )
                try:
                    logging.info(
                        "Starting read simulation pipeline for iteration %d.", sim_index
                    )
                    simulate_reads_pipeline(config, sim_fa_for_reads)
                    logging.info(
                        "Read simulation pipeline completed for iteration %d.",
                        sim_index,
                    )
                except Exception as e:
                    logging.error("Read simulation pipeline failed: %s", e)
                    sys.exit(1)

        # ----------------------------------------------------------------
        # NEW: Generate simulation statistics report for this iteration.
        iteration_end = time.time()
        vntr_coverage_stats = (
            {}
        )  # This can be filled with VNTR coverage data if available
        if dual_mutation_mode:
            normal_stats_report = generate_simulation_statistics(
                start_time=iteration_start,
                end_time=iteration_end,
                simulation_results=normal_results,
                config=config,
                mutation_info={"mutation_name": "normal"},
                vntr_coverage=vntr_coverage_stats,
            )
            mutated_stats_report = generate_simulation_statistics(
                start_time=iteration_start,
                end_time=iteration_end,
                simulation_results=mutated_results,
                config=config,
                mutation_info={"mutation_name": mutation_pair[1]},
                vntr_coverage=vntr_coverage_stats,
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
            simulation_results_for_stats = results
            mutation_info = {}
            if args.mutation_name:
                mutation_info = {
                    "mutation_name": args.mutation_name,
                    "mutation_targets": args.mutation_targets,
                }
            stats_report = generate_simulation_statistics(
                start_time=iteration_start,
                end_time=iteration_end,
                simulation_results=simulation_results_for_stats,
                config=config,
                mutation_info=mutation_info,
                vntr_coverage=vntr_coverage_stats,
            )
            stats_output_file = numbered_filename(
                out_dir, out_base, sim_index, "simulation_stats.json"
            )
            write_statistics_report(stats_report, stats_output_file)
        # ----------------------------------------------------------------

        sim_index += 1


if __name__ == "__main__":
    main()
