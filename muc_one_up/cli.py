# muc_one_up/cli.py
import argparse
import json
import logging
import random
import sys

from .mutate import apply_mutations
from .read_simulation import simulate_reads as simulate_reads_pipeline
from .simulate import simulate_diploid
from .translate import run_orf_finder_in_memory


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
    parser.add_argument(
        "--output",
        default="muc1_simulated.fa",
        help="Output FASTA filename."
    )
    parser.add_argument(
        "--structure-output",
        default=None,
        help=("Optional: Output text file describing the VNTR structure of "
              "each haplotype.")
    )
    parser.add_argument(
        "--num-haplotypes",
        default=2,
        type=int,
        help="Number of haplotypes to simulate. Typically 2 for diploid."
    )
    # multiple fixed lengths
    parser.add_argument(
        "--fixed-lengths",
        nargs="+",
        type=int,
        default=None,
        help=("One or more fixed lengths (in VNTR repeats). "
              "If one value is provided, it will apply to all haplotypes. "
              "If multiple, must match --num-haplotypes.")
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible simulations."
    )
    # mutation arguments
    parser.add_argument(
        "--mutation-name",
        default=None,
        help=("Name of the mutation to apply (must match an entry in "
              "config['mutations']).")
    )
    parser.add_argument(
        "--mutation-targets",
        nargs="+",
        default=None,
        help=("One or more 'haplotype_index,repeat_index' pairs, e.g. '1,5' '2,7'. "
              "If provided, each pair indicates the haplotype and repeat number "
              "at which to apply the named mutation once. "
              "If omitted, the mutation (if any) is applied once at a random allowed "
              "repeat.")
    )
    # orf-output argument
    parser.add_argument(
        "--orf-output",
        default=None,
        help=("If provided, run orfipy_core on the final haplotype sequences in "
              "memory, and write peptide FASTA to this file.")
    )
    # new: minimal protein length in AA
    parser.add_argument(
        "--orf-min-aa",
        type=int,
        default=100,
        help="Minimum peptide length (in amino acids) to report (default=100)."
    )
    # new: optional prefix filter. If used without a value, defaults to 'MTSSV'.
    parser.add_argument(
        "--orf-aa-prefix",
        nargs='?',
        const='MTSSV',
        default=None,
        help=("Filter resulting peptides to only those beginning with this prefix. "
              "If used without a value, defaults to 'MTSSV'. "
              "If omitted entirely, no prefix filter is applied.")
    )
    # NEW: read simulation flag.
    parser.add_argument(
        "--simulate-reads",
        action="store_true",
        help=("If provided, run the read simulation pipeline on the simulated FASTA. "
              "This pipeline will produce an aligned and indexed BAM and gzipped paired "
              "FASTQ files.")
    )
    # NEW: log-level parameter.
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL", "NONE"],
        help=("Set logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL, NONE). "
              "Default is INFO.")
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

    # Process --fixed-lengths logic
    fixed_lengths = None
    if args.fixed_lengths is not None:
        if len(args.fixed_lengths) == 1:
            fixed_lengths = [args.fixed_lengths[0]] * args.num_haplotypes
        elif len(args.fixed_lengths) == args.num_haplotypes:
            fixed_lengths = args.fixed_lengths
        else:
            logging.error(
                "--fixed-lengths must have either 1 value or N=--num-haplotypes"
            )
            sys.exit(1)

    # Simulate haplotypes
    try:
        results = simulate_diploid(
            config=config,
            num_haplotypes=args.num_haplotypes,
            fixed_lengths=fixed_lengths,
            seed=args.seed
        )
        logging.info("Haplotype simulation completed successfully.")
    except Exception as e:
        logging.error("Simulation failed: %s", e)
        sys.exit(1)

    # MUTATION LOGIC
    if args.mutation_name:
        logging.info("Applying mutation: %s", args.mutation_name)
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

    # Write final FASTA
    try:
        with open(args.output, "w") as out_fh:
            for i, (sequence, chain) in enumerate(results, start=1):
                out_fh.write(f">haplotype_{i}\n{sequence}\n")
        logging.info("FASTA output written to %s", args.output)
    except Exception as e:
        logging.error("Writing FASTA failed: %s", e)
        sys.exit(1)

    # Optionally write VNTR structure
    if args.structure_output:
        try:
            with open(args.structure_output, "w") as struct_fh:
                for i, (sequence, chain) in enumerate(results, start=1):
                    chain_str = "-".join(chain)
                    struct_fh.write(f"haplotype_{i}\t{chain_str}\n")
            logging.info("Structure file written to %s", args.structure_output)
        except Exception as e:
            logging.error("Writing structure file failed: %s", e)
            sys.exit(1)

    # ORF logic
    if args.orf_output:
        try:
            run_orf_finder_in_memory(
                results,
                output_pep=args.orf_output,
                orf_min_aa=args.orf_min_aa,
                required_prefix=args.orf_aa_prefix
            )
            logging.info("ORF finding completed; peptide FASTA written to %s",
                         args.orf_output)
        except Exception as e:
            logging.error("ORF finding failed: %s", e)
            sys.exit(1)

    # NEW: Run read simulation pipeline if requested.
    if args.simulate_reads:
        try:
            logging.info("Starting read simulation pipeline.")
            simulate_reads_pipeline(config, args.output)
            logging.info("Read simulation pipeline completed.")
        except Exception as e:
            logging.error("Read simulation pipeline failed: %s", e)
            sys.exit(1)


if __name__ == "__main__":
    main()
