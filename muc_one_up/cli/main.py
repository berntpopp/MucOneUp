"""
Main entry point for MucOneUp CLI.

Single Responsibility: Argument parsing, logging configuration, and top-level orchestration.
"""

import argparse
import logging

from ..version import __version__
from .config import (
    determine_simulation_mode,
    process_mutation_config,
    setup_configuration,
)
from .orchestration import run_single_simulation_iteration


def build_parser():
    """Build and return the argument parser."""
    parser = argparse.ArgumentParser(
        prog="MucOneUp", description="Tool to simulate MUC1 VNTR diploid references."
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s " + __version__,
        help="Print the version number and exit.",
    )
    parser.add_argument("--config", required=True, help="Path to JSON configuration file.")
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

    simulation_mode_group = parser.add_mutually_exclusive_group()
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
    parser.add_argument(
        "--simulate-reads",
        nargs="?",
        const="illumina",
        choices=["illumina", "ont"],
        help=(
            "Run the read simulation pipeline on the simulated FASTA. "
            "Specify 'illumina' for short reads (via reseq/WeSSim) or 'ont' for Oxford Nanopore long reads (via NanoSim). "
            "If used without a value, defaults to 'illumina'. "
            "This pipeline produces an aligned and indexed BAM and FASTQ files."
        ),
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL", "NONE"],
        help="Set logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL, NONE). Default is INFO.",
    )
    parser.add_argument(
        "--reference-assembly",
        choices=["hg19", "hg38"],
        help="Specify reference assembly (hg19 or hg38). Overrides the setting in config file.",
    )

    snp_group = parser.add_argument_group(
        "SNP Integration", "Options for integrating SNPs into simulated haplotypes"
    )
    snp_source_group = snp_group.add_mutually_exclusive_group()
    snp_source_group.add_argument(
        "--snp-input-file",
        type=str,
        help=(
            "Path to a TSV file defining specific SNPs to apply. Format: "
            "haplotype_index\tposition\tref_base\talt_base. "
            "Cannot be used with random SNP options."
        ),
    )
    snp_source_group.add_argument(
        "--random-snps",
        action="store_true",
        help="Enable random SNP generation. Cannot be used with --snp-input-file.",
    )
    snp_group.add_argument(
        "--random-snp-density",
        type=float,
        help=(
            "Target density of SNPs, specified as SNPs per 1000 bp "
            "(e.g., 0.5 for 1 SNP every 2000 bp). "
            "Required if --random-snps is used."
        ),
    )
    snp_group.add_argument(
        "--random-snp-output-file",
        type=str,
        help=(
            "Path to write the list of generated random SNPs in TSV format. "
            "Required if --random-snps is used."
        ),
    )
    snp_group.add_argument(
        "--random-snp-region",
        choices=["all", "constants_only", "vntr_only"],
        default="constants_only",
        help=(
            "Apply random SNPs to specified region: all, constants_only (default), "
            "or vntr_only. Requires --random-snps."
        ),
    )
    snp_group.add_argument(
        "--random-snp-haplotypes",
        choices=["all", "1", "2"],
        default="all",
        help=(
            "Apply random SNPs to specified haplotypes: all (default, applied "
            "independently), 1, or 2. Requires --random-snps."
        ),
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
        logging.basicConfig(level=level, format="%(asctime)s - %(levelname)s - %(message)s")
        logging.info("Logging configured at level: %s", level_str.upper())


def main():
    """
    Main entry point for the CLI.

    Orchestrates the simulation pipeline following SOLID principles.
    Each major concern is delegated to focused helper functions.

    Total lines: <100 (down from 801!)
    """
    # Parse arguments and configure logging
    parser = build_parser()
    args = parser.parse_args()
    configure_logging(args.log_level)

    # Load configuration and setup output
    config, out_dir, out_base = setup_configuration(args)

    # Determine simulation mode
    simulation_configs, predefined_chains, structure_mutation_info = determine_simulation_mode(
        args, config
    )

    # Process mutation configuration
    dual_mutation_mode, mutation_pair, mutation_name = process_mutation_config(
        args, structure_mutation_info
    )

    # Run simulations
    for sim_index, fixed_conf in enumerate(simulation_configs, start=1):
        run_single_simulation_iteration(
            args,
            config,
            out_dir,
            out_base,
            sim_index,
            fixed_conf,
            predefined_chains,
            dual_mutation_mode,
            mutation_pair,
            structure_mutation_info,
        )


if __name__ == "__main__":
    main()
