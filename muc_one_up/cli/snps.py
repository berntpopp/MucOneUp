"""
SNP integration functions for MucOneUp CLI.

Single Responsibility: Integrate SNPs from files or generate random SNPs.
DRY: Unified SNP integration for both dual and single mutation modes.
"""

import logging
import sys

from ..snp_integrator import (
    apply_snps_to_sequences,
    generate_random_snps,
    get_vntr_boundaries,
    parse_snp_file,
    write_snps_to_file,
)


def integrate_snps_unified(
    args,
    config,
    results: list[tuple[str, list[str]]],
    skip_reference_check: bool = False,
) -> tuple[list[tuple[str, list[str]]], dict]:
    """
    Unified SNP integration function - eliminates duplication.

    Single Responsibility: SNP integration (file-based or random).
    DRY: Used for both dual and single mutation modes.

    Args:
        args: Parsed command-line arguments
        config: Configuration dictionary
        results: List of (sequence, chain) tuples
        skip_reference_check: Skip reference base validation (for mutated sequences)

    Returns:
        Tuple of (modified_results, applied_snp_info)
    """
    applied_snp_info = {}
    snps_from_source = []
    generated_snp_list_for_output = []

    current_sequences = [seq for seq, chain in results]

    # Check for SNP integration options
    if args.snp_input_file:
        try:
            snps_from_source = parse_snp_file(args.snp_input_file)
            logging.info(f"Read {len(snps_from_source)} SNPs from {args.snp_input_file}")
        except Exception as e:
            logging.error(f"Failed to parse SNP input file: {e}. Skipping SNP integration.")
            snps_from_source = []

    elif args.random_snps:
        # Validate required parameters
        if not args.random_snp_density:
            logging.error("--random-snp-density is required when --random-snps is used")
            sys.exit(1)
        if not args.random_snp_output_file:
            logging.error("--random-snp-output-file is required when --random-snps is used")
            sys.exit(1)

        # Get VNTR boundaries if region filtering is used
        vntr_bounds = get_vntr_boundaries(results, config)

        # Generate random SNPs
        generated_snp_list_for_output = generate_random_snps(
            current_sequences,
            args.random_snp_density,
            args.random_snp_region,
            args.random_snp_haplotypes,
            vntr_bounds,
        )
        snps_from_source = generated_snp_list_for_output
        logging.info(f"Generated {len(snps_from_source)} random SNPs")

        # Write random SNPs to output file
        try:
            write_snps_to_file(generated_snp_list_for_output, args.random_snp_output_file)
            logging.info(
                f"Wrote {len(generated_snp_list_for_output)} SNPs to {args.random_snp_output_file}"
            )
        except Exception as e:
            logging.error(f"Failed to write generated SNP file: {e}")

    # Apply SNPs if we have any from either source
    if snps_from_source:
        modified_sequences, applied_snp_info = apply_snps_to_sequences(
            current_sequences, snps_from_source, skip_reference_check=skip_reference_check
        )

        # Replace original sequences with modified ones
        for i, (_seq, chain) in enumerate(results):
            results[i] = (modified_sequences[i], chain)

        logging.info(f"Applied {sum(len(v) for v in applied_snp_info.values())} SNPs to sequences")

    return results, applied_snp_info
