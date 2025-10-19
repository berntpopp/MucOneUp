#!/usr/bin/env python3
"""
snp_integrator.py

This module handles the integration of Single Nucleotide Polymorphisms
into simulated haplotype sequences. It supports:

    - Reading SNPs from a specified file in TSV format
    - Randomly generating SNPs based on user-defined parameters
    - Applying SNPs to sequences while validating reference bases
    - Tracking successfully applied SNPs for reporting

The SNP file format is TSV with the following columns:
    - haplotype_index (1-based): 1 or 2 for diploid
    - position (0-based): Position in the haplotype sequence
    - ref_base: Expected reference base at position (used for validation)
    - alt_base: Alternative base to introduce
"""

import logging
import random
from pathlib import Path
from typing import Any

# Valid DNA bases for SNP generation
DNA_BASES = ["A", "C", "G", "T"]


def parse_snp_file(filepath: str) -> list[dict[str, Any]]:
    """
    Parse a tab-separated SNP file and return a list of SNP definitions.

    The file format should have the following columns:
    - haplotype_index (1-based): 1 or 2 for diploid
    - position (0-based): Position in the haplotype sequence
    - ref_base: Expected reference base at position
    - alt_base: Alternative base to introduce

    Lines starting with # are treated as comments and ignored.

    Args:
        filepath (str): Path to the SNP definition file

    Returns:
        List[Dict[str, Any]]: List of SNP definitions as dictionaries

    Raises:
        ValueError: If the file format is invalid or contains invalid values
    """
    snps = []

    try:
        with Path(filepath).open() as f:
            for line_num, line in enumerate(f, 1):
                # Skip comments or empty lines
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                # Parse TSV fields
                fields = line.split("\t")
                if len(fields) != 4:
                    raise ValueError(
                        f"Line {line_num}: Expected 4 tab-separated fields, found {len(fields)}"
                    )

                # Validate and extract fields
                try:
                    hap_idx = int(fields[0])
                    if hap_idx not in [1, 2]:
                        raise ValueError(
                            f"Line {line_num}: Haplotype index must be 1 or 2, found {hap_idx}"
                        )

                    position = int(fields[1])
                    if position < 0:
                        raise ValueError(
                            f"Line {line_num}: Position must be a non-negative "
                            f"integer, found {position}"
                        )

                    ref_base = fields[2].upper()
                    if ref_base not in DNA_BASES:
                        raise ValueError(
                            f"Line {line_num}: Reference base must be one of "
                            f"{DNA_BASES}, found {ref_base}"
                        )

                    alt_base = fields[3].upper()
                    if alt_base not in DNA_BASES:
                        raise ValueError(
                            f"Line {line_num}: Alternative base must be one of "
                            f"{DNA_BASES}, found {alt_base}"
                        )

                    if ref_base == alt_base:
                        logging.warning(
                            f"Line {line_num}: Reference and alternative bases "
                            f"are identical ({ref_base})"
                        )

                except ValueError as e:
                    # Re-raise with more context
                    raise ValueError(f"Error parsing line {line_num}: {e}") from e

                # Add validated SNP definition
                snps.append(
                    {
                        "haplotype": hap_idx,
                        "position": position,
                        "ref_base": ref_base,
                        "alt_base": alt_base,
                    }
                )

        if not snps:
            logging.warning(f"No valid SNPs found in {filepath}")
        return snps

    except FileNotFoundError as e:
        raise ValueError(f"SNP file not found: {filepath}") from e
    except Exception as e:
        raise ValueError(f"Error reading SNP file {filepath}: {e!s}") from e


def get_vntr_boundaries(
    simulation_results: list[tuple[str, list[str]]], config: dict
) -> list[dict[str, int]]:
    """
    Calculate the boundaries of the VNTR region within each haplotype sequence.

    Args:
        simulation_results (List[Tuple[str, List[str]]]): Simulation results with sequences and chains
        config (dict): Configuration dictionary with constants

    Returns:
        List[Dict[str, int]]: List of dictionaries with 'start' and 'end' positions of VNTR regions
    """
    boundaries = []
    reference_assembly = config.get("reference_assembly", "hg38")

    for seq, _ in simulation_results:
        left_const = config["constants"][reference_assembly]["left"]
        right_const = config["constants"][reference_assembly]["right"]

        vntr_start = len(left_const)
        vntr_end = len(seq) - len(right_const)

        boundaries.append({"start": vntr_start, "end": vntr_end})

    return boundaries


def generate_random_snps(
    sequences: list[str],
    density_per_kb: float,
    region: str = "all",
    target_haplotypes: str = "all",
    vntr_boundaries: list[dict[str, int]] | None = None,
) -> list[dict[str, Any]]:
    """
    Generate random SNPs based on specified parameters.

    Args:
        sequences (List[str]): List of haplotype sequences
        density_per_kb (float): Target density of SNPs per 1000 bp
        region (str, optional): Region to apply SNPs to: "all", "constants_only", or "vntr_only"
        target_haplotypes (str, optional): Which haplotypes to target: "all", "1", or "2"
        vntr_boundaries (List[Dict[str, int]], optional): Boundaries of VNTR regions

    Returns:
        List[Dict[str, any]]: List of randomly generated SNP definitions
    """
    if not sequences:
        return []

    # Validate region parameter
    if region not in ["all", "constants_only", "vntr_only"]:
        raise ValueError(
            f"Invalid region parameter: {region}. Must be 'all', 'constants_only', or 'vntr_only'"
        )

    # Determine target haplotypes (0-based indices)
    target_hap_indices = []
    if target_haplotypes == "all":
        target_hap_indices = list(range(len(sequences)))
    elif target_haplotypes == "1":
        target_hap_indices = [0]
    elif target_haplotypes == "2" and len(sequences) > 1:
        target_hap_indices = [1]
    else:
        raise ValueError(
            f"Invalid target_haplotypes parameter: {target_haplotypes}. Must be 'all', '1', or '2'"
        )

    # Calculate target rate per base
    rate_per_base = density_per_kb / 1000.0

    # Generate SNPs for each target haplotype
    all_snps = []

    for hap_idx in target_hap_indices:
        if hap_idx >= len(sequences):
            continue

        sequence = sequences[hap_idx]
        seq_length = len(sequence)

        # Determine target regions based on region parameter
        valid_positions = set(range(seq_length))

        if region != "all" and vntr_boundaries and hap_idx < len(vntr_boundaries):
            boundaries = vntr_boundaries[hap_idx]
            vntr_start = boundaries["start"]
            vntr_end = boundaries["end"]

            if region == "constants_only":
                # Keep only positions in constant regions
                valid_positions = set(range(vntr_start)) | set(range(vntr_end, seq_length))
            elif region == "vntr_only":
                # Keep only positions in VNTR region
                valid_positions = set(range(vntr_start, vntr_end))

        # Calculate number of SNPs to generate
        target_count = int(len(valid_positions) * rate_per_base)
        logging.debug(
            f"Target SNP count for haplotype {hap_idx + 1}: {target_count} "
            f"(rate: {rate_per_base}, region: {region})"
        )

        if target_count <= 0:
            continue

        # Convert set to list for random sampling
        valid_pos_list = list(valid_positions)

        # Generate unique random positions
        if target_count >= len(valid_pos_list):
            # If target count exceeds available positions, use all positions
            positions = valid_pos_list
            logging.warning(
                f"Requested SNP count ({target_count}) exceeds available "
                f"positions ({len(valid_pos_list)})"
            )
        else:
            # Sample without replacement
            positions = random.sample(valid_pos_list, target_count)

        # Create SNP definitions
        for pos in positions:
            ref_base = sequence[pos].upper()

            # Skip if reference is not a standard base
            if ref_base not in DNA_BASES:
                continue

            # Choose a random alternative base different from reference
            alt_bases = [b for b in DNA_BASES if b != ref_base]
            alt_base = random.choice(alt_bases)

            all_snps.append(
                {
                    "haplotype": hap_idx + 1,  # Convert to 1-based for output
                    "position": pos,
                    "ref_base": ref_base,
                    "alt_base": alt_base,
                }
            )

    return all_snps


def apply_snps_to_sequences(
    haplotype_sequences: list[str],
    snps_to_apply: list[dict[str, Any]],
    skip_reference_check: bool = False,
) -> tuple[list[str], dict[int, list[dict[str, Any]]]]:
    """
    Apply a list of SNPs to the haplotype sequences.

    Args:
        haplotype_sequences (List[str]): Original haplotype sequences
        snps_to_apply (List[Dict[str, any]]): List of SNP definitions to apply
        skip_reference_check (bool): If True, skip checking if reference base matches
                                    (useful for mutated sequences where positions or bases
                                    may have changed)

    Returns:
        Tuple[List[str], Dict[int, List[Dict[str, Any]]]]:
            - Modified sequences
            - Dictionary mapping haplotype index (0-based) to successfully applied SNPs
    """
    if not haplotype_sequences or not snps_to_apply:
        return haplotype_sequences, {i: [] for i in range(len(haplotype_sequences))}

    # Convert sequences to lists for modification
    mutable_sequences = [list(seq) for seq in haplotype_sequences]

    # Track successfully applied SNPs
    applied_snps: dict[int, list[dict[str, Any]]] = {i: [] for i in range(len(haplotype_sequences))}

    # Apply each SNP
    for snp in snps_to_apply:
        # Convert 1-based haplotype index to 0-based
        hap_idx = snp["haplotype"] - 1

        # Validate haplotype index
        if hap_idx < 0 or hap_idx >= len(mutable_sequences):
            logging.warning(
                f"SNP specifies invalid haplotype index {snp['haplotype']} (0-based: {hap_idx})"
            )
            continue

        sequence = mutable_sequences[hap_idx]
        position = snp["position"]

        # Validate position
        if position < 0 or position >= len(sequence):
            logging.warning(
                f"SNP position {position} out of bounds for haplotype "
                f"{hap_idx + 1} (length: {len(sequence)})"
            )
            continue

        # Get reference base information
        actual_ref = sequence[position].upper()
        expected_ref = snp["ref_base"].upper()

        # Validate reference base (unless skip_reference_check is True)
        if not skip_reference_check and actual_ref != expected_ref:
            logging.warning(
                f"Reference base mismatch at position {position} in haplotype "
                f"{hap_idx + 1}: expected '{expected_ref}', found '{actual_ref}'"
            )
            continue

        # Apply the SNP
        sequence[position] = snp["alt_base"]

        # Record successful application
        applied_snps[hap_idx].append(snp)
        logging.debug(
            f"Applied SNP at position {position} in haplotype {hap_idx + 1}: "
            f"{expected_ref} -> {snp['alt_base']}"
        )

    # Convert lists back to strings
    modified_sequences = ["".join(seq) for seq in mutable_sequences]

    return modified_sequences, applied_snps


def write_snps_to_file(snps: list[dict[str, Any]], filepath: str) -> None:
    """
    Write a list of SNPs to a file in the standard TSV format.

    Args:
        snps (List[Dict[str, any]]): List of SNP definitions
        filepath (str): Output file path

    Raises:
        IOError: If file cannot be written
    """
    try:
        with Path(filepath).open("w") as f:
            f.write("# MucOneUp SNP File\n")
            f.write(
                "# Columns: haplotype_index (1-based)\tposition (0-based)\tref_base\talt_base\n"
            )

            for snp in snps:
                f.write(
                    f"{snp['haplotype']}\t{snp['position']}\t{snp['ref_base']}\t{snp['alt_base']}\n"
                )

        logging.info(f"Wrote {len(snps)} SNPs to {filepath}")
    except Exception as e:
        raise OSError(f"Failed to write SNP file {filepath}: {e!s}") from e
