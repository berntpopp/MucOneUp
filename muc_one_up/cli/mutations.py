"""
Mutation application functions for MucOneUp CLI.

Single Responsibility: Parse mutation targets and apply mutations to haplotypes.
"""

import logging
import random
import sys

from ..mutate import apply_mutations


def parse_mutation_targets(mutation_targets_list: list[str]) -> list[tuple[int, int]]:
    """
    Parse mutation target strings into tuples.

    Args:
        mutation_targets_list: List of "haplotype,repeat" strings

    Returns:
        List of (haplotype_idx, repeat_idx) tuples
    """
    mutation_positions = []
    for t in mutation_targets_list:
        try:
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
            logging.error("Invalid --mutation-targets format: '%s' (%s)", t, e)
            sys.exit(1)
    return mutation_positions


def find_random_mutation_target(
    results: list[tuple[str, list[str]]], config: dict, mutation_name: str
) -> list[tuple[int, int]]:
    """
    Find random valid targets for mutation.

    Args:
        results: List of (sequence, chain) tuples
        config: Configuration dictionary
        mutation_name: Name of mutation to apply

    Returns:
        List with single (haplotype_idx, repeat_idx) tuple
    """
    if "mutations" not in config or mutation_name not in config["mutations"]:
        raise ValueError(f"Mutation '{mutation_name}' not in config['mutations']")

    mut_def = config["mutations"][mutation_name]
    allowed_repeats = set(mut_def["allowed_repeats"])
    possible_targets = []

    for hap_idx, (_seq, chain) in enumerate(results, start=1):
        for rep_idx, sym in enumerate(chain, start=1):
            pure_sym = sym.replace("m", "")
            if pure_sym in allowed_repeats:
                possible_targets.append((hap_idx, rep_idx))

    if not possible_targets:
        raise ValueError(f"No repeats match 'allowed_repeats' for '{mutation_name}'")

    return [random.choice(possible_targets)]


def apply_mutation_pipeline(
    args, config, results, mutation_name, dual_mutation_mode, mutation_pair
) -> tuple[list, list | None, dict | None, list | None]:
    """
    Apply mutations based on configuration.

    Single Responsibility: Mutation application logic.

    Returns:
        Tuple of (results, mutated_results, mutated_units, mutation_positions)
    """
    mutated_results = None
    mutated_units = None
    mutation_positions = None

    if not mutation_name:
        return results, mutated_results, mutated_units, mutation_positions

    logging.info("Applying mutation: %s", mutation_name)

    if dual_mutation_mode:
        # Dual mode: keep normal, create mutated
        normal_results = results

        if args.mutation_targets:
            mutation_positions = parse_mutation_targets(args.mutation_targets)
        else:
            mutation_positions = find_random_mutation_target(results, config, mutation_pair[1])

        try:
            mutated_results, mutated_units = apply_mutations(
                config=config,
                results=[(seq, chain.copy()) for seq, chain in results],
                mutation_name=mutation_pair[1],
                targets=mutation_positions,
            )
            logging.info("Dual mutation applied for mutated version.")
            return normal_results, mutated_results, mutated_units, mutation_positions
        except Exception as e:
            logging.error("Mutation failed: %s", e)
            sys.exit(1)
    else:
        # Single mode: mutate in place
        if args.mutation_targets:
            mutation_positions = parse_mutation_targets(args.mutation_targets)
        else:
            try:
                mutation_positions = find_random_mutation_target(results, config, mutation_name)
            except Exception as e:
                logging.error("Random-target mutation failed: %s", e)
                sys.exit(1)

        try:
            results, mutated_units = apply_mutations(
                config=config,
                results=results,
                mutation_name=mutation_name,
                targets=mutation_positions,
            )
            logging.info("Mutation applied at targets: %s", mutation_positions)
            return results, mutated_results, mutated_units, mutation_positions
        except Exception as e:
            logging.error("Mutation failed: %s", e)
            sys.exit(1)
