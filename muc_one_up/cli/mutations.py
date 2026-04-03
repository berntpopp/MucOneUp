"""
Mutation application functions for MucOneUp CLI.

Single Responsibility: Parse mutation targets and apply mutations to haplotypes.
"""

import logging
import random as _random_module

from ..exceptions import MutationError, ValidationError
from ..mutate import apply_mutations
from ..type_defs import HaplotypeResult, MutationTarget


def parse_mutation_targets(
    mutation_targets_list: list[str | tuple[int, int]],
) -> list[MutationTarget]:
    """
    Parse mutation target strings into MutationTarget objects.

    Args:
        mutation_targets_list: List of "haplotype,repeat" strings

    Returns:
        List of MutationTarget objects

    Raises:
        ValidationError: If target format is invalid
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
            mutation_positions.append(MutationTarget(hap_i, rep_i))
        except Exception as e:
            raise ValidationError(f"Invalid --mutation-targets format: '{t}' ({e})") from e
    return mutation_positions


def find_random_mutation_target(
    results: list[HaplotypeResult],
    config: dict,
    mutation_name: str,
    rng: _random_module.Random | None = None,
) -> list[MutationTarget]:
    """
    Find random valid targets for mutation.

    Args:
        results: List of HaplotypeResult objects
        config: Configuration dictionary
        mutation_name: Name of mutation to apply

    Returns:
        List with single MutationTarget

    Raises:
        MutationError: If mutation not found or no valid targets exist
    """
    if "mutations" not in config or mutation_name not in config["mutations"]:
        raise MutationError(
            f"Mutation '{mutation_name}' not in config['mutations']. "
            f"Available mutations: {list(config.get('mutations', {}).keys())}"
        )

    mut_def = config["mutations"][mutation_name]
    allowed_repeats = set(mut_def["allowed_repeats"])
    possible_targets = []

    for hap_idx, hr in enumerate(results, start=1):
        for rep_idx, unit in enumerate(hr.chain, start=1):
            if unit.symbol in allowed_repeats:
                possible_targets.append(MutationTarget(hap_idx, rep_idx))

    if not possible_targets:
        raise MutationError(
            f"No repeats match 'allowed_repeats' for '{mutation_name}'. "
            f"Allowed repeats: {mut_def['allowed_repeats']}"
        )

    _rng = rng if rng is not None else _random_module
    return [_rng.choice(possible_targets)]


def apply_mutation_pipeline(
    args,
    config,
    results,
    mutation_name,
    dual_mutation_mode,
    mutation_pair,
    rng: _random_module.Random | None = None,
) -> tuple[
    list[HaplotypeResult], list[HaplotypeResult] | None, dict | None, list[MutationTarget] | None
]:
    """
    Apply mutations based on configuration.

    Single Responsibility: Mutation application logic.

    Args:
        args: Command-line arguments
        config: Configuration dictionary
        results: Haplotype results
        mutation_name: Name of mutation to apply
        dual_mutation_mode: Whether in dual mutation mode
        mutation_pair: Pair of mutation names (for dual mode)

    Returns:
        Tuple of (results, mutated_results, mutated_units, mutation_positions)

    Raises:
        MutationError: If mutation application fails
        ValidationError: If mutation targets are invalid
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
            mutation_positions = find_random_mutation_target(
                results, config, mutation_pair[1], rng=rng
            )

        try:
            mutated_results, mutated_units = apply_mutations(
                config=config,
                results=[HaplotypeResult(hr.sequence, list(hr.chain)) for hr in results],
                mutation_name=mutation_pair[1],
                targets=mutation_positions,
                rng=rng,
            )
            logging.info("Dual mutation applied for mutated version.")
            return normal_results, mutated_results, mutated_units, mutation_positions
        except Exception as e:
            raise MutationError(f"Dual mutation application failed: {e}") from e
    else:
        # Single mode: mutate in place
        if args.mutation_targets:
            mutation_positions = parse_mutation_targets(args.mutation_targets)
        else:
            mutation_positions = find_random_mutation_target(
                results, config, mutation_name, rng=rng
            )

        try:
            results, mutated_units = apply_mutations(
                config=config,
                results=results,
                mutation_name=mutation_name,
                targets=mutation_positions,
                rng=rng,
            )
            logging.info("Mutation applied at targets: %s", mutation_positions)
            return results, mutated_results, mutated_units, mutation_positions
        except Exception as e:
            raise MutationError(f"Mutation application failed: {e}") from e
