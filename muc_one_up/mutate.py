import random
import logging
from typing import Dict, List, Tuple, Any, Set


def validate_allowed_repeats(
    mutation_def: Dict[str, Any], config: Dict[str, Any]
) -> Set[str]:
    """
    Validate that the allowed_repeats in a mutation definition are valid repeat symbols.

    :param mutation_def: Mutation definition from config["mutations"][name].
    :param config: The entire config dict containing the "repeats" section.
    :return: Set of validated allowed repeat symbols.
    :raises ValueError: If any allowed repeat is not a valid repeat symbol.
    """
    allowed_repeats = set(mutation_def.get("allowed_repeats", []))
    valid_repeats = set(config["repeats"].keys())

    invalid_repeats = allowed_repeats - valid_repeats
    if invalid_repeats:
        valid_repeats_str = ", ".join(sorted(valid_repeats))
        invalid_repeats_str = ", ".join(sorted(invalid_repeats))
        raise ValueError(
            f"Invalid repeats in allowed_repeats: {invalid_repeats_str}. "
            f"Valid repeats are: {valid_repeats_str}"
        )

    return allowed_repeats


def apply_mutations(
    config: Dict[str, Any],
    results: List[Tuple[str, List[str]]],
    mutation_name: str,
    targets: List[Tuple[int, int]],
) -> Tuple[List[Tuple[str, List[str]]], Dict[int, List[Tuple[int, str]]]]:
    """
    Apply a single named mutation to one or more haplotypes at specific repeat indices,
    and record the mutated VNTR unit(s).

    :param config: The entire config dict.
    :param results: List of (sequence, chain) for each haplotype.
    :param mutation_name: Key in config["mutations"].
    :param targets: List of (haplotype_idx, repeat_idx) (1-based indexing).
    :return: Tuple (updated_results, mutated_units) where:
             - updated_results: List of (sequence, chain) for each haplotype after mutation.
             - mutated_units: Dict mapping haplotype index (1-based) to a list of tuples
               (repeat_index, mutated_unit_sequence).
    :raises ValueError: if configuration or target indices are invalid.
    """
    if "mutations" not in config:
        raise ValueError("No 'mutations' section in config; cannot apply mutations.")

    if mutation_name not in config["mutations"]:
        raise ValueError(
            f"Mutation '{mutation_name}' not found in config['mutations']."
        )

    mutation_def = config["mutations"][mutation_name]

    # Validate that allowed_repeats only contains valid repeat symbols
    try:
        allowed_repeats = validate_allowed_repeats(mutation_def, config)
    except ValueError as e:
        logging.error(f"In mutation '{mutation_name}': {e}")
        raise

    # Check if strict mode is enabled for this mutation
    strict_mode = mutation_def.get("strict_mode", False)
    changes = mutation_def["changes"]

    updated_results = list(results)
    mutated_units = (
        {}
    )  # key: haplotype (1-based), value: list of (repeat_index, mutated_unit_sequence)

    # Apply the mutation to each specified target.
    for hap_i, rep_i in targets:
        hap_index = hap_i - 1
        repeat_index = rep_i - 1

        if hap_index < 0 or hap_index >= len(updated_results):
            raise ValueError(
                f"Haplotype index {hap_i} out of range (1..{len(updated_results)})."
            )

        seq, chain = updated_results[hap_index]

        if repeat_index < 0 or repeat_index >= len(chain):
            raise ValueError(f"Repeat index {rep_i} out of range (1..{len(chain)}).")

        current_symbol = chain[repeat_index]

        # Handle the case where the current symbol is not in allowed_repeats
        current_symbol_clean = current_symbol.replace("m", "")
        if current_symbol_clean not in allowed_repeats:
            if not allowed_repeats:
                raise ValueError(
                    f"Mutation '{mutation_name}' has no allowed_repeats, cannot fix symbol."
                )

            # If in strict mode, raise an error instead of forcing a change
            if strict_mode:
                valid_repeats_str = ", ".join(sorted(allowed_repeats))
                raise ValueError(
                    f"Cannot apply mutation '{mutation_name}' at haplotype {hap_i}, "
                    f"repeat {rep_i}: Repeat symbol '{current_symbol_clean}' is not in "
                    f"allowed_repeats: {valid_repeats_str}. Set strict_mode=false to "
                    f"allow automatic conversion."
                )

            # In non-strict mode, force a change to a random allowed repeat
            new_symbol = random.choice(list(allowed_repeats))
            logging.warning(
                "Forcing change at haplotype %d, repeat %d: %s -> %s for mutation '%s'",
                hap_i,
                rep_i,
                current_symbol_clean,
                new_symbol,
                mutation_name,
            )
            chain[repeat_index] = new_symbol
            seq = rebuild_haplotype_sequence(chain, config)
            updated_results[hap_index] = (seq, chain)

        current_symbol = chain[repeat_index]

        # Apply the specified changes and capture the mutated unit.
        seq, chain, mutated_repeat = apply_changes_to_repeat(
            seq, chain, repeat_index, changes, config, mutation_name
        )

        # Mark the mutated repeat.
        chain[repeat_index] = chain[repeat_index] + "m"
        updated_results[hap_index] = (seq, chain)
        logging.info(
            "Applied mutation '%s' at haplotype %d, repeat %d",
            mutation_name,
            hap_i,
            rep_i,
        )

        # Record the mutated unit sequence.
        mutated_units.setdefault(hap_i, []).append((rep_i, mutated_repeat))

    return updated_results, mutated_units


def rebuild_haplotype_sequence(chain: List[str], config: Dict[str, Any]) -> str:
    """
    Rebuild the haplotype sequence from the chain of repeats and constant flanks.

    :param chain: List of repeat symbols (with possible appended 'm').
    :param config: Configuration dict containing 'constants' and 'repeats'.
    :return: Reassembled haplotype sequence.
    """
    left_const = config["constants"]["left"]
    right_const = config["constants"]["right"]
    repeats_dict = config["repeats"]

    seq = left_const
    for sym in chain:
        pure_sym = sym.replace("m", "")
        seq += repeats_dict[pure_sym]
    if chain[-1].replace("m", "") == "9":
        seq += right_const
    return seq


def apply_changes_to_repeat(
    seq: str,
    chain: List[str],
    repeat_index: int,
    changes: List[Dict[str, Any]],
    config: Dict[str, Any],
    mutation_name: str,
) -> Tuple[str, List[str], str]:
    """
    Modify the substring of 'seq' corresponding to chain[repeat_index] using the
    list of changes from the mutation definition.

    :param seq: Original haplotype sequence.
    :param chain: The repeat chain.
    :param repeat_index: Index (0-based) of the repeat to modify.
    :param changes: List of change dictionaries.
    :param config: Configuration dict.
    :param mutation_name: Name of the mutation.
    :return: Tuple (new_seq, chain, mutated_repeat) after applying the changes.
    """
    repeats_dict = config["repeats"]
    left_const_len = len(config["constants"]["left"])

    offset = left_const_len
    for i in range(repeat_index):
        pure_sym = chain[i].replace("m", "")
        offset += len(repeats_dict[pure_sym])

    current_symbol = chain[repeat_index].replace("m", "")
    repeat_seq = repeats_dict[current_symbol]
    start_pos = offset
    end_pos = offset + len(repeat_seq) - 1  # inclusive

    repeat_chars = list(repeat_seq)
    repeat_length = len(repeat_chars)

    for change in changes:
        ctype = change.get("type")
        start = change.get("start")
        end = change.get("end")
        insertion_str = change.get("sequence", "")

        if start is None or end is None or ctype is None:
            raise ValueError(
                f"Malformed change in mutation '{mutation_name}': missing fields."
            )

        start_idx = start - 1
        end_idx = end - 1

        if ctype == "insert":
            if not (0 <= start_idx <= repeat_length):
                raise ValueError(
                    f"Insert out of bounds in mutation '{mutation_name}': start={start}, repeat length={repeat_length}"
                )
            repeat_chars[start_idx:start_idx] = list(insertion_str)
        elif ctype == "delete":
            if start_idx < 0 or end_idx >= repeat_length or start_idx > end_idx:
                raise ValueError(
                    f"Delete out of bounds in mutation '{mutation_name}': start={start}, end={end}, repeat length={repeat_length}"
                )
            del repeat_chars[start_idx : end_idx + 1]
        elif ctype == "replace":
            if start_idx < 0 or end_idx >= repeat_length or start_idx > end_idx:
                raise ValueError(
                    f"Replace out of bounds in mutation '{mutation_name}': start={start}, end={end}, repeat length={repeat_length}"
                )
            repeat_chars[start_idx : end_idx + 1] = list(insertion_str)
        elif ctype == "delete_insert":
            # For delete_insert, we interpret the provided start and end as the boundaries that
            # are to be retained, deleting the bases strictly between them.
            if start_idx < 0 or end_idx >= repeat_length or start_idx >= end_idx:
                raise ValueError(
                    f"delete_insert out of bounds in mutation '{mutation_name}': start={start}, end={end}, repeat length={repeat_length}"
                )
            # Delete bases strictly between start_idx and end_idx.
            del repeat_chars[start_idx + 1 : end_idx]
            # Insert the provided sequence at position start_idx+1.
            repeat_chars[start_idx + 1 : start_idx + 1] = list(insertion_str)
        else:
            raise ValueError(
                f"Unknown mutation type '{ctype}' in mutation '{mutation_name}'"
            )

    mutated_repeat = "".join(repeat_chars)
    new_seq = seq[:start_pos] + mutated_repeat + seq[end_pos + 1 :]
    return new_seq, chain, mutated_repeat
