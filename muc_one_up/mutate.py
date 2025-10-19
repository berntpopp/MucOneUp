"""Mutation application engine for simulated VNTR haplotypes.

This module applies targeted mutations to repeat units within simulated
haplotypes. Supports four mutation types: insert, delete, replace, and
delete_insert. Mutations are validated against allowed_repeats and can
operate in strict mode (reject invalid targets) or permissive mode (auto-convert).

Key Functions:
    apply_mutations: Apply named mutation to specific haplotype positions
    validate_allowed_repeats: Ensure mutation targets are valid repeat symbols
    apply_changes_to_repeat: Execute insertion/deletion/replacement operations
    rebuild_haplotype_sequence: Reassemble sequence after mutations

Mutation Structure:
    Mutations are defined in config["mutations"][name] with:
    - allowed_repeats: Valid repeat symbols for this mutation
    - strict_mode: Boolean controlling auto-conversion behavior
    - changes: List of operations (type, start, end, sequence)

Example:
    Apply dupC mutation to haplotype 1, repeat position 25::

        from muc_one_up.mutate import apply_mutations

        config = load_config("config.json")
        results, mutated_units = apply_mutations(
            config=config,
            results=[(seq1, chain1), (seq2, chain2)],
            mutation_name="dupC",
            targets=[(1, 25)]  # 1-based indexing
        )

Notes:
    - Positions use 1-based indexing (matching biological conventions)
    - Mutated repeats are marked with 'm' suffix in chains
    - Strict mode prevents silent auto-conversions (recommended for reproducibility)

See Also:
    - config.py: CONFIG_SCHEMA for mutation definitions
    - simulate.py: Haplotype generation
"""

import logging
import random

from .type_defs import (
    ConfigDict,
    DNASequence,
    HaplotypeList,
    MutationDefinition,
    MutationName,
    MutationTargets,
    RepeatChain,
)


def validate_allowed_repeats(mutation_def: MutationDefinition, config: ConfigDict) -> set[str]:
    """Validate that allowed_repeats contains only valid repeat symbols.

    Args:
        mutation_def: Mutation definition from config["mutations"][name]
        config: Configuration dict containing the "repeats" section

    Returns:
        Set of validated allowed repeat symbols

    Raises:
        ValueError: If any allowed repeat is not a valid repeat symbol
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
    config: ConfigDict,
    results: HaplotypeList,
    mutation_name: MutationName,
    targets: MutationTargets,
) -> tuple[HaplotypeList, dict[int, list[tuple[int, str]]]]:
    """Apply named mutation to specific haplotype positions.

    Applies the mutation defined in config["mutations"][mutation_name] to each
    target position. Validates allowed_repeats, handles strict mode enforcement,
    and tracks mutated VNTR unit sequences for reporting.

    Args:
        config: Configuration dict containing mutations section
        results: List of (sequence, chain) tuples for each haplotype
        mutation_name: Key in config["mutations"] defining the mutation
        targets: List of (haplotype_idx, repeat_idx) tuples using 1-based indexing

    Returns:
        Tuple containing:
        - updated_results: Modified haplotypes with mutations applied
        - mutated_units: Dict mapping haplotype index (1-based) to list of
          (repeat_index, mutated_sequence) tuples

    Raises:
        ValueError: If mutation name not found, target indices invalid, or
                   strict mode rejects invalid repeat symbol

    Note:
        Mutated repeats are marked with 'm' suffix in the chain.
        Positions use 1-based indexing (haplotype and repeat).
    """
    if "mutations" not in config:
        raise ValueError("No 'mutations' section in config; cannot apply mutations.")

    if mutation_name not in config["mutations"]:
        raise ValueError(f"Mutation '{mutation_name}' not found in config['mutations'].")

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
    mutated_units: dict[
        int, list[tuple[int, str]]
    ] = {}  # key: haplotype (1-based), value: list of (repeat_index, mutated_unit_sequence)

    # Apply the mutation to each specified target.
    for hap_i, rep_i in targets:
        hap_index = hap_i - 1
        repeat_index = rep_i - 1

        if hap_index < 0 or hap_index >= len(updated_results):
            raise ValueError(f"Haplotype index {hap_i} out of range (1..{len(updated_results)}).")

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


def rebuild_haplotype_sequence(chain: RepeatChain, config: ConfigDict) -> DNASequence:
    """Rebuild haplotype sequence from repeat chain and flanking constants.

    Concatenates left constant + repeat units + right constant (if terminal repeat is 9).
    Strips mutation markers ('m') from symbols when looking up sequences.

    Args:
        chain: List of repeat symbols, possibly with 'm' suffix marking mutations
        config: Configuration dict containing 'constants' and 'repeats' sections

    Returns:
        Reassembled haplotype DNA sequence
    """
    reference_assembly = config.get("reference_assembly", "hg38")
    left_const = str(config["constants"][reference_assembly]["left"])
    right_const = str(config["constants"][reference_assembly]["right"])
    repeats_dict = config["repeats"]

    seq = left_const
    for sym in chain:
        pure_sym = sym.replace("m", "")
        seq += str(repeats_dict[pure_sym])
    if chain[-1].replace("m", "") == "9":
        seq += right_const
    return seq


def apply_changes_to_repeat(
    seq: DNASequence,
    chain: RepeatChain,
    repeat_index: int,
    changes: list[dict[str, int | str]],
    config: ConfigDict,
    mutation_name: MutationName,
) -> tuple[DNASequence, RepeatChain, DNASequence]:
    """Modify repeat unit sequence using mutation change operations.

    Applies a series of change operations (insert, delete, replace, delete_insert)
    to the repeat unit at the specified index. Operations use 1-based positioning
    within the repeat unit.

    Args:
        seq: Full haplotype sequence before mutation
        chain: Repeat chain identifying unit positions
        repeat_index: Index (0-based) of the repeat to modify
        changes: List of change dicts with 'type', 'start', 'end', 'sequence' keys
        config: Configuration dict for assembly constants and repeat lookups
        mutation_name: Name of mutation (for error messages)

    Returns:
        Tuple containing:
        - new_seq: Updated haplotype sequence with mutation applied
        - chain: Repeat chain (unmodified, marker added separately)
        - mutated_repeat: The modified repeat unit sequence

    Raises:
        ValueError: If change coordinates are out of bounds or type is unknown

    Note:
        Change operations:
        - insert: Insert sequence at position (inclusive)
        - delete: Delete bases from start to end (inclusive)
        - replace: Replace bases from start to end with sequence
        - delete_insert: Delete between boundaries, insert at boundary
    """
    repeats_dict = config["repeats"]
    reference_assembly = config.get("reference_assembly", "hg38")
    left_const_len = len(config["constants"][reference_assembly]["left"])

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
            raise ValueError(f"Malformed change in mutation '{mutation_name}': missing fields.")

        # Type narrowing: after None check, these must be int (from config validation)
        if not isinstance(start, int) or not isinstance(end, int):
            raise ValueError(f"Invalid types for start/end in mutation '{mutation_name}'")
        if not isinstance(insertion_str, str):
            insertion_str = str(insertion_str)

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
            raise ValueError(f"Unknown mutation type '{ctype}' in mutation '{mutation_name}'")

    mutated_repeat = "".join(repeat_chars)
    new_seq = seq[:start_pos] + mutated_repeat + seq[end_pos + 1 :]
    return new_seq, chain, mutated_repeat
