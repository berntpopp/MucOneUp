"""Mutation application engine for simulated VNTR haplotypes.

This module applies targeted mutations to repeat units within simulated
haplotypes. Supports four mutation types: insert, delete, replace, and
delete_insert. Mutations are validated against allowed_repeats and can
operate in strict mode (reject invalid targets) or permissive mode (auto-convert).

Key Functions:
    apply_mutations: Apply named mutation to specific haplotype positions
    validate_allowed_repeats: Ensure mutation targets are valid repeat symbols
    apply_changes_to_repeat: Execute insertion/deletion/replacement operations

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
            results=[HaplotypeResult(seq1, chain1), HaplotypeResult(seq2, chain2)],
            mutation_name="dupC",
            targets=[MutationTarget(1, 25)]  # 1-based indexing
        )

Notes:
    - Positions use 1-based indexing (matching biological conventions)
    - Mutated repeats are marked via RepeatUnit.mutated flag
    - Strict mode prevents silent auto-conversions (recommended for reproducibility)

See Also:
    - config.py: CONFIG_SCHEMA for mutation definitions
    - simulate.py: Haplotype generation
"""

import logging
import random as _random_module

from .type_defs import (
    AssemblyConstants,
    ConfigDict,
    DNASequence,
    HaplotypeResult,
    MutationDefinition,
    MutationName,
    MutationTarget,
    RepeatUnit,
)


def validate_allowed_repeats(mutation_def: MutationDefinition, valid_symbols: set[str]) -> set[str]:
    """Validate that allowed_repeats contains only valid repeat symbols.

    Args:
        mutation_def: Mutation definition from config["mutations"][name]
        valid_symbols: Set of valid repeat symbol strings

    Returns:
        Set of validated allowed repeat symbols

    Raises:
        ValueError: If any allowed repeat is not a valid repeat symbol
    """
    allowed_repeats = set(mutation_def.get("allowed_repeats", []))

    invalid_repeats = allowed_repeats - valid_symbols
    if invalid_repeats:
        valid_repeats_str = ", ".join(sorted(valid_symbols))
        invalid_repeats_str = ", ".join(sorted(invalid_repeats))
        raise ValueError(
            f"Invalid repeats in allowed_repeats: {invalid_repeats_str}. "
            f"Valid repeats are: {valid_repeats_str}"
        )

    return allowed_repeats


def apply_mutations(
    config: ConfigDict,
    results: list[HaplotypeResult],
    mutation_name: MutationName,
    targets: list[MutationTarget],
    rng: _random_module.Random | None = None,
) -> tuple[list[HaplotypeResult], dict[int, list[tuple[int, str]]]]:
    """Apply named mutation to specific haplotype positions.

    Applies the mutation defined in config["mutations"][mutation_name] to each
    target position. Validates allowed_repeats, handles strict mode enforcement,
    and tracks mutated VNTR unit sequences for reporting.

    Args:
        config: Configuration dict containing mutations section
        results: List of HaplotypeResult for each haplotype
        mutation_name: Key in config["mutations"] defining the mutation
        targets: List of MutationTarget using 1-based indexing

    Returns:
        Tuple containing:
        - updated_results: Modified haplotypes with mutations applied
        - mutated_units: Dict mapping haplotype index (1-based) to list of
          (repeat_index, mutated_sequence) tuples

    Raises:
        ValueError: If mutation name not found, target indices invalid, or
                   strict mode rejects invalid repeat symbol

    Note:
        Mutated repeats are marked via RepeatUnit.mutated flag.
        Positions use 1-based indexing (haplotype and repeat).
    """
    if "mutations" not in config:
        raise ValueError("No 'mutations' section in config; cannot apply mutations.")

    if mutation_name not in config["mutations"]:
        raise ValueError(f"Mutation '{mutation_name}' not found in config['mutations'].")

    mutation_def = config["mutations"][mutation_name]
    ref_assembly = config.get("reference_assembly", "hg38")
    constants: AssemblyConstants = config["constants"][ref_assembly]
    left_len = len(constants["left"])
    repeats_dict = config["repeats"]
    valid_symbols = set(repeats_dict.keys())

    # Validate that allowed_repeats only contains valid repeat symbols
    try:
        allowed_repeats = validate_allowed_repeats(mutation_def, valid_symbols)
    except ValueError as e:
        logging.error(f"In mutation '{mutation_name}': {e}")
        raise

    # Check if strict mode is enabled for this mutation
    strict_mode = mutation_def.get("strict_mode", False)
    changes = mutation_def["changes"]

    _rng = rng if rng is not None else _random_module

    updated_results = [
        HaplotypeResult(
            sequence=hr.sequence,
            chain=[RepeatUnit(ru.symbol, ru.mutated) for ru in hr.chain],
            _repeat_lengths=list(hr._repeat_lengths) if hr._repeat_lengths is not None else None,
        )
        for hr in results
    ]
    lengths_by_hap: list[list[int]] = []
    for hr in results:
        if hr._repeat_lengths is not None and len(hr._repeat_lengths) == len(hr.chain):
            lengths_by_hap.append(list(hr._repeat_lengths))
        else:
            lengths = [len(repeats_dict[ru.symbol]) for ru in hr.chain]
            if sum(1 for ru in hr.chain if ru.mutated) == 1:
                right_flank_len = (
                    len(constants.get("right", ""))
                    if (
                        not hr.chain
                        or hr.chain[-1].symbol == "9"
                        or hr.sequence.endswith(constants.get("right", ""))
                    )
                    else 0
                )
                delta = len(hr.sequence) - left_len - right_flank_len - sum(lengths)
                mut_idx = next(i for i, ru in enumerate(hr.chain) if ru.mutated)
                lengths[mut_idx] += delta
            lengths_by_hap.append(lengths)
    mutated_units: dict[
        int, list[tuple[int, str]]
    ] = {}  # key: haplotype (1-based), value: list of (repeat_index, mutated_unit_sequence)

    # Apply the mutation to each specified target.
    for target in targets:
        hap_i = target.haplotype_index
        rep_i = target.repeat_index
        h_i = hap_i - 1
        r_i = rep_i - 1

        if h_i < 0 or h_i >= len(updated_results):
            raise ValueError(f"Haplotype index {hap_i} out of range (1..{len(updated_results)}).")

        hr = updated_results[h_i]

        if r_i < 0 or r_i >= len(hr.chain):
            raise ValueError(f"Repeat index {rep_i} out of range (1..{len(hr.chain)}).")

        start = left_len + sum(lengths_by_hap[h_i][:r_i])
        stop = start + lengths_by_hap[h_i][r_i]

        # Handle the case where the current symbol is not in allowed_repeats
        current_symbol_clean = hr.chain[r_i].symbol
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
            new_symbol = _rng.choice(sorted(allowed_repeats))
            logging.warning(
                "Forcing change at haplotype %d, repeat %d: %s -> %s for mutation '%s'",
                hap_i,
                rep_i,
                current_symbol_clean,
                new_symbol,
                mutation_name,
            )
            hr.chain[r_i] = RepeatUnit(new_symbol)
            replacement_seq = repeats_dict[new_symbol]
            hr.sequence = hr.sequence[:start] + replacement_seq + hr.sequence[stop:]
            lengths_by_hap[h_i][r_i] = len(replacement_seq)
            stop = start + len(replacement_seq)

        # Apply the specified changes and capture the mutated unit.
        seq, chain, mutated_repeat = apply_changes_to_repeat(
            hr.sequence,
            hr.chain,
            r_i,
            changes,
            config,
            mutation_name,
            start=start,
            stop=stop,
        )

        # Mark the mutated repeat.
        chain[r_i] = RepeatUnit(chain[r_i].symbol, mutated=True)
        lengths_by_hap[h_i][r_i] = len(mutated_repeat)
        updated_results[h_i] = HaplotypeResult(
            sequence=seq, chain=chain, _repeat_lengths=list(lengths_by_hap[h_i])
        )
        logging.info(
            "Applied mutation '%s' at haplotype %d, repeat %d",
            mutation_name,
            hap_i,
            rep_i,
        )

        # Record the mutated unit sequence.
        hap_units = mutated_units.setdefault(hap_i, [])
        for idx_u, (existing_rep_i, _) in enumerate(hap_units):
            if existing_rep_i == rep_i:
                hap_units[idx_u] = (rep_i, mutated_repeat)
                break
        else:
            hap_units.append((rep_i, mutated_repeat))

    return updated_results, mutated_units


def apply_changes_to_repeat(
    seq: DNASequence,
    chain: list[RepeatUnit],
    repeat_index: int,
    changes: list[dict[str, int | str]],
    config: ConfigDict,
    mutation_name: MutationName,
    start: int,
    stop: int,
) -> tuple[DNASequence, list[RepeatUnit], DNASequence]:
    """Modify repeat unit sequence using mutation change operations.

    Applies a series of change operations (insert, delete, replace, delete_insert)
    to the repeat unit at the specified index. Operations use 1-based positioning
    within the repeat unit.

    Args:
        seq: Full haplotype sequence before mutation
        chain: Repeat chain as list of RepeatUnit
        repeat_index: Index (0-based) of the repeat to modify
        changes: List of change dicts with 'type', 'start', 'end', 'sequence' keys
        config: Configuration dict for assembly constants and repeat lookups
        mutation_name: Name of mutation (for error messages)
        start: 0-based start index of the repeat in seq
        stop: 0-based stop index of the repeat in seq

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
    current_repeat_seq = seq[start:stop]
    repeat_chars = list(current_repeat_seq)

    for change in changes:
        repeat_length = len(repeat_chars)
        ctype = change.get("type")
        change_start = change.get("start")
        change_end = change.get("end")
        insertion_str = change.get("sequence", "")

        if change_start is None or change_end is None or ctype is None:
            raise ValueError(f"Malformed change in mutation '{mutation_name}': missing fields.")

        # Type narrowing: after None check, these must be int (from config validation)
        if not isinstance(change_start, int) or not isinstance(change_end, int):
            raise ValueError(f"Invalid types for start/end in mutation '{mutation_name}'")
        if not isinstance(insertion_str, str):
            insertion_str = str(insertion_str)

        start_idx = change_start - 1
        end_idx = change_end - 1

        if ctype == "insert":
            if not (0 <= start_idx <= repeat_length):
                raise ValueError(
                    f"Insert out of bounds in mutation '{mutation_name}': start={change_start}, repeat length={repeat_length}"
                )
            repeat_chars[start_idx:start_idx] = list(insertion_str)
        elif ctype == "delete":
            if start_idx < 0 or end_idx >= repeat_length or start_idx > end_idx:
                raise ValueError(
                    f"Delete out of bounds in mutation '{mutation_name}': start={change_start}, end={change_end}, repeat length={repeat_length}"
                )
            del repeat_chars[start_idx : end_idx + 1]
        elif ctype == "replace":
            if start_idx < 0 or end_idx >= repeat_length or start_idx > end_idx:
                raise ValueError(
                    f"Replace out of bounds in mutation '{mutation_name}': start={change_start}, end={change_end}, repeat length={repeat_length}"
                )
            repeat_chars[start_idx : end_idx + 1] = list(insertion_str)
        elif ctype == "delete_insert":
            # For delete_insert, we interpret the provided start and end as the boundaries that
            # are to be retained, deleting the bases strictly between them.
            if start_idx < 0 or end_idx >= repeat_length or start_idx >= end_idx:
                raise ValueError(
                    f"delete_insert out of bounds in mutation '{mutation_name}': start={change_start}, end={change_end}, repeat length={repeat_length}"
                )
            # Delete bases strictly between start_idx and end_idx.
            del repeat_chars[start_idx + 1 : end_idx]
            # Insert the provided sequence at position start_idx+1.
            repeat_chars[start_idx + 1 : start_idx + 1] = list(insertion_str)
        else:
            raise ValueError(f"Unknown mutation type '{ctype}' in mutation '{mutation_name}'")

    mutated_repeat = "".join(repeat_chars)
    new_seq = seq[:start] + mutated_repeat + seq[stop:]
    return new_seq, chain, mutated_repeat
