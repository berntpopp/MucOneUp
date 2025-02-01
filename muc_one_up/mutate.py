# muc_one_up/mutate.py

import random
import logging


def apply_mutations(config, results, mutation_name, targets):
    """
    Apply a single named mutation to one or more haplotypes at specific repeat indices.

    :param config: The entire config dict.
    :param results: List of (sequence, chain) for each haplotype.
    :param mutation_name: Key in config["mutations"].
    :param targets: List of (haplotype_idx, repeat_idx) (1-based indexing).
    :return: Updated results (same structure).
    :raises ValueError: if configuration or target indices are invalid.
    """
    if "mutations" not in config:
        raise ValueError("No 'mutations' section in config; cannot apply mutations.")

    if mutation_name not in config["mutations"]:
        raise ValueError(
            f"Mutation '{mutation_name}' not found in config['mutations']."
        )

    mutation_def = config["mutations"][mutation_name]
    allowed_repeats = mutation_def["allowed_repeats"]
    changes = mutation_def["changes"]

    updated_results = list(results)

    # Apply the mutation to each specified target.
    for (hap_i, rep_i) in targets:
        hap_index = hap_i - 1
        repeat_index = rep_i - 1

        if hap_index < 0 or hap_index >= len(updated_results):
            raise ValueError(
                f"Haplotype index {hap_i} out of range (1..{len(updated_results)})."
            )

        seq, chain = updated_results[hap_index]

        if repeat_index < 0 or repeat_index >= len(chain):
            raise ValueError(
                f"Repeat index {rep_i} out of range (1..{len(chain)})."
            )

        current_symbol = chain[repeat_index]

        # Force a change if the current symbol is not allowed.
        if current_symbol not in allowed_repeats:
            if not allowed_repeats:
                raise ValueError(
                    f"Mutation '{mutation_name}' has no allowed_repeats, cannot fix symbol."
                )
            new_symbol = random.choice(allowed_repeats)
            logging.debug(
                "Forcing change at haplotype %d, repeat %d: %s -> %s",
                hap_i, rep_i, current_symbol, new_symbol
            )
            chain[repeat_index] = new_symbol
            seq = rebuild_haplotype_sequence(chain, config)
            updated_results[hap_index] = (seq, chain)

        current_symbol = chain[repeat_index]

        # Apply the specified changes.
        seq, chain = apply_changes_to_repeat(
            seq, chain, repeat_index, changes, config, mutation_name
        )

        # Mark the mutated repeat.
        chain[repeat_index] = chain[repeat_index] + "m"
        updated_results[hap_index] = (seq, chain)
        logging.info(
            "Applied mutation '%s' at haplotype %d, repeat %d",
            mutation_name, hap_i, rep_i
        )

    return updated_results


def rebuild_haplotype_sequence(chain, config):
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


def apply_changes_to_repeat(seq, chain, repeat_index, changes, config, mutation_name):
    """
    Modify the substring of 'seq' corresponding to chain[repeat_index] using the
    list of changes from the mutation definition.

    :param seq: Original haplotype sequence.
    :param chain: The repeat chain.
    :param repeat_index: Index (0-based) of the repeat to modify.
    :param changes: List of change dictionaries.
    :param config: Configuration dict.
    :param mutation_name: Name of the mutation.
    :return: Tuple (new_seq, chain) after applying the changes.
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
                    f"Insert out of bounds in mutation '{mutation_name}': "
                    f"start={start}, repeat length={repeat_length}"
                )
        elif ctype in ("delete", "replace"):
            if start_idx < 0 or end_idx >= repeat_length or start_idx > end_idx:
                raise ValueError(
                    f"{ctype.capitalize()} out of bounds in mutation '{mutation_name}': "
                    f"start={start}, end={end}, repeat length={repeat_length}"
                )
        else:
            raise ValueError(
                f"Unknown mutation type '{ctype}' in mutation '{mutation_name}'"
            )

        if ctype == "insert":
            repeat_chars[start_idx:start_idx] = list(insertion_str)
        elif ctype == "delete":
            del repeat_chars[start_idx:end_idx + 1]
        elif ctype == "replace":
            repeat_chars[start_idx:end_idx + 1] = list(insertion_str)

    mutated_repeat = "".join(repeat_chars)
    new_seq = seq[:start_pos] + mutated_repeat + seq[end_pos + 1:]
    return new_seq, chain
