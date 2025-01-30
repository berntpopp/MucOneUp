# muc_one_up/mutate.py

import random

def apply_mutations(config, results, mutation_name, targets):
    """
    Apply a single named mutation to one or more haplotypes at specific repeat indices.

    :param config: the entire config dict
    :param results: list of (sequence, chain) for each haplotype
    :param mutation_name: key in config["mutations"]
    :param targets: list of (haplotype_idx, repeat_idx), 1-based
    :return: updated results (same structure)
    """
    if "mutations" not in config:
        raise ValueError("No 'mutations' section in config; cannot apply mutations.")

    if mutation_name not in config["mutations"]:
        raise ValueError(f"Mutation '{mutation_name}' not found in config['mutations'].")

    mutation_def = config["mutations"][mutation_name]
    allowed_repeats = mutation_def["allowed_repeats"]
    changes = mutation_def["changes"]  # a list of dicts describing insert/delete/replace

    updated_results = list(results)  # make a copy if you prefer

    # Each target => apply mutation once
    for (hap_i, rep_i) in targets:
        # Convert 1-based to 0-based
        hap_index = hap_i - 1
        repeat_index = rep_i - 1

        if hap_index < 0 or hap_index >= len(updated_results):
            raise ValueError(f"Haplotype index {hap_i} out of range (1..{len(updated_results)}).")
        
        seq, chain = updated_results[hap_index]

        if repeat_index < 0 or repeat_index >= len(chain):
            raise ValueError(f"Repeat index {rep_i} out of range (1..{len(chain)}).")

        current_symbol = chain[repeat_index]

        # 1) If current_symbol not in allowed_repeats => forcibly change
        #    to a random choice among allowed_repeats (instead of always [0]).
        if current_symbol not in allowed_repeats:
            if not allowed_repeats:
                raise ValueError(
                    f"Mutation '{mutation_name}' has no allowed_repeats, cannot fix symbol."
                )
            new_symbol = random.choice(allowed_repeats)
            chain[repeat_index] = new_symbol

            # Rebuild haplotype sequence with that new repeat symbol
            seq = rebuild_haplotype_sequence(chain, config)
            updated_results[hap_index] = (seq, chain)
        
        # Refresh current_symbol (it may have changed)
        current_symbol = chain[repeat_index]

        # 2) Actually apply the changes in 'changes' to just that repeat substring
        seq, chain = apply_changes_to_repeat(
            seq, chain, repeat_index, changes, config, mutation_name
        )

        # 3) Mark the mutated repeat in the chain with 'm'
        chain[repeat_index] = chain[repeat_index] + "m"

        # Final update for this haplotype
        updated_results[hap_index] = (seq, chain)

    return updated_results


def rebuild_haplotype_sequence(chain, config):
    """
    Rebuild entire haplotype sequence from chain of repeats
    plus left & right constants. 
    Used if we forcibly changed a repeat symbol.
    """
    left_const = config["constants"]["left"]
    right_const = config["constants"]["right"]
    repeats_dict = config["repeats"]

    seq = left_const
    for sym in chain:
        # If user appended "m", the actual underlying repeat is sym.replace("m","")
        pure_sym = sym.replace("m", "")
        seq += repeats_dict[pure_sym]

    # If the chain ends with '9' (minus any 'm'), we attach the right const
    if chain[-1].replace("m","") == "9":
        seq += right_const

    return seq


def apply_changes_to_repeat(seq, chain, repeat_index, changes, config, mutation_name):
    """
    Modifies the substring of 'seq' that corresponds to chain[repeat_index].
    Uses the 'changes' list from the mutation definition.

    Returns (new_seq, new_chain).
    """
    repeats_dict = config["repeats"]
    left_const_len = len(config["constants"]["left"])

    # 1) find offset in 'seq' where this repeat begins
    offset = left_const_len
    for i in range(repeat_index):
        pure_sym = chain[i].replace("m","")
        offset += len(repeats_dict[pure_sym])

    # Identify the original repeat substring
    current_symbol = chain[repeat_index].replace("m","")
    repeat_seq = repeats_dict[current_symbol]
    start_pos = offset
    end_pos = offset + len(repeat_seq) - 1  # inclusive

    # We'll manipulate a list of chars
    repeat_chars = list(repeat_seq)
    repeat_length = len(repeat_chars)

    for change in changes:
        ctype = change.get("type")
        start = change.get("start")
        end = change.get("end")
        insertion_str = change.get("sequence", "")

        if start is None or end is None or ctype is None:
            raise ValueError(f"Malformed change in mutation '{mutation_name}': missing fields.")

        # Convert to 0-based (Python) indices
        start_idx = start - 1
        end_idx = end - 1

        # Check coordinate validity. 
        # For inserts, it may be valid to insert at start_idx == repeat_length
        # but for delete/replace, we typically need start_idx <= end_idx < repeat_length.
        if ctype == "insert":
            # Insert can happen anywhere from 0..repeat_length
            if not (0 <= start_idx <= repeat_length):
                raise ValueError(
                    f"Insert out of bounds in mutation '{mutation_name}': "
                    f"start={start}, repeat length={repeat_length}"
                )

        elif ctype in ("delete", "replace"):
            # We require 0 <= start_idx <= end_idx < repeat_length
            if start_idx < 0 or end_idx >= repeat_length or start_idx > end_idx:
                raise ValueError(
                    f"{ctype.capitalize()} out of bounds in mutation '{mutation_name}': "
                    f"start={start}, end={end}, repeat length={repeat_length}"
                )
        else:
            raise ValueError(f"Unknown mutation type '{ctype}' in mutation '{mutation_name}'")

        # Perform the actual edit
        if ctype == "insert":
            # Insert insertion_str at start_idx
            repeat_chars[start_idx:start_idx] = list(insertion_str)

        elif ctype == "delete":
            del repeat_chars[start_idx : end_idx + 1]

        elif ctype == "replace":
            repeat_chars[start_idx : end_idx + 1] = list(insertion_str)

    # Rebuild the mutated repeat as a string
    mutated_repeat = "".join(repeat_chars)

    # Reassemble the full haplotype:
    # everything before start_pos + mutated_repeat + everything after end_pos
    new_seq = seq[:start_pos] + mutated_repeat + seq[end_pos+1:]

    return (new_seq, chain)
