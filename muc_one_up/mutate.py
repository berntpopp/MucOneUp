# muc_one_up/mutate.py

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

        # 1) If current_symbol not in allowed_repeats => forcibly change the symbol
        #    to something that is allowed. We'll pick the first allowed symbol, or you can define logic.
        if current_symbol not in allowed_repeats:
            if not allowed_repeats:
                raise ValueError(
                    f"Mutation '{mutation_name}' has no allowed_repeats, cannot fix symbol."
                )
            new_symbol = allowed_repeats[0]
            # update the chain
            chain[repeat_index] = new_symbol
            # We also need to reconstruct the entire haplotype sequence
            # with that new repeat symbol. 
            seq = rebuild_haplotype_sequence(chain, config)  
            # update in results
            updated_results[hap_index] = (seq, chain)
        
        # Now re-fetch current_symbol, because it might have changed
        current_symbol = chain[repeat_index]

        # 2) Actually apply the changes in 'changes' to just that repeat substring
        seq, chain = apply_changes_to_repeat(
            seq, chain, repeat_index, changes, config, mutation_name
        )

        # 3) Mark the mutated repeat in the chain with 'm'
        chain[repeat_index] = chain[repeat_index] + "m"

        # update results
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
    # We do not forcibly put "m" in the symbol here because "m" is a marker for structure
    # But if the chain symbol literally changed, e.g. "Xm", we can handle it carefully
    for sym in chain:
        # If user appended "m", the actual underlying repeat is sym.replace("m", "")
        pure_sym = sym.replace("m","")
        seq += repeats_dict[pure_sym]
    # but watch out for if you have forced final 6->7->8->9->END logic in your pipeline...
    # This is a simple approach. 
    # If your code always appended right_const, do that here:
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

    current_symbol = chain[repeat_index].replace("m","")
    repeat_seq = repeats_dict[current_symbol]
    start_pos = offset
    end_pos = offset + len(repeat_seq) - 1

    # Convert the repeat substring to a list for easy mutation
    repeat_chars = list(repeat_seq)

    for change in changes:
        ctype = change["type"]  # "insert", "delete", or "replace"
        start = change["start"]  # 1-based
        end = change["end"]      # 1-based

        # Convert to 0-based for the substring
        start_idx = start - 1
        end_idx = end - 1

        if ctype == "insert":
            # Insert 'sequence' at start_idx
            insertion_str = change.get("sequence", "")
            # Example: "1-2 G" => Insert G *between* base 1 and 2 => effectively at index 1
            # but if end=2, you might interpret it as "after the first base"? 
            # Here we keep it simple: insert at start_idx
            repeat_chars[start_idx:start_idx] = list(insertion_str)

        elif ctype == "delete":
            # Delete from start_idx to end_idx (inclusive)
            del repeat_chars[start_idx : end_idx + 1]

        elif ctype == "replace":
            # Replace from start_idx to end_idx with change["sequence"]
            repl_str = change.get("sequence", "")
            repeat_chars[start_idx : end_idx + 1] = list(repl_str)

        else:
            raise ValueError(f"Unknown mutation type '{ctype}' in mutation '{mutation_name}'")

    # Rejoin the mutated repeat
    mutated_repeat = "".join(repeat_chars)

    # Now we must assemble the entire haplotype again:
    #  seq[0:start_pos] + mutated_repeat + seq[end_pos+1:]
    new_seq = seq[:start_pos] + mutated_repeat + seq[end_pos+1:]

    # Return updated
    return (new_seq, chain)
