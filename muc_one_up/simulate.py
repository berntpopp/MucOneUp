# muc_one_up/simulate.py

import random
from .distribution import sample_repeat_count
from .probabilities import pick_next_repeat


def simulate_diploid(config, num_haplotypes=2, fixed_length=None, seed=None):
    """
    Main entry point to simulate multiple haplotypes using the config and parameters.
    Returns a list of (haplotype_sequence, haplotype_structure_chain).
    """
    if seed is not None:
        random.seed(seed)

    haplotypes = []
    for _ in range(num_haplotypes):
        repeat_count = fixed_length if fixed_length is not None else sample_repeat_count(config["length_model"])
        seq, chain = simulate_single_haplotype(config, repeat_count)
        haplotypes.append((seq, chain))

    return haplotypes


def simulate_single_haplotype(config, target_length):
    """
    Builds a single haplotype by chaining repeats
    according to probabilities, respecting a target length.

    Returns:
        assembled_seq (str): full DNA sequence (left flank + repeats + right flank)
        repeat_chain (list of str): the sequence of repeat labels used
    """
    left_const = config["constants"]["left"]
    right_const = config["constants"]["right"]

    repeats_dict = config["repeats"]
    probabilities = config["probabilities"]

    # We'll build the final DNA sequence in this variable
    assembled_seq = left_const

    # Keep track of which repeat symbols we add
    repeat_chain = []

    current_symbol = "1"  # Example assumption: always start from repeat '1'
    total_repeats = 0

    while True:
        # Add the current symbol to the chain
        repeat_chain.append(current_symbol)

        # Add the sequence for the current repeat
        assembled_seq += repeats_dict[current_symbol]
        total_repeats += 1

        # If we've reached or exceeded the target_length, we force or bias an END
        if total_repeats >= target_length:
            next_symbol = pick_next_repeat(probabilities, current_symbol, force_end=True)
        else:
            next_symbol = pick_next_repeat(probabilities, current_symbol)

        if next_symbol == "END":
            assembled_seq += right_const
            break
        else:
            current_symbol = next_symbol

    return assembled_seq, repeat_chain
