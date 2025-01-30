# muc_one_up/simulate.py

import random
from .distribution import sample_repeat_count
from .probabilities import pick_next_repeat

def simulate_diploid(config, num_haplotypes=2, fixed_lengths=None, seed=None):
    """
    Main entry point to simulate multiple haplotypes using the config and parameters.
    Returns a list of (haplotype_sequence, haplotype_structure_chain).

    :param config: dict with repeats, probabilities, constants, length_model
    :param num_haplotypes: number of haplotypes to generate
    :param fixed_lengths: list of integers (same length as num_haplotypes), or None
    :param seed: seed for reproducibility
    """
    if seed is not None:
        random.seed(seed)

    haplotypes = []
    for i in range(num_haplotypes):
        # Decide which length to use for this haplotype
        if fixed_lengths is not None:
            repeat_count = fixed_lengths[i]
        else:
            # fallback: random distribution
            repeat_count = sample_repeat_count(config["length_model"])

        seq, chain = simulate_single_haplotype(config, repeat_count)
        haplotypes.append((seq, chain))

    return haplotypes


def simulate_single_haplotype(config, target_length):
    """
    Builds a single haplotype by chaining repeats
    according to probabilities, respecting a target length.
    The last 4 repeats must be 6->7->8->9 (or 6p->7->8->9).
    """
    left_const = config["constants"]["left"]
    right_const = config["constants"]["right"]

    repeats_dict = config["repeats"]
    probabilities = config["probabilities"]

    assembled_seq = left_const
    repeat_chain = []

    current_symbol = "1"  # Typically start from repeat '1'
    total_repeats = 0

    while True:
        repeat_chain.append(current_symbol)
        assembled_seq += repeats_dict[current_symbol]
        total_repeats += 1

        # If we've reached the target length minus 4, we forcibly pick 6 or 6p,
        # then manually chain 7->8->9->END in the next steps.
        if total_repeats == target_length - 4:
            # 1) Choose between '6' or '6p'. You can do a 50/50,
            #    or read from a parameter if you like.
            next_symbol = random.choice(["6", "6p"])

            # Add this forced repeat
            repeat_chain.append(next_symbol)
            assembled_seq += repeats_dict[next_symbol]
            total_repeats += 1

            # 2) Next is always 7
            repeat_chain.append("7")
            assembled_seq += repeats_dict["7"]
            total_repeats += 1

            # 3) Next is always 8
            repeat_chain.append("8")
            assembled_seq += repeats_dict["8"]
            total_repeats += 1

            # 4) Next is always 9
            repeat_chain.append("9")
            assembled_seq += repeats_dict["9"]
            total_repeats += 1

            # Finally, attach right flank
            assembled_seq += right_const
            break

        # If we've reached or exceeded target_length (and haven't forced the last 4 yet)
        # we forcibly end. This covers edge cases (like if target_length < 4).
        if total_repeats >= target_length:
            assembled_seq += right_const
            break

        # Otherwise pick next_symbol via probabilities
        next_symbol = pick_next_repeat(probabilities, current_symbol)

        if next_symbol == "END":
            assembled_seq += right_const
            break
        else:
            current_symbol = next_symbol

    return assembled_seq, repeat_chain
