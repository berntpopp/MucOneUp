# muc_one_up/simulate.py

import random
import logging
from .distribution import sample_repeat_count


def pick_next_symbol_no_end(probabilities, current_symbol):
    """
    Pick the next symbol from the probability table for the current symbol,
    excluding forbidden symbols (6, 6p, 9, END).

    :param probabilities: Dict with probability options.
    :param current_symbol: The current repeat symbol.
    :return: Next symbol (str) or None if no valid symbol is available.
    """
    if current_symbol not in probabilities:
        logging.debug("Current symbol '%s' not in probabilities.", current_symbol)
        return None

    next_options = probabilities[current_symbol]
    forbidden = {"6", "6p", "9", "END"}
    filtered = [(sym, prob) for (sym, prob) in next_options.items() if sym not in forbidden]
    if not filtered:
        logging.debug("No valid next symbol available from '%s'.", current_symbol)
        return None

    symbols, weights = zip(*filtered)
    chosen = random.choices(symbols, weights=weights, k=1)[0]
    logging.debug("Picked next symbol '%s' from '%s'.", chosen, current_symbol)
    return chosen


def simulate_diploid(config, num_haplotypes=2, fixed_lengths=None, seed=None):
    """
    Simulate multiple haplotypes using the provided config and parameters.

    :param config: Configuration dict.
    :param num_haplotypes: Number of haplotypes to simulate.
    :param fixed_lengths: Optional list of fixed repeat lengths.
    :param seed: Optional random seed.
    :return: List of tuples (haplotype_sequence, repeat_chain).
    """
    if seed is not None:
        random.seed(seed)
        logging.debug("Random seed set to %d", seed)

    haplotypes = []
    for i in range(num_haplotypes):
        if fixed_lengths is not None:
            target_length = fixed_lengths[i]
        else:
            target_length = sample_repeat_count(config["length_model"])
        logging.info("Simulating haplotype %d with target length %d", i + 1, target_length)
        seq, chain = simulate_single_haplotype(config, target_length)
        haplotypes.append((seq, chain))
    return haplotypes


def simulate_single_haplotype(config, target_length, min_length=10):
    """
    Build a single haplotype by chaining repeats according to probabilities,
    while respecting a target length.

    The final 4 repeats are forced to be 6/6p -> 7 -> 8 -> 9 and then the right
    flank is appended.

    :param config: Configuration dict.
    :param target_length: Desired total number of repeats.
    :param min_length: Minimum allowed repeats.
    :return: Tuple (assembled_seq, repeat_chain).
    :raises ValueError: if target_length is below the minimum.
    """
    if target_length < min_length:
        raise ValueError(
            f"Requested target_length={target_length} but minimum is {min_length}. Aborting."
        )

    left_const = config["constants"]["left"]
    right_const = config["constants"]["right"]
    probabilities = config["probabilities"]
    repeats_dict = config["repeats"]

    assembled_seq = left_const
    repeat_chain = []
    current_symbol = "1"
    total_repeats = 0
    final_block_start = target_length - 4

    while True:
        repeat_chain.append(current_symbol)
        assembled_seq += repeats_dict[current_symbol]
        total_repeats += 1

        if total_repeats == final_block_start:
            forced_6 = random.choice(["6", "6p"])
            logging.debug("Forcing terminal block: %s, 7, 8, 9", forced_6)
            repeat_chain.append(forced_6)
            assembled_seq += repeats_dict[forced_6]
            total_repeats += 1

            repeat_chain.append("7")
            assembled_seq += repeats_dict["7"]
            total_repeats += 1

            repeat_chain.append("8")
            assembled_seq += repeats_dict["8"]
            total_repeats += 1

            repeat_chain.append("9")
            assembled_seq += repeats_dict["9"]
            total_repeats += 1

            assembled_seq += right_const
            break

        if total_repeats >= target_length:
            assembled_seq += right_const
            break

        next_symbol = pick_next_symbol_no_end(probabilities, current_symbol)
        if next_symbol is None:
            logging.debug("No valid next symbol; appending right constant.")
            assembled_seq += right_const
            break

        current_symbol = next_symbol

    logging.info("Simulated haplotype with %d repeats.", total_repeats)
    return assembled_seq, repeat_chain
