# muc_one_up/simulate.py
import random
from .distribution import sample_repeat_count

# You can place this function in probabilities.py if you prefer.
def pick_next_symbol_no_end(probabilities, current_symbol):
    """
    Like pick_next_repeat, but excludes 6, 6p, 9, or END from possible picks
    (i.e. we don't want to end early).
    Returns:
      symbol (str) or None if no valid next symbol is available.
    """
    if current_symbol not in probabilities:
        return None

    next_options = probabilities[current_symbol]

    # Filter out 6, 6p, 9, END
    forbidden = {"6", "6p", "9", "END"}
    filtered = [(sym, prob) for (sym, prob) in next_options.items() if sym not in forbidden]
    if not filtered:
        return None

    symbols, weights = zip(*filtered)
    return random.choices(symbols, weights=weights, k=1)[0]


def simulate_diploid(config, num_haplotypes=2, fixed_lengths=None, seed=None):
    """
    Main entry point to simulate multiple haplotypes using the config and parameters.
    Returns a list of (haplotype_sequence, haplotype_structure_chain).
    """
    if seed is not None:
        random.seed(seed)

    haplotypes = []

    for i in range(num_haplotypes):
        if fixed_lengths is not None:
            target_length = fixed_lengths[i]
        else:
            target_length = sample_repeat_count(config["length_model"])

        seq, chain = simulate_single_haplotype(config, target_length)
        haplotypes.append((seq, chain))

    return haplotypes


def simulate_single_haplotype(config, target_length, min_length=10):
    """
    Builds a single haplotype by chaining repeats
    according to probabilities, respecting a target length.

    The final 4 repeats are forced to be 6/6p -> 7 -> 8 -> 9 -> right flank.
    We do not allow picking these "end repeats" (6, 6p, 9, END) before that time.

    Returns: (assembled_seq, repeat_chain)
    """
    # 1) Check minimal length
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

    # The repeat at which we forcibly place 6->7->8->9
    final_block_start = target_length - 4

    while True:
        repeat_chain.append(current_symbol)
        assembled_seq += repeats_dict[current_symbol]
        total_repeats += 1

        # If we've arrived at the point where only 4 repeats remain,
        # forcibly place them: [6 or 6p] -> 7 -> 8 -> 9
        if total_repeats == final_block_start:
            # Force 6 or 6p (random)
            forced_6 = random.choice(["6", "6p"])
            repeat_chain.append(forced_6)
            assembled_seq += repeats_dict[forced_6]
            total_repeats += 1

            # Then 7
            repeat_chain.append("7")
            assembled_seq += repeats_dict["7"]
            total_repeats += 1

            # Then 8
            repeat_chain.append("8")
            assembled_seq += repeats_dict["8"]
            total_repeats += 1

            # Then 9
            repeat_chain.append("9")
            assembled_seq += repeats_dict["9"]
            total_repeats += 1

            # Attach right flank, done
            assembled_seq += right_const
            break

        # If for some reason we exceed final_block_start, just end
        if total_repeats >= target_length:
            assembled_seq += right_const
            break

        # Otherwise pick next repeat from a filtered distribution
        next_symbol = pick_next_symbol_no_end(probabilities, current_symbol)
        if next_symbol is None:
            # No valid next symbol => forcibly end
            assembled_seq += right_const
            break

        current_symbol = next_symbol

    return assembled_seq, repeat_chain
