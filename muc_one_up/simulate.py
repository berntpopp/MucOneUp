# muc_one_up/simulate.py

import logging
import random

from .distribution import sample_repeat_count
from .type_defs import (
    ConfigDict,
    DNASequence,
    Haplotype,
    HaplotypeList,
    MutationTargets,
    ProbabilitiesDict,
    RepeatChain,
)


def pick_next_symbol_no_end(
    probabilities: ProbabilitiesDict,
    current_symbol: str,
) -> str | None:
    """Pick next symbol excluding forbidden terminal symbols.

    Excludes forbidden symbols (6, 6p, 9, END) to prevent premature termination
    of the repeat chain.

    Args:
        probabilities: Probability table for state transitions
        current_symbol: The current repeat symbol

    Returns:
        Next symbol as string, or None if no valid symbol is available
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

    symbols, weights = zip(*filtered, strict=True)
    chosen = random.choices(symbols, weights=weights, k=1)[0]
    logging.debug("Picked next symbol '%s' from '%s'.", chosen, current_symbol)
    return str(chosen) if chosen is not None else None


def assemble_haplotype_from_chain(chain: RepeatChain, config: ConfigDict) -> DNASequence:
    """Assemble complete haplotype sequence from repeat chain.

    Concatenates left constant + repeat units + right constant to form
    the final haplotype sequence.

    Args:
        chain: List of repeat symbols to assemble (e.g., ['1', '2', '7', '8', '9'])
        config: Configuration with repeats and constants

    Returns:
        Assembled haplotype DNA sequence

    Raises:
        ValueError: If repeat symbol not found in config
        KeyError: If reference assembly constants not found
    """
    reference_assembly = config.get("reference_assembly", "hg38")
    left_const = config["constants"][reference_assembly]["left"]
    right_const = config["constants"][reference_assembly]["right"]
    repeats_dict = config["repeats"]

    # Start with the left constant
    assembled_seq = left_const

    # Add each repeat unit according to the chain
    for symbol in chain:
        # If symbol has mutation marker ('m'), we still use the base symbol for lookup
        base_symbol = symbol.rstrip("m")
        if base_symbol not in repeats_dict:
            raise ValueError(f"Symbol '{base_symbol}' not found in config repeats.")

        assembled_seq += str(repeats_dict[base_symbol])

    # Add right constant
    assembled_seq += right_const

    logging.info(f"Assembled haplotype with {len(chain)} repeats.")
    return assembled_seq  # type: ignore[no-any-return]


def simulate_from_chains(
    predefined_chains: list[RepeatChain],
    config: ConfigDict,
    mutation_name: str | None = None,
    mutation_targets: MutationTargets | None = None,
) -> HaplotypeList:
    """Simulate haplotypes from predefined repeat chains.

    Takes predefined repeat chains and assembles them into haplotype sequences.
    Optionally marks positions with mutation markers ('m' suffix) if
    mutation_targets are specified.

    Args:
        predefined_chains: List of repeat chains, one per haplotype
        config: Configuration dictionary
        mutation_name: Optional mutation name (for logging)
        mutation_targets: Optional list of (haplotype_idx, repeat_idx) tuples
                         using 1-based indexing

    Returns:
        List of (haplotype_sequence, repeat_chain) tuples

    Note:
        mutation_targets use 1-based indexing for both haplotype and repeat positions
    """
    haplotypes = []

    for i, chain in enumerate(predefined_chains):
        # Create a copy of the chain to avoid modifying the original
        working_chain = chain.copy()

        # Apply mutations if specified in the mutation_targets
        if mutation_name and mutation_targets:
            for hap_idx, repeat_idx in mutation_targets:
                # Check if this mutation applies to the current haplotype
                if hap_idx == i + 1:  # User provides 1-indexed haplotype numbers
                    # Convert from 1-indexed to 0-indexed for the repeat position
                    zero_based_idx = repeat_idx - 1
                    # Ensure the repeat index is valid
                    if 0 <= zero_based_idx < len(working_chain):
                        # Only add mutation marker if it doesn't already have one
                        if not working_chain[zero_based_idx].endswith("m"):
                            logging.info(
                                f"Applying mutation '{mutation_name}' to haplotype {i+1}, "
                                f"repeat position {repeat_idx} (0-based: {zero_based_idx})"
                            )
                            working_chain[zero_based_idx] = working_chain[zero_based_idx] + "m"
                    else:
                        logging.warning(
                            f"Mutation target repeat index {repeat_idx} (0-based: {zero_based_idx}) "
                            f"is out of range for haplotype {i+1} (length {len(working_chain)})"
                        )

        logging.info(
            f"Assembling haplotype {i+1} from predefined chain with {len(working_chain)} repeats"
        )
        seq = assemble_haplotype_from_chain(working_chain, config)
        haplotypes.append((seq, working_chain))

    return haplotypes


def simulate_diploid(
    config: ConfigDict,
    num_haplotypes: int = 2,
    fixed_lengths: list[int] | None = None,
    seed: int | None = None,
) -> HaplotypeList:
    """Simulate multiple haplotypes with configurable parameters.

    Generates the specified number of haplotypes, either with fixed repeat lengths
    or by sampling from the configured length distribution.

    Args:
        config: Configuration dictionary with repeats, probabilities, and length model
        num_haplotypes: Number of haplotypes to simulate (default: 2 for diploid)
        fixed_lengths: Optional list of fixed repeat lengths for each haplotype
        seed: Optional random seed for reproducibility

    Returns:
        List of (haplotype_sequence, repeat_chain) tuples

    Raises:
        IndexError: If fixed_lengths is provided but has wrong length
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


def simulate_single_haplotype(
    config: ConfigDict,
    target_length: int,
    min_length: int = 10,
) -> Haplotype:
    """Build single haplotype by chaining repeats probabilistically.

    Chains repeats according to configured probabilities to reach target length.
    Enforces canonical terminal block: 6/6p -> 7 -> 8 -> 9 at the end.

    Args:
        config: Configuration with repeats, probabilities, and constants
        target_length: Desired total number of repeats
        min_length: Minimum allowed repeats (default: 10)

    Returns:
        Tuple of (assembled_sequence, repeat_chain)

    Raises:
        ValueError: If target_length is below min_length
    """
    if target_length < min_length:
        raise ValueError(
            f"Requested target_length={target_length} but minimum is {min_length}. Aborting."
        )

    reference_assembly = config.get("reference_assembly", "hg38")
    left_const = config["constants"][reference_assembly]["left"]
    right_const = config["constants"][reference_assembly]["right"]
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
