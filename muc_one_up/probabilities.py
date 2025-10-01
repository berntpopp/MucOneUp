# muc_one_up/probabilities.py

import logging
import random

from .type_defs import ProbabilitiesDict


def pick_next_repeat(
    probabilities: ProbabilitiesDict,
    current_symbol: str,
    force_end: bool = False,
) -> str:
    """Randomly pick the next symbol from the probability table.

    If force_end=True, bias the pick in favor of 'END' if available.

    Args:
        probabilities: Probability table mapping current state to next state probabilities
        current_symbol: The current repeat symbol
        force_end: Boolean flag to force end state

    Returns:
        The next symbol as a string

    Example:
        >>> probs = {"1": {"2": 0.7, "7": 0.3}, "2": {"7": 1.0}}
        >>> next_sym = pick_next_repeat(probs, "1")
        >>> next_sym in ["2", "7"]
        True
    """
    if current_symbol not in probabilities:
        logging.debug("Current symbol '%s' not in probabilities; returning 'END'.", current_symbol)
        return "END"

    next_options = probabilities[current_symbol]

    if force_end and "END" in next_options:
        logging.debug("Force_end=True and 'END' available; returning 'END'.")
        return "END"

    items = list(next_options.items())
    symbols, weights = zip(*items, strict=True)
    chosen = random.choices(symbols, weights=weights, k=1)[0]
    logging.debug("Picked next symbol '%s' from current symbol '%s'.", chosen, current_symbol)
    return str(chosen)  # Type assertion: random.choices returns list of keys
