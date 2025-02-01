# muc_one_up/probabilities.py

import random
import logging


def pick_next_repeat(probabilities, current_symbol, force_end=False):
    """
    Randomly pick the next symbol from the probability table for the current symbol.

    If force_end=True, bias the pick in favor of 'END' if available.

    :param probabilities: Dict representing the probability table.
    :param current_symbol: The current repeat symbol.
    :param force_end: Boolean flag to force end.
    :return: The next symbol as a string.
    """
    if current_symbol not in probabilities:
        logging.debug("Current symbol '%s' not in probabilities; returning 'END'.", current_symbol)
        return "END"

    next_options = probabilities[current_symbol]

    if force_end and "END" in next_options:
        logging.debug("Force_end=True and 'END' available; returning 'END'.")
        return "END"

    items = list(next_options.items())
    symbols, weights = zip(*items)
    chosen = random.choices(symbols, weights=weights, k=1)[0]
    logging.debug("Picked next symbol '%s' from current symbol '%s'.", chosen, current_symbol)
    return chosen
