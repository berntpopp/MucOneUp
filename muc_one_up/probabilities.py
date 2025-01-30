# muc_one_up/probabilities.py
import random

def pick_next_repeat(probabilities, current_symbol, force_end=False):
    """
    Given a probability table (like config["probabilities"]) and the current_symbol,
    randomly picks the next symbol.
    
    If force_end=True, we bias the pick in favor of 'END' if available,
    or any path that leads quickly to 'END'.
    """
    if current_symbol not in probabilities:
        # fallback: if it's not in the table, just force 'END'
        return "END"

    next_options = probabilities[current_symbol]

    # If we must force an end but there's no direct 'END', consider rewriting logic or
    # picking something that leads quickly to 'END'. For now we do a simple approach:
    if force_end and "END" in next_options:
        return "END"

    # Weighted random choice among next_options
    # e.g. { "X": 0.8, "5": 0.2 }
    items = list(next_options.items())  # [(symbol, prob), ...]
    symbols, weights = zip(*items)
    return random.choices(symbols, weights=weights, k=1)[0]
