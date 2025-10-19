import random

from muc_one_up.probabilities import pick_next_repeat


def test_pick_next_repeat_simple():
    """Test pick_next_repeat with a simple probability dictionary."""
    probabilities = {
        "A": {"B": 0.5, "C": 0.5},
    }
    # Mock random so we can control the outcome
    random.seed(1234)
    next_sym = pick_next_repeat(probabilities, "A")
    # Because of the seed, you get the same random choice
    # It's either "B" or "C". Let's see if we just do an assert in the set:
    assert next_sym in ("B", "C")


def test_pick_next_repeat_force_end():
    """If force_end=True and 'END' is available, we pick it."""
    probabilities = {"X": {"END": 0.9, "A": 0.1}}
    next_sym = pick_next_repeat(probabilities, "X", force_end=True)
    assert next_sym == "END"


def test_pick_next_repeat_missing_symbol():
    """If current_symbol not in probabilities, fallback is 'END'."""
    probabilities = {"A": {"B": 1.0}}
    next_sym = pick_next_repeat(probabilities, "Z")
    assert next_sym == "END"
