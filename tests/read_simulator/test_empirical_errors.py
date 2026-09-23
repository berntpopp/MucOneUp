"""Tests for the empirical long-read error channel."""

from __future__ import annotations

import random

import pytest

from muc_one_up.read_simulator.empirical_errors import EmpiricalErrorModel, apply_errors

_SEQ_RNG = random.Random(0)
RNG_SEQ = "".join(_SEQ_RNG.choice("ACGT") for _ in range(20000))
MODEL = EmpiricalErrorModel(
    mismatch_rate=0.006,
    insertion_rate=0.004,
    deletion_rate=0.003,
    insertion_len_pmf={1: 0.85, 2: 0.15},
    deletion_len_pmf={1: 0.8, 2: 0.2},
    read_error_sigma=0.0,
)


def test_zero_rates_are_identity() -> None:
    clean = EmpiricalErrorModel(0.0, 0.0, 0.0, {1: 1.0}, {1: 1.0}, read_error_sigma=0.0)
    read, qual = apply_errors("ACGTTTTGCA", clean, random.Random(1))
    assert read == "ACGTTTTGCA" and len(qual) == len(read)


def test_rates_match_configuration() -> None:
    from collections import Counter

    events: Counter[str] = Counter()
    for i in range(10):
        apply_errors(RNG_SEQ, MODEL, random.Random(i), events)
    exposed = events["exposed"]
    assert exposed < 10 * len(RNG_SEQ)  # homopolymer runs are excluded
    for name, rate in (("mismatch", 0.006), ("insertion", 0.004), ("deletion", 0.003)):
        assert events[name] / exposed == pytest.approx(rate, rel=0.15)


def test_homopolymer_runs_are_untouched() -> None:
    seq = ("ACGT" + "C" * 7) * 500
    noisy = EmpiricalErrorModel(0.05, 0.05, 0.05, {1: 1.0}, {1: 1.0}, read_error_sigma=0.0)
    read, _ = apply_errors(seq, noisy, random.Random(2))
    assert read.count("C" * 7) >= 450  # runs themselves never receive channel errors


def test_quality_tracks_read_error() -> None:
    _, good = apply_errors(RNG_SEQ, MODEL, random.Random(3))
    worse = EmpiricalErrorModel(0.03, 0.02, 0.02, {1: 1.0}, {1: 1.0}, read_error_sigma=0.0)
    _, bad = apply_errors(RNG_SEQ, worse, random.Random(3))
    mean_q = lambda q: sum(ord(c) - 33 for c in q) / len(q)  # noqa: E731
    assert mean_q(good) > mean_q(bad) + 3


def test_deterministic_per_seed() -> None:
    assert apply_errors(RNG_SEQ, MODEL, random.Random(9)) == apply_errors(
        RNG_SEQ, MODEL, random.Random(9)
    )


@pytest.mark.parametrize(
    "data,match",
    [
        ({"mismatch_rate": 0.5}, "rate"),
        ({"insertion_len_pmf": {"1": 0.5}}, "sum"),
        ({"surprise": 1}, "unknown"),
    ],
)
def test_from_dict_validation(data: dict, match: str) -> None:
    base = {
        "mismatch_rate": 0.006,
        "insertion_rate": 0.004,
        "deletion_rate": 0.003,
        "insertion_len_pmf": {"1": 1.0},
        "deletion_len_pmf": {"1": 1.0},
    }
    with pytest.raises(ValueError, match=match):
        EmpiricalErrorModel.from_dict({**base, **data})
