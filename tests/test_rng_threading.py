"""Tests for explicit RNG threading through domain functions."""

from __future__ import annotations

import random

from muc_one_up.distribution import sample_repeat_count
from muc_one_up.probabilities import pick_next_repeat


class TestDistributionRNG:
    def test_sample_repeat_count_accepts_rng(self):
        """sample_repeat_count should accept an rng parameter."""
        rng = random.Random(42)
        length_model = {
            "distribution": "normal",
            "min_repeats": 20,
            "max_repeats": 60,
            "mean_repeats": 40,
        }
        result = sample_repeat_count(length_model, rng=rng)
        assert 20 <= result <= 60

    def test_sample_repeat_count_reproducible_with_rng(self):
        """Same RNG seed should produce same results."""
        length_model = {
            "distribution": "normal",
            "min_repeats": 20,
            "max_repeats": 60,
            "mean_repeats": 40,
        }
        results_a = [sample_repeat_count(length_model, rng=random.Random(42)) for _ in range(5)]
        results_b = [sample_repeat_count(length_model, rng=random.Random(42)) for _ in range(5)]
        assert results_a == results_b

    def test_sample_repeat_count_defaults_to_global_random(self):
        """Without rng parameter, should still work (backward compat)."""
        length_model = {
            "distribution": "normal",
            "min_repeats": 20,
            "max_repeats": 60,
            "mean_repeats": 40,
        }
        result = sample_repeat_count(length_model)
        assert 20 <= result <= 60


class TestProbabilitiesRNG:
    def test_pick_next_repeat_accepts_rng(self):
        """pick_next_repeat should accept an rng parameter."""
        rng = random.Random(42)
        probs = {"1": {"2": 0.5, "3": 0.5}, "2": {}, "3": {}}
        result = pick_next_repeat(probs, "1", rng=rng)
        assert result in ("2", "3")

    def test_pick_next_repeat_reproducible_with_rng(self):
        """Same RNG seed should produce same sequence."""
        probs = {"1": {"2": 0.5, "3": 0.5}, "2": {}, "3": {}}
        results_a = [pick_next_repeat(probs, "1", rng=random.Random(42)) for _ in range(10)]
        results_b = [pick_next_repeat(probs, "1", rng=random.Random(42)) for _ in range(10)]
        assert results_a == results_b

    def test_pick_next_repeat_defaults_to_global_random(self):
        """Without rng parameter, should still work (backward compat)."""
        probs = {"1": {"2": 1.0}, "2": {}}
        result = pick_next_repeat(probs, "1")
        assert result == "2"
