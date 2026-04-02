"""Tests for explicit RNG threading through domain functions."""

from __future__ import annotations

import random

from muc_one_up.distribution import sample_repeat_count
from muc_one_up.probabilities import pick_next_repeat
from muc_one_up.simulate import simulate_diploid


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
        """Same RNG seed should produce same sequence of results."""
        length_model = {
            "distribution": "normal",
            "min_repeats": 20,
            "max_repeats": 60,
            "mean_repeats": 40,
        }
        rng_a = random.Random(42)
        results_a = [sample_repeat_count(length_model, rng=rng_a) for _ in range(5)]
        rng_b = random.Random(42)
        results_b = [sample_repeat_count(length_model, rng=rng_b) for _ in range(5)]
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
        """Same RNG seed should produce same sequence of results."""
        probs = {"1": {"2": 0.5, "3": 0.5}, "2": {}, "3": {}}
        rng_a = random.Random(42)
        results_a = [pick_next_repeat(probs, "1", rng=rng_a) for _ in range(10)]
        rng_b = random.Random(42)
        results_b = [pick_next_repeat(probs, "1", rng=rng_b) for _ in range(10)]
        assert results_a == results_b

    def test_pick_next_repeat_defaults_to_global_random(self):
        """Without rng parameter, should still work (backward compat)."""
        probs = {"1": {"2": 1.0}, "2": {}}
        result = pick_next_repeat(probs, "1")
        assert result == "2"


class TestSimulateRNG:
    def test_simulate_diploid_accepts_rng(self):
        """simulate_diploid should accept an rng parameter."""
        import inspect

        sig = inspect.signature(simulate_diploid)
        assert "rng" in sig.parameters

    def test_simulate_diploid_reproducible_with_rng(self):
        """Same RNG produces identical haplotypes."""
        config = {
            "repeats": {
                "1": "GCCACCACCTCCAACTCCT",
                "2": "GCCACCACC",
                "3": "TCCAACTCCT",
                "6": "TCCGACTCCT",
                "6p": "TCCAACTCCT",
                "7": "TCTAGGACCTCCTAGCTCCCAGAACTCCT",
                "8": "TCTAGGACCTCCTAGCTCCTAGCTCCT",
                "9": "TCTAGGACCTAGCTCCT",
            },
            "constants": {
                "hg38": {
                    "left": "ATGGCCCCATCTCTCACCGTCTCGGTCATCTCCTTGATG",
                    "right": "CAATGGTGTCTTGGGTAGCTTCGTCACGGTTTTCCAG",
                }
            },
            "probabilities": {
                "1": {"2": 1.0},
                "2": {"3": 1.0},
                "3": {"6": 0.5, "6p": 0.5},
                "6": {"7": 1.0},
                "6p": {"7": 1.0},
                "7": {"8": 1.0},
                "8": {"9": 1.0},
                "9": {},
            },
            "length_model": {
                "distribution": "normal",
                "min_repeats": 20,
                "max_repeats": 40,
                "mean_repeats": 30,
            },
            "reference_assembly": "hg38",
        }
        rng1 = random.Random(42)
        rng2 = random.Random(42)
        r1 = simulate_diploid(config, seed=None, rng=rng1)
        r2 = simulate_diploid(config, seed=None, rng=rng2)
        assert r1[0].chain == r2[0].chain  # same chains
        assert r1[1].chain == r2[1].chain
