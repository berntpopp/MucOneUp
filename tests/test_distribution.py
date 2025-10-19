"""Tests for muc_one_up.distribution module.

Tests cover:
- Sampling from normal distribution
- Respecting min/max bounds
- Handling unsupported distributions
- Reproducibility with random seed
"""

import random

import pytest

from muc_one_up.distribution import sample_repeat_count


@pytest.mark.unit
class TestSampleRepeatCount:
    """Tests for sample_repeat_count function."""

    def test_sample_normal_distribution(self):
        """Test sampling from normal distribution."""
        length_model = {
            "distribution": "normal",
            "min_repeats": 10,
            "max_repeats": 50,
            "mean_repeats": 30,
        }

        # Sample multiple times
        samples = [sample_repeat_count(length_model) for _ in range(100)]

        # All samples should be within bounds
        assert all(10 <= s <= 50 for s in samples)

        # Samples should have some variance (not all the same)
        assert len(set(samples)) > 1

    def test_sample_respects_min_bound(self):
        """Test that samples respect minimum bound."""
        length_model = {
            "distribution": "normal",
            "min_repeats": 20,
            "max_repeats": 30,
            "mean_repeats": 25,
        }

        # Sample many times to ensure bound is never violated
        samples = [sample_repeat_count(length_model) for _ in range(200)]

        assert all(s >= 20 for s in samples)

    def test_sample_respects_max_bound(self):
        """Test that samples respect maximum bound."""
        length_model = {
            "distribution": "normal",
            "min_repeats": 20,
            "max_repeats": 30,
            "mean_repeats": 25,
        }

        # Sample many times to ensure bound is never violated
        samples = [sample_repeat_count(length_model) for _ in range(200)]

        assert all(s <= 30 for s in samples)

    def test_sample_with_tight_bounds(self):
        """Test sampling with very tight bounds (min close to max)."""
        length_model = {
            "distribution": "normal",
            "min_repeats": 29,
            "max_repeats": 31,
            "mean_repeats": 30,
        }

        samples = [sample_repeat_count(length_model) for _ in range(50)]

        # All samples should be 29, 30, or 31
        assert all(s in {29, 30, 31} for s in samples)

    def test_sample_with_single_value(self):
        """Test sampling when min equals max (deterministic)."""
        length_model = {
            "distribution": "normal",
            "min_repeats": 25,
            "max_repeats": 25,
            "mean_repeats": 25,
        }

        samples = [sample_repeat_count(length_model) for _ in range(50)]

        # All samples should be exactly 25
        assert all(s == 25 for s in samples)

    def test_sample_is_reproducible_with_seed(self):
        """Test that sampling is reproducible when using random seed."""
        length_model = {
            "distribution": "normal",
            "min_repeats": 10,
            "max_repeats": 50,
            "mean_repeats": 30,
        }

        # First run with seed
        random.seed(42)
        samples1 = [sample_repeat_count(length_model) for _ in range(10)]

        # Second run with same seed
        random.seed(42)
        samples2 = [sample_repeat_count(length_model) for _ in range(10)]

        # Should be identical
        assert samples1 == samples2

    def test_sample_different_seeds_give_different_results(self):
        """Test that different seeds give different samples."""
        length_model = {
            "distribution": "normal",
            "min_repeats": 10,
            "max_repeats": 50,
            "mean_repeats": 30,
        }

        # Sample with two different seeds
        random.seed(42)
        samples1 = [sample_repeat_count(length_model) for _ in range(20)]

        random.seed(123)
        samples2 = [sample_repeat_count(length_model) for _ in range(20)]

        # Very unlikely to be identical
        assert samples1 != samples2

    def test_sample_mean_near_target(self):
        """Test that sample mean approximates target mean."""
        length_model = {
            "distribution": "normal",
            "min_repeats": 10,
            "max_repeats": 50,
            "mean_repeats": 30,
        }

        # Sample many times
        random.seed(42)
        samples = [sample_repeat_count(length_model) for _ in range(1000)]

        # Average should be close to 30 (within 10% tolerance)
        avg = sum(samples) / len(samples)
        assert 27 <= avg <= 33

    def test_sample_unsupported_distribution_returns_mean(self):
        """Test that unsupported distribution type returns mean."""
        length_model = {
            "distribution": "uniform",  # Not implemented
            "min_repeats": 10,
            "max_repeats": 50,
            "mean_repeats": 30,
        }

        result = sample_repeat_count(length_model)

        # Should return mean_repeats
        assert result == 30

    def test_sample_unknown_distribution_returns_mean(self):
        """Test that unknown distribution type returns mean."""
        length_model = {
            "distribution": "poisson",  # Not implemented
            "min_repeats": 10,
            "max_repeats": 50,
            "mean_repeats": 35,
        }

        result = sample_repeat_count(length_model)

        # Should return mean_repeats
        assert result == 35

    def test_sample_default_distribution_is_normal(self):
        """Test that omitted distribution defaults to normal."""
        length_model = {
            # No 'distribution' key
            "min_repeats": 10,
            "max_repeats": 50,
            "mean_repeats": 30,
        }

        # Should use normal distribution (default)
        samples = [sample_repeat_count(length_model) for _ in range(100)]

        # Should have variance (not deterministic)
        assert len(set(samples)) > 1
        assert all(10 <= s <= 50 for s in samples)


@pytest.mark.bioinformatics
class TestDistributionBioinformatics:
    """Bioinformatics-specific tests for repeat count sampling."""

    def test_sample_realistic_vntr_lengths(self):
        """Test sampling realistic MUC1 VNTR repeat counts."""
        # Typical MUC1 VNTR: 20-120 repeats, mean ~60
        length_model = {
            "distribution": "normal",
            "min_repeats": 20,
            "max_repeats": 120,
            "mean_repeats": 60,
        }

        random.seed(42)
        samples = [sample_repeat_count(length_model) for _ in range(100)]

        # Check biologically realistic properties
        assert all(20 <= s <= 120 for s in samples)

        # Mean should be near 60
        avg = sum(samples) / len(samples)
        assert 55 <= avg <= 65

        # Should have some diversity
        assert len(set(samples)) >= 20  # At least 20 different values

    def test_sample_short_vntr_range(self):
        """Test sampling for short VNTR regions."""
        length_model = {
            "distribution": "normal",
            "min_repeats": 5,
            "max_repeats": 15,
            "mean_repeats": 10,
        }

        samples = [sample_repeat_count(length_model) for _ in range(100)]

        assert all(5 <= s <= 15 for s in samples)
        assert len(set(samples)) > 1

    def test_sample_long_vntr_range(self):
        """Test sampling for long VNTR regions."""
        length_model = {
            "distribution": "normal",
            "min_repeats": 100,
            "max_repeats": 200,
            "mean_repeats": 150,
        }

        samples = [sample_repeat_count(length_model) for _ in range(100)]

        assert all(100 <= s <= 200 for s in samples)

        # Should still have reasonable spread
        assert max(samples) - min(samples) > 20


@pytest.mark.unit
class TestDistributionEdgeCases:
    """Edge case tests for distribution sampling."""

    def test_sample_with_mean_below_range(self):
        """Test sampling when mean is below the min-max range."""
        length_model = {
            "distribution": "normal",
            "min_repeats": 50,
            "max_repeats": 100,
            "mean_repeats": 30,  # Below min_repeats
        }

        # Should still sample within bounds
        samples = [sample_repeat_count(length_model) for _ in range(50)]

        assert all(50 <= s <= 100 for s in samples)

    def test_sample_with_mean_above_range(self):
        """Test sampling when mean is above the min-max range."""
        length_model = {
            "distribution": "normal",
            "min_repeats": 10,
            "max_repeats": 30,
            "mean_repeats": 50,  # Above max_repeats
        }

        # Should still sample within bounds
        samples = [sample_repeat_count(length_model) for _ in range(50)]

        assert all(10 <= s <= 30 for s in samples)

    def test_sample_with_large_range(self):
        """Test sampling with very large min-max range."""
        length_model = {
            "distribution": "normal",
            "min_repeats": 1,
            "max_repeats": 1000,
            "mean_repeats": 500,
        }

        samples = [sample_repeat_count(length_model) for _ in range(100)]

        assert all(1 <= s <= 1000 for s in samples)

        # Should have substantial spread
        assert max(samples) - min(samples) > 100
