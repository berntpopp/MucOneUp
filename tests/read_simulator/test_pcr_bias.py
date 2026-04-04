"""Tests for PCR length bias model."""

import pytest

from muc_one_up.read_simulator.pcr_bias import PCRBiasModel


class TestPCRBiasModelDeterministic:
    """Tests for deterministic PCR bias computation."""

    def test_equal_length_alleles_get_equal_coverage(self):
        model = PCRBiasModel.from_preset("default")
        n1, n2 = model.compute_coverage_split(1000, 3000, 3000)
        assert n1 == 500
        assert n2 == 500

    def test_shorter_allele_gets_more_reads(self):
        model = PCRBiasModel.from_preset("default")
        n1, n2 = model.compute_coverage_split(1000, 2000, 5000)
        assert n1 > n2
        assert n1 + n2 == 1000

    def test_total_coverage_preserved(self):
        model = PCRBiasModel.from_preset("default")
        n1, n2 = model.compute_coverage_split(1000, 1500, 6000)
        assert n1 + n2 == 1000

    def test_no_bias_preset_always_equal(self):
        model = PCRBiasModel.from_preset("no_bias")
        n1, n2 = model.compute_coverage_split(1000, 1500, 6000)
        assert n1 == 500
        assert n2 == 500

    def test_symmetric_with_swapped_alleles(self):
        model = PCRBiasModel.from_preset("default")
        n1_a, n2_a = model.compute_coverage_split(1000, 2000, 4000)
        n1_b, n2_b = model.compute_coverage_split(1000, 4000, 2000)
        assert n1_a == n2_b
        assert n2_a == n1_b

    def test_odd_total_coverage_rounds_correctly(self):
        model = PCRBiasModel.from_preset("default")
        n1, n2 = model.compute_coverage_split(101, 2000, 4000)
        assert n1 + n2 == 101


class TestPCRBiasModelPresets:
    """Tests for preset loading."""

    def test_default_preset_loads(self):
        model = PCRBiasModel.from_preset("default")
        assert model.e_max == pytest.approx(0.95)
        assert model.cycles == 25

    def test_no_bias_preset_loads(self):
        model = PCRBiasModel.from_preset("no_bias")
        assert model.e_max == pytest.approx(1.0)
        assert model.alpha == pytest.approx(0.0)

    def test_invalid_preset_raises(self):
        with pytest.raises(ValueError, match="Unknown PCR bias preset"):
            PCRBiasModel.from_preset("nonexistent")

    def test_preset_with_overrides(self):
        model = PCRBiasModel.from_preset("default", cycles=30)
        assert model.cycles == 30
        assert model.e_max == pytest.approx(0.95)


class TestPCRBiasModelFromParams:
    """Tests for explicit parameter construction."""

    def test_from_params(self):
        model = PCRBiasModel.from_params(
            e_max=0.90,
            alpha=0.0002,
            cycles=20,
        )
        assert model.e_max == pytest.approx(0.90)
        assert model.alpha == pytest.approx(0.0002)
        assert model.cycles == 20

    def test_stochastic_flag(self):
        model = PCRBiasModel.from_params(
            e_max=0.95,
            alpha=0.0001,
            cycles=25,
            stochastic=True,
        )
        assert model.stochastic is True


class TestPCRBiasModelStochastic:
    """Tests for stochastic PCR bias mode."""

    def test_stochastic_reproducible_with_seed(self):
        model = PCRBiasModel.from_params(
            e_max=0.95,
            alpha=0.0001,
            cycles=25,
            stochastic=True,
        )
        r1 = model.compute_coverage_split(1000, 2000, 5000, seed=42)
        r2 = model.compute_coverage_split(1000, 2000, 5000, seed=42)
        assert r1 == r2

    def test_stochastic_different_seeds_differ(self):
        model = PCRBiasModel.from_params(
            e_max=0.95,
            alpha=0.0001,
            cycles=25,
            stochastic=True,
        )
        r1 = model.compute_coverage_split(1000, 2000, 5000, seed=42)
        r2 = model.compute_coverage_split(1000, 2000, 5000, seed=99)
        assert r1 != r2

    def test_stochastic_preserves_total(self):
        model = PCRBiasModel.from_params(
            e_max=0.95,
            alpha=0.0001,
            cycles=25,
            stochastic=True,
        )
        n1, n2 = model.compute_coverage_split(1000, 2000, 5000, seed=42)
        assert n1 + n2 == 1000

    def test_stochastic_shorter_allele_tends_higher(self):
        """Over many runs, shorter allele should average more reads."""
        model = PCRBiasModel.from_params(
            e_max=0.95,
            alpha=0.0003,
            cycles=25,
            stochastic=True,
        )
        ratios = []
        for seed in range(50):
            n1, n2 = model.compute_coverage_split(1000, 2000, 5000, seed=seed)
            ratios.append(n1 / (n1 + n2))
        avg_ratio = sum(ratios) / len(ratios)
        assert avg_ratio > 0.5


class TestPCRBiasModelFromConfig:
    """Tests for constructing from config dict."""

    def test_from_config_with_preset(self):
        config = {"preset": "default"}
        model = PCRBiasModel.from_config(config)
        assert model.e_max == pytest.approx(0.95)

    def test_from_config_with_preset_and_overrides(self):
        config = {"preset": "default", "cycles": 30}
        model = PCRBiasModel.from_config(config)
        assert model.cycles == 30

    def test_from_config_without_preset(self):
        config = {"e_max": 0.90, "alpha": 0.0002, "cycles": 20}
        model = PCRBiasModel.from_config(config)
        assert model.e_max == pytest.approx(0.90)

    def test_from_config_empty_uses_default_preset(self):
        config = {}
        model = PCRBiasModel.from_config(config)
        assert model.e_max == pytest.approx(0.95)
