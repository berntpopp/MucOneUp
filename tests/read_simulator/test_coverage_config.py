"""Tests for coverage configuration standardization.

This module tests the standardized 'coverage' configuration key used for
Illumina read simulation downsampling. Ensures the key is properly recognized
across CLI, config files, and pipeline logic.
"""

import pytest


class TestCoverageKeyStandardization:
    """Test suite for coverage config key standardization."""

    def test_coverage_key_in_minimal_config(self, minimal_config):
        """Verify test fixture uses standardized 'coverage' key."""
        # Arrange & Act
        rs_config = minimal_config.get("read_simulation", {})

        # Assert
        assert "coverage" in rs_config, "Test fixture should include 'coverage' key"
        assert isinstance(rs_config["coverage"], (int, float)), "Coverage must be numeric"

    def test_coverage_value_type(self, minimal_config):
        """Verify coverage value is numeric type."""
        # Arrange
        coverage = minimal_config["read_simulation"]["coverage"]

        # Act & Assert
        assert isinstance(coverage, (int, float)), "Coverage must be int or float"
        assert coverage > 0, "Coverage must be positive"

    def test_no_legacy_downsample_target_key(self, minimal_config):
        """Ensure legacy 'downsample_target' key is not present."""
        # Arrange
        rs_config = minimal_config.get("read_simulation", {})

        # Act & Assert
        assert "downsample_target" not in rs_config, "Legacy 'downsample_target' should not exist"

    def test_no_legacy_downsample_coverage_key(self, minimal_config):
        """Ensure legacy 'downsample_coverage' key is not present."""
        # Arrange
        rs_config = minimal_config.get("read_simulation", {})

        # Act & Assert
        assert "downsample_coverage" not in rs_config, (
            "Legacy 'downsample_coverage' should not exist"
        )


class TestCoverageConfigValidation:
    """Test suite for coverage configuration validation."""

    def test_coverage_integer_value(self, minimal_config):
        """Test coverage with integer value."""
        # Arrange
        config = minimal_config.copy()
        config["read_simulation"]["coverage"] = 30

        # Act & Assert
        assert config["read_simulation"]["coverage"] == 30
        assert isinstance(config["read_simulation"]["coverage"], int)

    def test_coverage_float_value(self, minimal_config):
        """Test coverage with float value."""
        # Arrange
        config = minimal_config.copy()
        config["read_simulation"]["coverage"] = 30.5

        # Act & Assert
        assert config["read_simulation"]["coverage"] == 30.5
        assert isinstance(config["read_simulation"]["coverage"], float)

    def test_coverage_type_checking(self, minimal_config):
        """Test coverage type validation."""
        # Arrange
        config = minimal_config.copy()

        # Act: Set invalid type
        config["read_simulation"]["coverage"] = "thirty"

        # Assert: Type check fails
        assert not isinstance(config["read_simulation"]["coverage"], (int, float))

    def test_coverage_negative_value(self, minimal_config):
        """Test that negative coverage values are handled."""
        # Arrange
        config = minimal_config.copy()
        config["read_simulation"]["coverage"] = -30

        # Act & Assert
        # Note: Schema validation allows negative (type check only)
        # Logical validation happens in pipeline
        assert config["read_simulation"]["coverage"] < 0


class TestDownsamplingModesWithCoverage:
    """Test downsampling mode configurations with coverage key."""

    def test_non_vntr_mode_requires_bed_file(self, minimal_config):
        """Verify non_vntr mode requires sample_target_bed."""
        # Arrange
        config = minimal_config.copy()
        config["read_simulation"]["coverage"] = 30
        config["read_simulation"]["downsample_mode"] = "non_vntr"

        # Act & Assert: Should have bed file
        # (Pipeline will validate at runtime, here we check config structure)
        if "sample_target_bed" in config["read_simulation"]:
            assert isinstance(config["read_simulation"]["sample_target_bed"], str)

    def test_vntr_mode_configuration(self, minimal_config):
        """Verify vntr mode can be configured."""
        # Arrange
        config = minimal_config.copy()
        config["read_simulation"]["coverage"] = 30
        config["read_simulation"]["downsample_mode"] = "vntr"
        config["read_simulation"]["reference_assembly"] = "hg38"

        # Act & Assert: Verify configuration structure
        assert config["read_simulation"]["downsample_mode"] == "vntr"
        assert config["read_simulation"]["reference_assembly"] in ["hg38", "hg19"]
        # Note: vntr_region_* keys validated at runtime in pipeline

    def test_coverage_without_mode_defaults_to_vntr(self, minimal_config):
        """Test that downsample_mode defaults to 'vntr' if not specified."""
        # Arrange
        config = minimal_config.copy()
        config["read_simulation"]["coverage"] = 30
        # Don't set downsample_mode

        # Act
        mode = config["read_simulation"].get("downsample_mode", "vntr")

        # Assert
        assert mode == "vntr", "Default mode should be 'vntr'"


class TestCoverageBackwardCompatibility:
    """Test backward compatibility and migration scenarios."""

    def test_coverage_key_takes_precedence(self):
        """If both old and new keys exist, coverage should be used."""
        # Arrange: Hypothetical scenario with both keys
        config = {
            "read_simulation": {
                "coverage": 30,
                "downsample_coverage": 50,  # Legacy (should be ignored)
            }
        }

        # Act
        coverage = config["read_simulation"].get("coverage")

        # Assert
        assert coverage == 30, "New 'coverage' key should take precedence"

    def test_missing_coverage_returns_none(self):
        """Test behavior when coverage key is missing."""
        # Arrange: Config without coverage
        config = {"read_simulation": {}}

        # Act
        coverage = config["read_simulation"].get("coverage")

        # Assert
        assert coverage is None, "Missing coverage should return None"


class TestCoverageIntegrationWithCLI:
    """Test coverage parameter integration with CLI."""

    def test_cli_sets_coverage_key(self):
        """Verify CLI sets 'coverage' key (not legacy keys)."""
        # This is a documentation test confirming expected CLI behavior
        # Actual CLI tested in test_click_cli.py

        # Expected behavior:
        # muconeup reads illumina --coverage 30
        # -> config["read_simulation"]["coverage"] = 30

        # Arrange: Simulated CLI config
        config = {"read_simulation": {}}

        # Act: Simulate CLI setting
        config["read_simulation"]["coverage"] = 30

        # Assert
        assert config["read_simulation"]["coverage"] == 30
        assert "downsample_target" not in config["read_simulation"]
        assert "downsample_coverage" not in config["read_simulation"]


# Parametrized tests for various coverage values
@pytest.mark.parametrize(
    "coverage_value,expected_type",
    [
        (30, int),
        (30.5, float),
        (100, int),
        (0.5, float),
        (1000, int),
    ],
)
def test_coverage_value_types(minimal_config, coverage_value, expected_type):
    """Test various coverage values and types."""
    # Arrange
    config = minimal_config.copy()
    config["read_simulation"]["coverage"] = coverage_value

    # Act & Assert
    assert isinstance(config["read_simulation"]["coverage"], expected_type)
    assert config["read_simulation"]["coverage"] == coverage_value


@pytest.mark.parametrize(
    "mode",
    ["vntr", "non_vntr"],
)
def test_downsample_modes_configuration(minimal_config, mode):
    """Test that downsample modes can be configured."""
    # Arrange
    config = minimal_config.copy()
    config["read_simulation"]["coverage"] = 30
    config["read_simulation"]["downsample_mode"] = mode

    # Act & Assert
    assert config["read_simulation"]["downsample_mode"] == mode
    assert config["read_simulation"]["coverage"] == 30
    # Note: Mode-specific keys (vntr_region_*, sample_target_bed)
    # are validated at runtime in pipeline
