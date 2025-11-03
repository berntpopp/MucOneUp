"""Tests for muc_one_up.config_fingerprint module.

Tests cover:
- Config canonicalization (granular filtering)
- Config fingerprinting (RFC 8785 compliance)
- Error handling (non-serializable objects)
- Edge cases (circular refs, deep nesting, unicode, etc.)
"""

import pytest

from muc_one_up.config_fingerprint import (
    PATH_KEYS,
    canonicalize_config,
    compute_config_fingerprint,
)


@pytest.mark.unit
class TestCanonicalizeConfig:
    """Tests for canonicalize_config function."""

    def test_removes_tools_section_entirely(self):
        """Test that entire 'tools' section is removed."""
        config = {
            "repeats": {"X": "ACGACT"},
            "tools": {"bwa": "/usr/bin/bwa", "samtools": "/usr/bin/samtools"},
        }
        canonical = canonicalize_config(config)
        assert "repeats" in canonical
        assert "tools" not in canonical

    def test_removes_path_keys_from_nested_dicts(self):
        """Test that path keys are removed from nested dicts."""
        config = {
            "read_simulation": {
                "coverage": 100,
                "human_reference": "/path/to/hg38.fa",  # Should be excluded
                "threads": 4,
            }
        }
        canonical = canonicalize_config(config)
        assert "read_simulation" in canonical
        assert canonical["read_simulation"]["coverage"] == 100
        assert canonical["read_simulation"]["threads"] == 4
        assert "human_reference" not in canonical["read_simulation"]

    def test_preserves_scientific_parameters(self):
        """Test that scientific params are preserved in fingerprint."""
        config = {
            "seed": 42,
            "read_simulation": {"coverage": 100, "threads": 4},
            "nanosim_params": {
                "coverage": 50,
                "min_read_length": 1500,
                "training_data_path": "/path/to/model",  # Excluded
            },
        }
        canonical = canonicalize_config(config)

        # Scientific params preserved
        assert canonical["seed"] == 42
        assert canonical["read_simulation"]["coverage"] == 100
        assert canonical["read_simulation"]["threads"] == 4
        assert canonical["nanosim_params"]["coverage"] == 50
        assert canonical["nanosim_params"]["min_read_length"] == 1500

        # Paths excluded
        assert "training_data_path" not in canonical["nanosim_params"]

    def test_handles_deeply_nested_structures(self):
        """Test filtering in deeply nested structures."""
        config = {
            "outer": {
                "middle": {
                    "inner": {"value": 123, "human_reference": "/path"},  # Path excluded
                    "keep": "this",
                }
            }
        }
        canonical = canonicalize_config(config)
        assert canonical["outer"]["middle"]["inner"]["value"] == 123
        assert "human_reference" not in canonical["outer"]["middle"]["inner"]
        assert canonical["outer"]["middle"]["keep"] == "this"

    def test_handles_lists_in_config(self):
        """Test that lists are preserved."""
        config = {
            "repeats": {"X": "ACGACT"},
            "targets": [1, 2, 3],
            "nested_list": [{"a": 1}, {"b": 2}],
        }
        canonical = canonicalize_config(config)
        assert canonical["targets"] == [1, 2, 3]
        assert canonical["nested_list"] == [{"a": 1}, {"b": 2}]

    def test_handles_empty_config(self):
        """Test with empty config."""
        canonical = canonicalize_config({})
        assert canonical == {}

    def test_handles_none_values(self):
        """Test config with None values."""
        config = {"seed": None, "repeats": {"X": "ACGACT"}}
        canonical = canonicalize_config(config)
        assert canonical["seed"] is None
        assert "repeats" in canonical

    def test_circular_reference_raises_error(self):
        """Test handling of circular references."""
        config = {"self": None}
        config["self"] = config  # Circular reference
        with pytest.raises(RecursionError):
            canonicalize_config(config)


@pytest.mark.unit
class TestComputeConfigFingerprint:
    """Tests for compute_config_fingerprint function."""

    def test_deterministic(self):
        """Test same config produces same fingerprint."""
        config = {"repeats": {"X": "ACGACT"}, "seed": 42}
        fp1 = compute_config_fingerprint(config)
        fp2 = compute_config_fingerprint(config)
        assert fp1 == fp2

    def test_unique_configs_different_fingerprints(self):
        """Test different configs produce different fingerprints."""
        config1 = {"repeats": {"X": "ACGACT"}, "seed": 42}
        config2 = {"repeats": {"X": "ACGACT"}, "seed": 43}
        fp1 = compute_config_fingerprint(config1)
        fp2 = compute_config_fingerprint(config2)
        assert fp1 != fp2

    def test_fingerprint_format(self):
        """Test fingerprint format."""
        config = {"repeats": {"X": "ACGACT"}}
        fp = compute_config_fingerprint(config)
        assert fp.startswith("sha256:")
        assert len(fp) == 71  # "sha256:" (7) + 64 hex chars

    def test_key_order_independence(self):
        """Test that key order doesn't affect fingerprint (RFC 8785)."""
        config1 = {"a": 1, "b": 2, "c": 3}
        config2 = {"c": 3, "a": 1, "b": 2}
        fp1 = compute_config_fingerprint(config1)
        fp2 = compute_config_fingerprint(config2)
        assert fp1 == fp2  # RFC 8785 sorts keys

    def test_nested_structures(self):
        """Test with nested config."""
        config = {
            "repeats": {"X": "ACGACT", "A": "TCGACT"},
            "probabilities": {"X": {"A": 0.5, "B": 0.5}},
            "seed": 42,
        }
        fp = compute_config_fingerprint(config)
        assert fp.startswith("sha256:")

    def test_different_coverage_different_fingerprint(self):
        """Test different coverage produces different fingerprint."""
        config1 = {"read_simulation": {"coverage": 30}}
        config2 = {"read_simulation": {"coverage": 100}}
        fp1 = compute_config_fingerprint(config1)
        fp2 = compute_config_fingerprint(config2)
        assert fp1 != fp2, "Different coverage must produce different fingerprints!"

    def test_same_coverage_different_paths_same_fingerprint(self):
        """Test that different paths don't affect fingerprint (filtered)."""
        config1 = {"read_simulation": {"coverage": 100, "human_reference": "/path/1/hg38.fa"}}
        config2 = {"read_simulation": {"coverage": 100, "human_reference": "/path/2/hg38.fa"}}
        fp1 = compute_config_fingerprint(config1)
        fp2 = compute_config_fingerprint(config2)
        # Should be identical (paths excluded)
        assert fp1 == fp2

    def test_handles_non_serializable_config(self):
        """Test graceful handling of non-JSON-serializable config."""
        config = {"data": object()}  # Non-serializable
        fp = compute_config_fingerprint(config)
        assert fp.startswith("error:fingerprint_failed")

    def test_handles_circular_reference_gracefully(self):
        """Test graceful handling of circular references."""
        config = {"self": None}
        config["self"] = config  # Circular reference
        fp = compute_config_fingerprint(config)
        assert fp.startswith("error:fingerprint_failed")

    def test_handles_infinity(self):
        """Test config with infinity values."""
        config = {"value": float("inf")}
        fp = compute_config_fingerprint(config)
        # rfc8785 may handle this or raise - either is acceptable
        assert fp.startswith("sha256:") or fp.startswith("error:")

    def test_handles_nan(self):
        """Test config with NaN values."""
        config = {"value": float("nan")}
        fp = compute_config_fingerprint(config)
        # rfc8785 may handle this or raise - either is acceptable
        assert fp.startswith("sha256:") or fp.startswith("error:")

    def test_handles_unicode(self):
        """Test unicode characters are handled correctly."""
        config = {"description": "MUC1 VNTR — μm resolution"}
        fp = compute_config_fingerprint(config)
        assert fp.startswith("sha256:")

    def test_handles_empty_nested_dicts(self):
        """Test config with empty nested dicts."""
        config1 = {"outer": {"inner": {}}}
        config2 = {"outer": {}}
        fp1 = compute_config_fingerprint(config1)
        fp2 = compute_config_fingerprint(config2)
        # Should be different (structure matters)
        assert fp1 != fp2

    def test_empty_config_produces_valid_fingerprint(self):
        """Test empty config."""
        config = {}
        fp = compute_config_fingerprint(config)
        assert fp.startswith("sha256:")

    def test_float_precision_affects_fingerprint(self):
        """Test that float precision differences are detected."""
        config1 = {"value": 0.1 + 0.2}  # 0.30000000000000004
        config2 = {"value": 0.3}
        fp1 = compute_config_fingerprint(config1)
        fp2 = compute_config_fingerprint(config2)
        # Should be different (this is correct!)
        assert fp1 != fp2


@pytest.mark.unit
class TestPathKeysConstant:
    """Tests for PATH_KEYS constant."""

    def test_path_keys_is_frozenset(self):
        """Test PATH_KEYS is immutable."""
        assert isinstance(PATH_KEYS, frozenset)

    def test_path_keys_contains_expected_keys(self):
        """Test PATH_KEYS contains expected path keys."""
        expected_keys = {
            "human_reference",
            "training_data_path",
            "model_file",
            "bwa",
            "samtools",
        }
        assert expected_keys.issubset(PATH_KEYS)

    def test_path_keys_not_empty(self):
        """Test PATH_KEYS is not empty."""
        assert len(PATH_KEYS) > 0
