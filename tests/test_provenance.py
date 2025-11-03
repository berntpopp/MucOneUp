"""Tests for muc_one_up.provenance module.

Tests cover:
- Software version retrieval
- Seed extraction from multiple config locations
- Timestamp formatting (ISO 8601)
- Command-line capture and sanitization
- Provenance metadata collection
- Error handling and graceful degradation
"""

import sys
import time
from unittest import mock

import pytest

from muc_one_up.provenance import (
    ENABLE_PROVENANCE,
    collect_provenance_metadata,
    extract_seed,
    format_timestamp_iso8601,
    get_command_line,
    get_software_version,
)
from muc_one_up.version import __version__


@pytest.mark.unit
class TestGetSoftwareVersion:
    """Tests for get_software_version function."""

    def test_returns_version_string(self):
        """Test that version is returned as string."""
        version = get_software_version()
        assert isinstance(version, str)
        assert version == __version__

    def test_semantic_versioning_format(self):
        """Test version follows semantic versioning."""
        version = get_software_version()
        parts = version.split(".")
        assert len(parts) >= 3  # major.minor.patch


@pytest.mark.unit
class TestExtractSeed:
    """Tests for extract_seed function."""

    def test_extract_from_top_level(self):
        """Test seed extraction from config root."""
        config = {"seed": 42}
        seed = extract_seed(config)
        assert seed == 42

    def test_extract_from_read_simulation(self):
        """Test seed extraction from read_simulation section."""
        config = {"read_simulation": {"seed": 123}}
        seed = extract_seed(config)
        assert seed == 123

    def test_extract_from_nanosim_params(self):
        """Test seed extraction from nanosim_params section."""
        config = {"nanosim_params": {"seed": 456}}
        seed = extract_seed(config)
        assert seed == 456

    def test_extract_from_pacbio_params(self):
        """Test seed extraction from pacbio_params section."""
        config = {"pacbio_params": {"seed": 789}}
        seed = extract_seed(config)
        assert seed == 789

    def test_priority_order_top_level_wins(self):
        """Test that top-level seed takes priority."""
        config = {
            "seed": 42,
            "read_simulation": {"seed": 123},
            "nanosim_params": {"seed": 456},
        }
        seed = extract_seed(config)
        assert seed == 42  # Top-level wins

    def test_priority_order_read_simulation_over_nanosim(self):
        """Test read_simulation takes priority over nanosim_params."""
        config = {
            "read_simulation": {"seed": 123},
            "nanosim_params": {"seed": 456},
        }
        seed = extract_seed(config)
        assert seed == 123  # read_simulation wins

    def test_missing_seed_returns_none(self):
        """Test graceful handling of missing seed."""
        config = {"repeats": {"X": "ACGACT"}}
        seed = extract_seed(config)
        assert seed is None

    def test_none_seed_value_skipped(self):
        """Test that None seed values are skipped."""
        config = {"seed": None, "read_simulation": {"seed": 123}}
        seed = extract_seed(config)
        assert seed == 123  # Skips None, finds next

    def test_handles_malformed_config(self):
        """Test handling of malformed config."""
        config = {"read_simulation": "not_a_dict"}
        seed = extract_seed(config)
        assert seed is None  # Graceful degradation

    def test_empty_config(self):
        """Test with empty config."""
        seed = extract_seed({})
        assert seed is None


@pytest.mark.unit
class TestFormatTimestampIso8601:
    """Tests for format_timestamp_iso8601 function."""

    def test_format_with_integer_timestamp(self):
        """Test formatting with integer Unix timestamp."""
        ts = 1730628930  # 2024-11-03 10:15:30 UTC (corrected timestamp)
        iso = format_timestamp_iso8601(ts)
        assert iso.startswith("2024-11-03T10:15:30")
        assert "+00:00" in iso or "Z" in iso  # UTC timezone

    def test_format_with_float_timestamp(self):
        """Test formatting with float Unix timestamp."""
        ts = 1730628930.123456  # 2024-11-03 10:15:30.123456 UTC (corrected timestamp)
        iso = format_timestamp_iso8601(ts)
        assert ".123456" in iso  # Microseconds included
        assert "2024-11-03" in iso
        assert "+00:00" in iso  # UTC timezone

    def test_microseconds_included(self):
        """Test microseconds are included in output."""
        ts = 1730628930.123456
        iso = format_timestamp_iso8601(ts)
        assert ".123456" in iso

    def test_timezone_offset_included(self):
        """Test timezone offset is included."""
        ts = 1730628930.0
        iso = format_timestamp_iso8601(ts)
        assert "+00:00" in iso  # UTC offset

    def test_handles_zero_timestamp(self):
        """Test formatting Unix epoch (0)."""
        iso = format_timestamp_iso8601(0.0)
        assert "1970-01-01" in iso

    def test_handles_recent_timestamp(self):
        """Test formatting recent timestamp."""
        ts = time.time()
        iso = format_timestamp_iso8601(ts)
        assert isinstance(iso, str)
        assert "T" in iso  # ISO 8601 format

    def test_negative_timestamp_works(self):
        """Test negative timestamps work (dates before 1970-01-01)."""
        # Negative timestamps are valid in Python 3.x
        iso = format_timestamp_iso8601(-1.0)
        assert "1969-12-31" in iso
        assert isinstance(iso, str)

    def test_very_large_timestamp_raises_error(self):
        """Test very large timestamp raises error (year 9999+)."""
        with pytest.raises((ValueError, OSError)):
            format_timestamp_iso8601(9999999999999.0)


@pytest.mark.unit
class TestGetCommandLine:
    """Tests for get_command_line function."""

    def test_returns_non_empty_string(self):
        """Test command line is non-empty."""
        cmd = get_command_line(sanitize=False)
        assert isinstance(cmd, str)
        assert len(cmd) > 0

    def test_contains_program_name(self):
        """Test command line contains program name."""
        cmd = get_command_line(sanitize=False)
        # Should contain pytest (or python)
        assert "pytest" in cmd or "python" in cmd

    def test_sanitizes_api_keys(self):
        """Test API keys are redacted from command line."""
        original_argv = sys.argv
        try:
            sys.argv = ["muconeup", "--api-key", "secret123", "simulate"]
            cmd = get_command_line(sanitize=True)
            assert "secret123" not in cmd
            assert "***REDACTED***" in cmd
            assert "--api-key" in cmd  # Flag preserved
        finally:
            sys.argv = original_argv

    def test_sanitizes_inline_secrets(self):
        """Test inline secrets are redacted."""
        original_argv = sys.argv
        try:
            sys.argv = ["muconeup", "--token=abc123", "run"]
            cmd = get_command_line(sanitize=True)
            assert "abc123" not in cmd
            assert "--token=***REDACTED***" in cmd
        finally:
            sys.argv = original_argv

    def test_sanitizes_multiple_secrets(self):
        """Test multiple secrets are redacted."""
        original_argv = sys.argv
        try:
            sys.argv = [
                "muconeup",
                "--api-key",
                "key123",
                "--password",
                "pass456",
                "run",
            ]
            cmd = get_command_line(sanitize=True)
            assert "key123" not in cmd
            assert "pass456" not in cmd
            assert cmd.count("***REDACTED***") == 2
        finally:
            sys.argv = original_argv

    def test_preserves_non_secret_args(self):
        """Test non-secret args are preserved."""
        original_argv = sys.argv
        try:
            sys.argv = ["muconeup", "--config", "config.json", "--seed", "42"]
            cmd = get_command_line(sanitize=True)
            assert "--config" in cmd
            assert "config.json" in cmd
            assert "--seed" in cmd
            assert "42" in cmd
        finally:
            sys.argv = original_argv

    def test_without_sanitization(self):
        """Test command line without sanitization."""
        original_argv = sys.argv
        try:
            sys.argv = ["muconeup", "--api-key", "secret123"]
            cmd = get_command_line(sanitize=False)
            assert "secret123" in cmd  # Not redacted
        finally:
            sys.argv = original_argv


@pytest.mark.unit
class TestCollectProvenanceMetadata:
    """Tests for collect_provenance_metadata function."""

    def test_all_fields_present(self):
        """Test all required fields are present."""
        config = {"repeats": {"X": "ACGACT"}, "seed": 42}
        start = time.time()
        end = start + 1.0

        provenance = collect_provenance_metadata(config, start, end)

        assert "software_version" in provenance
        assert "config_fingerprint" in provenance
        assert "seed" in provenance
        assert "start_time" in provenance
        assert "end_time" in provenance
        assert "duration_seconds" in provenance
        assert "command_line" in provenance

    def test_correct_values(self):
        """Test values are correct."""
        config = {"repeats": {"X": "ACGACT"}, "seed": 42}
        start = 1730649330.0
        end = 1730649331.5

        provenance = collect_provenance_metadata(config, start, end)

        assert provenance["software_version"] == __version__
        assert provenance["seed"] == 42
        assert provenance["duration_seconds"] == 1.5
        assert provenance["config_fingerprint"].startswith("sha256:")

    def test_missing_seed_handled(self):
        """Test graceful handling of missing seed."""
        config = {"repeats": {"X": "ACGACT"}}
        start = time.time()
        end = start + 1.0

        provenance = collect_provenance_metadata(config, start, end)

        assert provenance["seed"] is None
        assert "config_fingerprint" in provenance  # Other fields still work

    def test_duration_calculation(self):
        """Test duration is calculated correctly."""
        config = {}
        start = 1000.0
        end = 1005.5

        provenance = collect_provenance_metadata(config, start, end)

        assert provenance["duration_seconds"] == 5.5

    def test_negative_duration_allowed(self):
        """Test negative duration (end < start) is allowed."""
        config = {}
        start = 1000.0
        end = 999.0  # Wrong order!

        provenance = collect_provenance_metadata(config, start, end)

        assert provenance["duration_seconds"] == -1.0  # Negative but doesn't crash

    def test_start_time_is_iso8601(self):
        """Test start_time is ISO 8601 formatted."""
        config = {}
        start = 1730649330.0
        end = start + 1.0

        provenance = collect_provenance_metadata(config, start, end)

        assert "T" in provenance["start_time"]  # ISO 8601 format
        assert "+00:00" in provenance["start_time"]  # UTC timezone

    def test_end_time_is_iso8601(self):
        """Test end_time is ISO 8601 formatted."""
        config = {}
        start = 1730649330.0
        end = start + 1.0

        provenance = collect_provenance_metadata(config, start, end)

        assert "T" in provenance["end_time"]
        assert "+00:00" in provenance["end_time"]

    def test_command_line_sanitized(self):
        """Test command line is sanitized."""
        original_argv = sys.argv
        try:
            sys.argv = ["muconeup", "--api-key", "secret123", "run"]
            config = {}
            start = time.time()
            end = start + 1.0

            provenance = collect_provenance_metadata(config, start, end)

            assert "secret123" not in provenance["command_line"]
            assert "***REDACTED***" in provenance["command_line"]
        finally:
            sys.argv = original_argv

    def test_handles_config_fingerprint_error(self):
        """Test graceful handling of config fingerprint errors."""
        config = {"data": object()}  # Non-serializable
        start = time.time()
        end = start + 1.0

        provenance = collect_provenance_metadata(config, start, end)

        # Should have error sentinel, not crash
        assert provenance["config_fingerprint"].startswith("error:")
        assert "software_version" in provenance  # Other fields still work

    def test_respects_feature_flag_disabled(self):
        """Test provenance collection can be disabled via feature flag."""
        with mock.patch.dict("os.environ", {"MUCONEUP_ENABLE_PROVENANCE": "false"}):
            # Need to reload module to pick up env var
            from importlib import reload

            from muc_one_up import provenance as prov_module

            reload(prov_module)

            config = {}
            start = time.time()
            end = start + 1.0

            result = prov_module.collect_provenance_metadata(config, start, end)

            # Should return empty dict when disabled
            assert result == {}

    def test_feature_flag_enabled_by_default(self):
        """Test provenance collection enabled by default."""
        # Don't set env var, should default to enabled
        assert ENABLE_PROVENANCE is True  # Default


@pytest.mark.unit
class TestFeatureFlag:
    """Tests for ENABLE_PROVENANCE feature flag."""

    def test_feature_flag_recognizes_true_values(self):
        """Test feature flag recognizes true values."""
        for value in ["true", "1", "yes", "on", "True", "YES", "ON"]:
            with mock.patch.dict("os.environ", {"MUCONEUP_ENABLE_PROVENANCE": value}):
                from importlib import reload

                from muc_one_up import provenance as prov_module

                reload(prov_module)
                assert prov_module.ENABLE_PROVENANCE is True

    def test_feature_flag_recognizes_false_values(self):
        """Test feature flag recognizes false values."""
        for value in ["false", "0", "no", "off", "False", "NO", "OFF"]:
            with mock.patch.dict("os.environ", {"MUCONEUP_ENABLE_PROVENANCE": value}):
                from importlib import reload

                from muc_one_up import provenance as prov_module

                reload(prov_module)
                assert prov_module.ENABLE_PROVENANCE is False
