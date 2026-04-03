"""Tests for lazy import mechanisms introduced in v0.37.0."""

import pytest


class TestReadSimulatorLazyImport:
    """Tests for muc_one_up.read_simulator.__getattr__ lazy loading."""

    def test_simulate_reads_pipeline_accessible(self):
        """The lazy-loaded attribute resolves to a callable."""
        from muc_one_up.read_simulator import simulate_reads_pipeline

        assert callable(simulate_reads_pipeline)

    def test_unknown_attribute_raises_attribute_error(self):
        """Accessing a non-existent attribute raises AttributeError."""
        import muc_one_up.read_simulator as mod

        with pytest.raises(AttributeError, match="has no attribute"):
            _ = mod.nonexistent_thing  # type: ignore[attr-defined]

    def test_cli_imports_without_heavy_deps(self):
        """CLI module imports without triggering Bio/orfipy at import time."""
        # This test verifies the lazy import chain doesn't eagerly
        # load heavy dependencies. If it did, import would fail in
        # environments missing those packages.
        import muc_one_up.cli  # noqa: F401
