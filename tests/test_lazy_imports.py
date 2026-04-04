"""Tests for lazy import mechanisms introduced in v0.37.0."""

import contextlib
import importlib
import sys

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


def test_read_simulation_import_does_not_load_amplicon():
    """Importing read_simulation must not transitively load amplicon_pipeline."""
    mod_name = "muc_one_up.read_simulator.amplicon_pipeline"
    sys.modules.pop(mod_name, None)

    import muc_one_up.read_simulation

    importlib.reload(muc_one_up.read_simulation)

    assert mod_name not in sys.modules, (
        f"{mod_name} was loaded at import time — _get_simulator_map() is still eager"
    )


def test_amplicon_selection_loads_amplicon_module():
    """Selecting simulator='amplicon' must load the amplicon pipeline module."""
    mod_name = "muc_one_up.read_simulator.amplicon_pipeline"
    sys.modules.pop(mod_name, None)

    from muc_one_up.read_simulation import simulate_reads

    config = {"read_simulation": {"simulator": "amplicon"}}
    with contextlib.suppress(Exception):
        simulate_reads(config, "/fake/input.fa")

    assert mod_name in sys.modules, f"{mod_name} was NOT loaded after selecting amplicon simulator"


class TestGetSimulatorDispatch:
    """Verify _get_simulator returns a callable for each valid backend."""

    @pytest.mark.parametrize("sim_type", ["illumina", "ont", "pacbio", "amplicon"])
    def test_valid_simulator_returns_callable(self, sim_type):
        from muc_one_up.read_simulation import _get_simulator

        func = _get_simulator(sim_type)
        assert callable(func)

    def test_unknown_simulator_raises_value_error(self):
        from muc_one_up.read_simulation import _get_simulator

        with pytest.raises(ValueError, match="Unknown simulator"):
            _get_simulator("nonexistent")
