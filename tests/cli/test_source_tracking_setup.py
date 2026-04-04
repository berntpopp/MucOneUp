"""Tests for source tracker construction helper."""

import pytest

from muc_one_up.type_defs import HaplotypeResult, RepeatUnit


@pytest.fixture
def minimal_config():
    """Minimal config for source tracker construction."""
    return {
        "repeats": {"X": "XXXXX"},
        "constants": {"hg38": {"left": "TTTT", "right": "GGGG"}},
        "reference_assembly": "hg38",
    }


@pytest.fixture
def mock_results():
    """Minimal simulation results."""
    return [
        HaplotypeResult(sequence="TTTTXXXXXGGGG", chain=[RepeatUnit("X")]),
        HaplotypeResult(sequence="TTTTXXXXXGGGG", chain=[RepeatUnit("X")]),
    ]


class TestBuildSourceTrackers:
    def test_non_dual_mode_returns_tracker_and_none(self, tmp_path, minimal_config, mock_results):
        """Non-dual mode: returns (tracker, None), writes one coord map."""
        from muc_one_up.cli.source_tracking_setup import build_source_trackers

        tracker, mut_tracker = build_source_trackers(
            config=minimal_config,
            results=mock_results,
            mutated_results=None,
            mutation_positions=None,
            mutation_name=None,
            applied_snp_info_normal=None,
            applied_snp_info_mut=None,
            reference_assembly=None,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            dual_mutation_mode=False,
        )

        assert tracker is not None
        assert mut_tracker is None
        coord_file = tmp_path / "test.001.repeat_coordinates.tsv"
        assert coord_file.exists()

    def test_dual_mode_returns_two_trackers(self, tmp_path, minimal_config, mock_results):
        """Dual mode: returns (tracker, mut_tracker), writes two coord maps."""
        from muc_one_up.cli.source_tracking_setup import build_source_trackers

        tracker, mut_tracker = build_source_trackers(
            config=minimal_config,
            results=mock_results,
            mutated_results=mock_results,
            mutation_positions=None,
            mutation_name="testMut",
            applied_snp_info_normal=None,
            applied_snp_info_mut=None,
            reference_assembly=None,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            dual_mutation_mode=True,
        )

        assert tracker is not None
        assert mut_tracker is not None
        assert (tmp_path / "test.001.normal.repeat_coordinates.tsv").exists()
        assert (tmp_path / "test.001.mut.repeat_coordinates.tsv").exists()

    def test_reference_assembly_defaults_to_hg38(self, tmp_path, minimal_config, mock_results):
        """When reference_assembly is None, falls back to config then hg38."""
        from muc_one_up.cli.source_tracking_setup import build_source_trackers

        tracker, _ = build_source_trackers(
            config=minimal_config,
            results=mock_results,
            mutated_results=None,
            mutation_positions=None,
            mutation_name=None,
            applied_snp_info_normal=None,
            applied_snp_info_mut=None,
            reference_assembly=None,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            dual_mutation_mode=False,
        )

        assert tracker is not None
