"""Tests for muc_one_up.cli.analysis module.

Covers: _run_toxic_detection, run_orf_prediction, run_read_simulation,
write_simulation_statistics.
"""

import json
import time
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# ============================================================================
# _run_toxic_detection tests
# ============================================================================


class TestRunToxicDetection:
    """Tests for the _run_toxic_detection helper."""

    def test_extracts_config_parameters(self, minimal_config):
        """Custom toxic_protein_detection values are forwarded to scan_orf_fasta."""
        from muc_one_up.cli.analysis import _run_toxic_detection

        minimal_config["toxic_protein_detection"] = {
            "consensus_motif": "CUSTOM",
            "identity_threshold": 0.9,
            "key_residues": ["A", "B"],
            "expected_repeat_count": 5,
            "weights": {"repeat": 0.7, "composition": 0.3},
            "toxic_cutoff": 0.6,
        }

        with patch(
            "muc_one_up.toxic_protein_detector.scan_orf_fasta", return_value={}
        ) as mock_scan:
            _run_toxic_detection("/fake/orfs.fa", minimal_config)

        mock_scan.assert_called_once()
        kwargs = mock_scan.call_args
        assert kwargs[1]["consensus"] == "CUSTOM"
        assert kwargs[1]["identity_threshold"] == 0.9
        assert kwargs[1]["key_residues"] == ["A", "B"]
        assert kwargs[1]["expected_repeat_count"] == 5
        assert kwargs[1]["w_repeat"] == 0.7
        assert kwargs[1]["w_composition"] == 0.3
        assert kwargs[1]["toxic_detection_cutoff"] == 0.6

    def test_uses_defaults_when_toxic_config_missing(self, minimal_config):
        """Default detection parameters used when toxic_protein_detection absent."""
        from muc_one_up.cli.analysis import _run_toxic_detection

        # Ensure no toxic config
        minimal_config.pop("toxic_protein_detection", None)

        with patch(
            "muc_one_up.toxic_protein_detector.scan_orf_fasta", return_value={}
        ) as mock_scan:
            _run_toxic_detection("/fake/orfs.fa", minimal_config)

        kwargs = mock_scan.call_args[1]
        assert kwargs["consensus"] == "RCHLGPGHQAGPGLHR"
        assert kwargs["identity_threshold"] == 0.8
        assert kwargs["toxic_detection_cutoff"] == 0.5

    def test_uses_correct_assembly_constants(self, minimal_config):
        """Left and right constants from the correct reference assembly are passed."""
        from muc_one_up.cli.analysis import _run_toxic_detection

        with patch(
            "muc_one_up.toxic_protein_detector.scan_orf_fasta", return_value={}
        ) as mock_scan:
            _run_toxic_detection("/fake/orfs.fa", minimal_config)

        call_args = mock_scan.call_args
        expected_left = minimal_config["constants"]["hg38"]["left"]
        expected_right = minimal_config["constants"]["hg38"]["right"]
        assert call_args[1]["left_const"] == expected_left
        assert call_args[1]["right_const"] == expected_right


# ============================================================================
# run_orf_prediction tests
# ============================================================================


class TestRunOrfPrediction:
    """Tests for run_orf_prediction."""

    def test_skips_when_output_orfs_false(self, mock_simulation_options, minimal_config):
        """No ORF work when output_orfs is False."""
        from muc_one_up.cli.analysis import run_orf_prediction

        mock_simulation_options.output_orfs = False

        with patch("muc_one_up.translate.run_orf_finder_in_memory") as mock_orf:
            run_orf_prediction(
                mock_simulation_options,
                minimal_config,
                str(Path(mock_simulation_options.out_dir)),
                "test",
                1,
                [],
                None,
                False,
            )
            mock_orf.assert_not_called()

    def test_single_mode_calls_orf_finder_once(
        self, mock_simulation_options, minimal_config, sample_haplotype_results, tmp_path
    ):
        """Single mode: one ORF call, one toxic detection, one stats file."""
        from muc_one_up.cli.analysis import run_orf_prediction

        mock_simulation_options.output_orfs = True
        out_dir = str(tmp_path / "out")
        Path(out_dir).mkdir()

        with (
            patch("muc_one_up.translate.run_orf_finder_in_memory") as mock_orf,
            patch("muc_one_up.cli.analysis._run_toxic_detection", return_value={"score": 0.5}),
        ):
            run_orf_prediction(
                mock_simulation_options,
                minimal_config,
                out_dir,
                "test",
                1,
                sample_haplotype_results,
                None,
                False,
            )

        mock_orf.assert_called_once()
        # Stats file should have been written
        stats_files = list(Path(out_dir).glob("*.orf_stats.txt"))
        assert len(stats_files) == 1

    def test_dual_mode_calls_orf_finder_twice(
        self, mock_simulation_options, minimal_config, sample_haplotype_results, tmp_path
    ):
        """Dual mode: two ORF calls, two toxic detections, two stats files."""
        from muc_one_up.cli.analysis import run_orf_prediction

        mock_simulation_options.output_orfs = True
        out_dir = str(tmp_path / "out")
        Path(out_dir).mkdir()

        with (
            patch("muc_one_up.translate.run_orf_finder_in_memory") as mock_orf,
            patch("muc_one_up.cli.analysis._run_toxic_detection", return_value={"score": 0.5}),
        ):
            run_orf_prediction(
                mock_simulation_options,
                minimal_config,
                out_dir,
                "test",
                1,
                sample_haplotype_results,
                sample_haplotype_results,
                True,
            )

        assert mock_orf.call_count == 2
        stats_files = list(Path(out_dir).glob("*.orf_stats.txt"))
        assert len(stats_files) == 2

    def test_orf_failure_raises_file_operation_error(
        self, mock_simulation_options, minimal_config, sample_haplotype_results, tmp_path
    ):
        """FileOperationError raised when ORF finder fails."""
        from muc_one_up.cli.analysis import run_orf_prediction
        from muc_one_up.exceptions import FileOperationError

        mock_simulation_options.output_orfs = True
        out_dir = str(tmp_path / "out")
        Path(out_dir).mkdir()

        with (
            patch(
                "muc_one_up.translate.run_orf_finder_in_memory",
                side_effect=RuntimeError("ORF crash"),
            ),
            pytest.raises(FileOperationError, match="ORF finding failed"),
        ):
            run_orf_prediction(
                mock_simulation_options,
                minimal_config,
                out_dir,
                "test",
                1,
                sample_haplotype_results,
                None,
                False,
            )

    def test_toxic_import_error_raises_file_operation_error(
        self, mock_simulation_options, minimal_config, sample_haplotype_results, tmp_path
    ):
        """FileOperationError raised when toxic_protein_detector can't be imported."""
        from muc_one_up.cli.analysis import run_orf_prediction
        from muc_one_up.exceptions import FileOperationError

        mock_simulation_options.output_orfs = True
        out_dir = str(tmp_path / "out")
        Path(out_dir).mkdir()

        with (
            patch("muc_one_up.translate.run_orf_finder_in_memory"),
            patch(
                "muc_one_up.cli.analysis._run_toxic_detection",
                side_effect=ImportError("no module"),
            ),
            pytest.raises(FileOperationError, match="import"),
        ):
            run_orf_prediction(
                mock_simulation_options,
                minimal_config,
                out_dir,
                "test",
                1,
                sample_haplotype_results,
                None,
                False,
            )


# ============================================================================
# run_read_simulation tests
# ============================================================================


class TestRunReadSimulation:
    """Tests for run_read_simulation."""

    def test_skips_when_simulate_reads_none(self, mock_simulation_options, minimal_config):
        """No pipeline call when simulate_reads is None."""
        from muc_one_up.cli.analysis import run_read_simulation

        mock_simulation_options.simulate_reads = None

        with patch("muc_one_up.read_simulation.simulate_reads") as mock_pipeline:
            run_read_simulation(mock_simulation_options, minimal_config, ".", "test", 1, False)
            mock_pipeline.assert_not_called()

    def test_single_mode_calls_pipeline_once(
        self, mock_simulation_options, minimal_config, tmp_path
    ):
        """Single mode: pipeline called once."""
        from muc_one_up.cli.analysis import run_read_simulation

        mock_simulation_options.simulate_reads = "illumina"
        out_dir = str(tmp_path / "out")
        Path(out_dir).mkdir()

        with patch("muc_one_up.read_simulation.simulate_reads") as mock_pipeline:
            run_read_simulation(mock_simulation_options, minimal_config, out_dir, "test", 1, False)

        mock_pipeline.assert_called_once()

    def test_dual_mode_calls_pipeline_twice(
        self, mock_simulation_options, minimal_config, tmp_path
    ):
        """Dual mode: pipeline called for normal and mutated."""
        from muc_one_up.cli.analysis import run_read_simulation

        mock_simulation_options.simulate_reads = "illumina"
        out_dir = str(tmp_path / "out")
        Path(out_dir).mkdir()

        with patch("muc_one_up.read_simulation.simulate_reads") as mock_pipeline:
            run_read_simulation(mock_simulation_options, minimal_config, out_dir, "test", 1, True)

        assert mock_pipeline.call_count == 2

    def test_sets_simulator_in_config(self, mock_simulation_options, minimal_config, tmp_path):
        """Config read_simulation.simulator is set to the requested type."""
        from muc_one_up.cli.analysis import run_read_simulation

        mock_simulation_options.simulate_reads = "ont"
        out_dir = str(tmp_path / "out")
        Path(out_dir).mkdir()

        with patch("muc_one_up.read_simulation.simulate_reads"):
            run_read_simulation(mock_simulation_options, minimal_config, out_dir, "test", 1, False)

        assert minimal_config["read_simulation"]["simulator"] == "ont"

    def test_failure_raises_read_simulation_error(
        self, mock_simulation_options, minimal_config, tmp_path
    ):
        """ReadSimulationError raised when pipeline fails."""
        from muc_one_up.cli.analysis import run_read_simulation
        from muc_one_up.exceptions import ReadSimulationError

        mock_simulation_options.simulate_reads = "illumina"
        out_dir = str(tmp_path / "out")
        Path(out_dir).mkdir()

        with (
            patch(
                "muc_one_up.read_simulation.simulate_reads",
                side_effect=RuntimeError("pipeline crash"),
            ),
            pytest.raises(ReadSimulationError),
        ):
            run_read_simulation(mock_simulation_options, minimal_config, out_dir, "test", 1, False)

    def test_dual_mode_passes_source_trackers(
        self, mock_simulation_options, minimal_config, tmp_path
    ):
        """Dual mode forwards source_tracker and source_tracker_mut."""
        from muc_one_up.cli.analysis import run_read_simulation

        mock_simulation_options.simulate_reads = "illumina"
        out_dir = str(tmp_path / "out")
        Path(out_dir).mkdir()

        tracker_normal = MagicMock()
        tracker_mut = MagicMock()

        with patch("muc_one_up.read_simulation.simulate_reads") as mock_pipeline:
            run_read_simulation(
                mock_simulation_options,
                minimal_config,
                out_dir,
                "test",
                1,
                True,
                source_tracker=tracker_normal,
                source_tracker_mut=tracker_mut,
            )

        # First call (normal) should use tracker_normal
        assert mock_pipeline.call_args_list[0][1]["source_tracker"] == tracker_normal
        # Second call (mutated) should use tracker_mut
        assert mock_pipeline.call_args_list[1][1]["source_tracker"] == tracker_mut

    def test_dual_mode_falls_back_to_source_tracker_when_mut_none(
        self, mock_simulation_options, minimal_config, tmp_path
    ):
        """Dual mode uses source_tracker for mutated when source_tracker_mut is None."""
        from muc_one_up.cli.analysis import run_read_simulation

        mock_simulation_options.simulate_reads = "illumina"
        out_dir = str(tmp_path / "out")
        Path(out_dir).mkdir()

        tracker = MagicMock()

        with patch("muc_one_up.read_simulation.simulate_reads") as mock_pipeline:
            run_read_simulation(
                mock_simulation_options,
                minimal_config,
                out_dir,
                "test",
                1,
                True,
                source_tracker=tracker,
                source_tracker_mut=None,
            )

        # Second call should fall back to tracker
        assert mock_pipeline.call_args_list[1][1]["source_tracker"] == tracker


# ============================================================================
# write_simulation_statistics tests
# ============================================================================


class TestWriteSimulationStatistics:
    """Tests for write_simulation_statistics."""

    def test_single_mode_writes_one_file(
        self, mock_simulation_options, minimal_config, sample_haplotype_results, tmp_path
    ):
        """Single mode produces one stats JSON."""
        from muc_one_up.cli.analysis import write_simulation_statistics

        out_dir = str(tmp_path / "out")
        Path(out_dir).mkdir()

        t_start = time.time()
        t_end = t_start + 1.0

        write_simulation_statistics(
            mock_simulation_options,
            minimal_config,
            out_dir,
            "test",
            1,
            t_start,
            t_end,
            sample_haplotype_results,
            None,
            False,
            None,
            None,
            None,
        )

        stats_files = list(Path(out_dir).glob("*.simulation_stats.json"))
        assert len(stats_files) == 1

        data = json.loads(stats_files[0].read_text())
        assert "haplotypes" in data or "runtime_seconds" in data or "provenance" in data

    def test_dual_mode_writes_two_files(
        self, mock_simulation_options, minimal_config, sample_haplotype_results, tmp_path
    ):
        """Dual mode produces normal + mut stats JSONs."""
        from muc_one_up.cli.analysis import write_simulation_statistics

        out_dir = str(tmp_path / "out")
        Path(out_dir).mkdir()

        t_start = time.time()
        t_end = t_start + 1.0

        write_simulation_statistics(
            mock_simulation_options,
            minimal_config,
            out_dir,
            "test",
            1,
            t_start,
            t_end,
            sample_haplotype_results,
            sample_haplotype_results,
            True,
            ("normal", "dupC"),
            None,
            None,
        )

        stats_files = list(Path(out_dir).glob("*.simulation_stats.json"))
        assert len(stats_files) == 2

        names = {f.name for f in stats_files}
        assert any("normal" in n for n in names)
        assert any("mut" in n for n in names)

    def test_provenance_failure_does_not_crash(
        self, mock_simulation_options, minimal_config, sample_haplotype_results, tmp_path
    ):
        """Graceful degradation when provenance collection fails."""
        from muc_one_up.cli.analysis import write_simulation_statistics

        out_dir = str(tmp_path / "out")
        Path(out_dir).mkdir()

        t_start = time.time()
        t_end = t_start + 1.0

        with patch(
            "muc_one_up.cli.analysis.collect_provenance_metadata",
            side_effect=RuntimeError("provenance crash"),
        ):
            # Should not raise
            write_simulation_statistics(
                mock_simulation_options,
                minimal_config,
                out_dir,
                "test",
                1,
                t_start,
                t_end,
                sample_haplotype_results,
                None,
                False,
                None,
                None,
                None,
            )

        stats_files = list(Path(out_dir).glob("*.simulation_stats.json"))
        assert len(stats_files) == 1

    def test_single_mode_includes_mutation_info(
        self, mock_simulation_options, minimal_config, sample_haplotype_results, tmp_path
    ):
        """Mutation info included in stats when mutation_name set."""
        from muc_one_up.cli.analysis import write_simulation_statistics

        mock_simulation_options.mutation_name = "dupC"
        mock_simulation_options.mutation_targets = ["1,3"]
        out_dir = str(tmp_path / "out")
        Path(out_dir).mkdir()

        t_start = time.time()
        t_end = t_start + 1.0

        write_simulation_statistics(
            mock_simulation_options,
            minimal_config,
            out_dir,
            "test",
            1,
            t_start,
            t_end,
            sample_haplotype_results,
            None,
            False,
            None,
            None,
            None,
        )

        stats_files = list(Path(out_dir).glob("*.simulation_stats.json"))
        data = json.loads(stats_files[0].read_text())
        # Should contain mutation information somewhere in the stats
        stats_text = json.dumps(data)
        assert "dupC" in stats_text or "mutation" in stats_text
