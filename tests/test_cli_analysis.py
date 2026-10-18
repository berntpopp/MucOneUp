"""Tests for muc_one_up.cli.analysis module.

Covers: write_simulation_statistics, run_orf_analysis_standalone,
and _translate_dna_flanks.
"""

import json
import time
from pathlib import Path
from unittest.mock import patch

from muc_one_up.cli.analysis import (
    _is_dna,
    _translate_dna_flanks,
    run_orf_analysis_standalone,
    write_simulation_statistics,
)


class TestTranslateDnaFlanks:
    """Tests for DNA flank detection and peptide translation helper."""

    def test_is_dna_helper(self):
        assert _is_dna("ATGAAACCC") is True
        assert _is_dna("ACGTN") is True
        assert _is_dna("MKPGF") is False
        assert _is_dna("A") is False
        assert _is_dna(None) is False

    def test_non_dna_strings_passed_as_is(self):
        """Peptide strings are not translated and returned directly."""
        left, right = _translate_dna_flanks("MKP", "GF")
        assert left == "MKP"
        assert right == "GF"

    def test_none_constants_return_none(self):
        """None constants return None."""
        left, right = _translate_dna_flanks(None, None)
        assert left is None
        assert right is None

    def test_dna_constants_translated(self):
        """DNA sequences are translated into candidate peptides."""
        left, right = _translate_dna_flanks("ATGAAACCC", "GGGTTTTAA")
        assert left == "MKP"
        assert right == "GF"

    def test_dna_constants_matched_against_orf_output(self, tmp_path: Path):
        """Active frame is selected by matching predicted ORF sequences."""
        orf_file = tmp_path / "test.orfs.fa"
        right_dna = "AGTTTAA"
        orf_file.write_text(">seq1\nMKPV\n")

        left, right = _translate_dna_flanks("ATGAAACCC", right_dna, orf_file)
        assert left == "MKP"
        assert right == "V"

    def test_orfipy_header_coordinates_strip_left_flank(self, tmp_path: Path):
        """When orfipy header specifies start < len(left_const), prefix AA count is computed."""
        orf_file = tmp_path / "test.orfs.fa"
        left_dna = "ATGAAACCCTTT"
        orf_file.write_text(">seq1 [3-30](+) type:complete\nKPFAAA\n")

        left, _ = _translate_dna_flanks(left_dna, None, orf_file)
        assert left == "KPF"

    def test_multi_orf_dna_constants_matched_both_frames(self, tmp_path: Path):
        """When multiple ORFs in the FASTA have different matched frames/coordinates,
        _translate_dna_flanks returns lists of matched candidate flanks.
        """
        orf_file = tmp_path / "multi.orfs.fa"
        left_dna = "ATGAAACCCTTT"
        right_dna = "AGTTTAA"
        orf_file.write_text(
            ">seq1 [3-30](+) type:complete\nKPFAAAAV\n>seq2 [0-30](+) type:complete\nMKPFAAAASL\n"
        )
        left, right = _translate_dna_flanks(left_dna, right_dna, orf_file)
        assert isinstance(left, list)
        assert set(left) == {"KPF", "MKPF"}
        assert isinstance(right, list)
        assert set(right) == {"V", "SL"}

    def test_multi_orf_order_independent_toxicity_scores(self, tmp_path: Path):
        """Verifies that two ORFs in different reading frames / coordinates produce
        identical scan_orf_fasta toxicity scores regardless of record order.
        """
        from muc_one_up.toxic_protein_detector import scan_orf_fasta

        left_dna = "ATGAAACCCTTT"
        right_dna = "AGTTTAA"
        motif = "RCHLGPGHQAGPGLHR"

        seq1 = f"KPF{motif * 3}V"
        seq2 = f"MKPF{motif * 2}SL"

        f1 = tmp_path / "order1.fa"
        f1.write_text(
            f">orf1 [3-99](+) type:complete\n{seq1}\n>orf2 [0-99](+) type:complete\n{seq2}\n"
        )
        f2 = tmp_path / "order2.fa"
        f2.write_text(
            f">orf2 [0-99](+) type:complete\n{seq2}\n>orf1 [3-99](+) type:complete\n{seq1}\n"
        )

        left1, right1 = _translate_dna_flanks(left_dna, right_dna, f1)
        left2, right2 = _translate_dna_flanks(left_dna, right_dna, f2)

        res1 = scan_orf_fasta(str(f1), left_const=left1, right_const=right1, consensus=motif)
        res2 = scan_orf_fasta(str(f2), left_const=left2, right_const=right2, consensus=motif)

        assert res1["orf1 [3-99](+) type:complete"] == res2["orf1 [3-99](+) type:complete"]
        assert res1["orf2 [0-99](+) type:complete"] == res2["orf2 [0-99](+) type:complete"]
        assert res1["orf1 [3-99](+) type:complete"]["repeat_count"] == 3
        assert res2["orf2 [0-99](+) type:complete"]["repeat_count"] == 2


class TestRunOrfAnalysisStandalone:
    """Tests for run_orf_analysis_standalone."""

    def test_run_orf_analysis_translates_dna_constants(self, tmp_path: Path):
        """Standalone ORF analysis translates DNA constants before calling scan_orf_fasta."""
        fasta_file = tmp_path / "input.fa"
        fasta_file.write_text(">seq1\nATGAAACCCGGGTTTTAA\n")

        orf_output = tmp_path / "test.orfs.fa"
        orf_output.write_text(">seq1\nMKPGF\n")

        with (
            patch("muc_one_up.read_simulator.utils.common_utils.run_command"),
            patch(
                "muc_one_up.toxic_protein_detector.scan_orf_fasta", return_value={"score": 0.1}
            ) as mock_scan,
        ):
            run_orf_analysis_standalone(
                input_fasta=str(fasta_file),
                out_dir=str(tmp_path),
                out_base="test",
                orf_min_aa=1,
                orf_aa_prefix=None,
                left_const="ATGAAACCC",
                right_const="GGGTTTTAA",
            )

        mock_scan.assert_called_once()
        call_kwargs = mock_scan.call_args[1]
        assert call_kwargs["left_const"] == "MKP"
        assert call_kwargs["right_const"] == "GF"

        stats_file = tmp_path / "test.orf_stats.json"
        assert stats_file.exists()


class TestWriteSimulationStatistics:
    """Tests for write_simulation_statistics."""

    def test_single_mode_writes_one_file(
        self, mock_simulation_options, minimal_config, sample_haplotype_results, tmp_path
    ):
        """Single mode produces one stats JSON."""
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
        out_dir = str(tmp_path / "out")
        Path(out_dir).mkdir()

        t_start = time.time()
        t_end = t_start + 1.0

        with patch(
            "muc_one_up.cli.analysis.collect_provenance_metadata",
            side_effect=RuntimeError("provenance crash"),
        ):
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
        stats_text = json.dumps(data)
        assert "dupC" in stats_text or "mutation" in stats_text

    def test_passes_mutated_units_to_simulation_statistics_and_companion(
        self, mock_simulation_options, minimal_config, sample_haplotype_results, tmp_path
    ):
        """mutated_units is passed to generate_simulation_statistics and companion fields."""
        out_dir = str(tmp_path / "out")
        Path(out_dir).mkdir()

        mutated_units = {1: [(2, "ATCG_MUT")]}
        mock_simulation_options.track_read_source = True
        mock_simulation_options.mutation_name = "dupC"

        with (
            patch("muc_one_up.cli.analysis.generate_simulation_statistics") as mock_gen,
            patch("muc_one_up.cli.analysis._add_companion_fields") as mock_comp,
        ):
            mock_gen.return_value = {}
            write_simulation_statistics(
                mock_simulation_options,
                minimal_config,
                out_dir,
                "test",
                1,
                0.0,
                1.0,
                sample_haplotype_results,
                None,
                False,
                None,
                None,
                None,
                mutated_units=mutated_units,
            )

        assert mock_gen.call_args[1]["mutated_units"] == mutated_units
        assert mock_comp.call_args[1]["mutated_units"] == mutated_units
