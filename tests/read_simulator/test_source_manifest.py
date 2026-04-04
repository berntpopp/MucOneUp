"""Tests for source_manifest stage (Illumina pipeline stage 12a).

Following project testing principles:
- Mock at the module level where names are imported into the stage module
- Test OUR code's logic (argument passing, filtering, error handling)
- NOT testing underlying parsers/trackers themselves
"""

from __future__ import annotations

from unittest.mock import MagicMock

from muc_one_up.read_simulator.stages.source_manifest import generate_read_manifest

MODULE = "muc_one_up.read_simulator.stages.source_manifest"


class TestGenerateReadManifest:
    """Tests for generate_read_manifest."""

    def test_returns_none_when_no_tracker(self, tmp_path):
        """source_tracker=None → returns None immediately."""
        result = generate_read_manifest(
            sidecar_path="some_sidecar.tsv",
            final_bam="some.bam",
            input_fa="some.fa",
            source_tracker=None,
            output_dir=tmp_path,
            output_base="sample",
        )
        assert result is None

    def test_returns_none_when_no_sidecar(self, tmp_path):
        """sidecar_path=None → returns None immediately."""
        tracker = MagicMock()
        result = generate_read_manifest(
            sidecar_path=None,
            final_bam="some.bam",
            input_fa="some.fa",
            source_tracker=tracker,
            output_dir=tmp_path,
            output_base="sample",
        )
        assert result is None

    def test_generates_manifest(self, mocker, tmp_path):
        """Happy-path: mock parse_illumina_reads and pysam; verify manifest written."""
        # Create a real FASTA with 2 sequences
        fasta_path = tmp_path / "input.fa"
        fasta_path.write_text(">haplotype_1\nACGTACGTACGT\n>haplotype_2\nTGCATGCATGCA\n")

        # Create a sidecar file that exists on disk (so unlink doesn't fail)
        sidecar_path = tmp_path / "fragment_origins.tsv"
        sidecar_path.write_text("fragment_index\tchrom\tfstart\tfend\tstrand\n")

        # Mock parse_illumina_reads
        fake_origin = MagicMock()
        fake_origin.read_id = "read1"
        mock_parse = mocker.patch(
            f"{MODULE}.parse_illumina_reads",
            return_value=[fake_origin],
        )

        # Mock pysam: surviving read ids include "read1"
        mock_pysam = mocker.patch(f"{MODULE}.pysam")
        mock_bam = MagicMock()
        mock_bam.__iter__ = MagicMock(return_value=iter([]))
        mock_pysam.AlignmentFile.return_value.__enter__ = MagicMock(return_value=mock_bam)
        mock_pysam.AlignmentFile.return_value.__exit__ = MagicMock(return_value=False)

        # Build a read with query_name for pysam iteration
        mock_read = MagicMock()
        mock_read.query_name = "read1"
        mock_bam.__iter__ = MagicMock(return_value=iter([mock_read]))

        # source_tracker mock
        annotated = [MagicMock()]
        tracker = MagicMock()
        tracker.annotate_reads.return_value = annotated

        result = generate_read_manifest(
            sidecar_path=str(sidecar_path),
            final_bam="final.bam",
            input_fa=str(fasta_path),
            source_tracker=tracker,
            output_dir=tmp_path,
            output_base="sample",
        )

        # Verify result is a non-None string containing "read_manifest"
        assert result is not None
        assert "read_manifest" in result

        # Verify parse_illumina_reads was called with correct seq_name_to_haplotype
        mock_parse.assert_called_once()
        call_kwargs = mock_parse.call_args
        seq_map = call_kwargs.kwargs.get(
            "seq_name_to_haplotype", call_kwargs.args[2] if len(call_kwargs.args) > 2 else None
        )
        assert seq_map is not None
        assert "haplotype_1" in seq_map
        assert "haplotype_2" in seq_map

        # Verify tracker.write_manifest was called
        tracker.write_manifest.assert_called_once()

        # Verify sidecar file was cleaned up
        assert not sidecar_path.exists()
