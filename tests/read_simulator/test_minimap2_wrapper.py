"""Tests for generic minimap2 wrapper - focuses on OUR preset handling logic.

Following Phase 2 testing principles:
- Mock ONLY at system boundary (run_command, convert_sam_to_bam, sort_and_index_bam)
- Test OUR code's logic (preset selection, command construction, file management)
- NOT testing minimap2 itself

These tests verify that our wrapper correctly:
1. Constructs minimap2 commands with technology-specific presets
2. Handles ONT, PacBio CLR, and PacBio HiFi presets
3. Performs SAM→BAM conversion and sorting
4. Validates input/output files
5. Cleans up intermediate files
"""

import pytest

from muc_one_up.exceptions import FileOperationError
from muc_one_up.read_simulator.constants import (
    MINIMAP2_PRESET_ONT,
    MINIMAP2_PRESET_PACBIO_CLR,
    MINIMAP2_PRESET_PACBIO_HIFI,
)
from muc_one_up.read_simulator.wrappers.minimap2_wrapper import (
    align_reads_with_minimap2,
    validate_preset,
)


class TestAlignReadsWithMinimap2:
    """Test align_reads_with_minimap2 preset handling and command construction."""

    def test_constructs_command_with_ont_preset(self, mocker, tmp_path):
        """Test minimap2 command with map-ont preset."""
        # Arrange
        ref_fa = tmp_path / "ref.fa"
        reads_fastq = tmp_path / "reads.fastq"
        output_bam = tmp_path / "aligned.bam"

        ref_fa.write_text(">chr1\nACGT\n")
        reads_fastq.write_text("@read1\nACGT\n+\nIIII\n")

        # Mock dependencies
        mock_run_cmd = mocker.patch(
            "muc_one_up.read_simulator.wrappers.minimap2_wrapper.run_command"
        )
        mock_convert = mocker.patch(
            "muc_one_up.read_simulator.wrappers.minimap2_wrapper.convert_sam_to_bam"
        )
        mock_sort = mocker.patch(
            "muc_one_up.read_simulator.wrappers.minimap2_wrapper.sort_and_index_bam"
        )

        # Create SAM file when run_command is called
        def create_sam(*args, **kwargs):
            sam_path = tmp_path / "aligned.sam"
            sam_path.write_text("SAM data")

        mock_run_cmd.side_effect = create_sam

        # Create BAM when convert is called
        def create_bam(*args, **kwargs):
            output_bam.write_text("BAM data")
            return str(output_bam)

        mock_convert.side_effect = create_bam
        mock_sort.return_value = str(output_bam)

        # Act
        result = align_reads_with_minimap2(
            minimap2_cmd="minimap2",
            samtools_cmd="samtools",
            reference=str(ref_fa),
            reads_fastq=str(reads_fastq),
            output_bam=str(output_bam),
            preset=MINIMAP2_PRESET_ONT,
            threads=4,
        )

        # Assert: Command includes ONT preset
        cmd = mock_run_cmd.call_args[0][0]
        assert "-ax" in cmd
        assert "map-ont" in cmd
        assert str(ref_fa) in cmd
        assert str(reads_fastq) in cmd

        # Verify conversion and sorting were called
        assert mock_convert.called
        assert mock_sort.called
        assert result == str(output_bam)

    def test_constructs_command_with_pacbio_hifi_preset(self, mocker, tmp_path):
        """Test minimap2 command with map-hifi preset."""
        # Arrange
        ref_fa = tmp_path / "ref.fa"
        reads_fastq = tmp_path / "reads.fastq"
        output_bam = tmp_path / "aligned.bam"

        ref_fa.write_text(">chr1\nACGT\n")
        reads_fastq.write_text("@read1\nACGT\n+\nIIII\n")

        # Mock dependencies
        mock_run_cmd = mocker.patch(
            "muc_one_up.read_simulator.wrappers.minimap2_wrapper.run_command"
        )
        mock_convert = mocker.patch(
            "muc_one_up.read_simulator.wrappers.minimap2_wrapper.convert_sam_to_bam"
        )
        mock_sort = mocker.patch(
            "muc_one_up.read_simulator.wrappers.minimap2_wrapper.sort_and_index_bam"
        )

        # Create files
        def create_sam(*args, **kwargs):
            sam_path = tmp_path / "aligned.sam"
            sam_path.write_text("SAM data")

        mock_run_cmd.side_effect = create_sam

        def create_bam(*args, **kwargs):
            output_bam.write_text("BAM data")
            return str(output_bam)

        mock_convert.side_effect = create_bam
        mock_sort.return_value = str(output_bam)

        # Act
        align_reads_with_minimap2(
            minimap2_cmd="minimap2",
            samtools_cmd="samtools",
            reference=str(ref_fa),
            reads_fastq=str(reads_fastq),
            output_bam=str(output_bam),
            preset=MINIMAP2_PRESET_PACBIO_HIFI,  # HiFi preset
        )

        # Assert: Command includes HiFi preset
        cmd = mock_run_cmd.call_args[0][0]
        assert "-ax" in cmd
        assert "map-hifi" in cmd

    def test_validates_reference_exists(self, tmp_path):
        """Test that missing reference file is detected."""
        # Arrange
        reads_fastq = tmp_path / "reads.fastq"
        reads_fastq.write_text("@read1\nACGT\n+\nIIII\n")
        missing_ref = tmp_path / "nonexistent.fa"

        # Act & Assert
        with pytest.raises(FileOperationError, match="Reference genome not found"):
            align_reads_with_minimap2(
                minimap2_cmd="minimap2",
                samtools_cmd="samtools",
                reference=str(missing_ref),  # Doesn't exist!
                reads_fastq=str(reads_fastq),
                output_bam=str(tmp_path / "out.bam"),
            )

    def test_validates_fastq_exists(self, tmp_path):
        """Test that missing FASTQ file is detected."""
        # Arrange
        ref_fa = tmp_path / "ref.fa"
        ref_fa.write_text(">chr1\nACGT\n")
        missing_fastq = tmp_path / "nonexistent.fastq"

        # Act & Assert
        with pytest.raises(FileOperationError, match="Input FASTQ file not found"):
            align_reads_with_minimap2(
                minimap2_cmd="minimap2",
                samtools_cmd="samtools",
                reference=str(ref_fa),
                reads_fastq=str(missing_fastq),  # Doesn't exist!
                output_bam=str(tmp_path / "out.bam"),
            )

    def test_cleans_up_intermediate_sam(self, mocker, tmp_path):
        """Test that intermediate SAM file is cleaned up."""
        # Arrange
        ref_fa = tmp_path / "ref.fa"
        reads_fastq = tmp_path / "reads.fastq"
        output_bam = tmp_path / "aligned.bam"
        output_sam = tmp_path / "aligned.sam"

        ref_fa.write_text(">chr1\nACGT\n")
        reads_fastq.write_text("@read1\nACGT\n+\nIIII\n")

        # Mock dependencies
        mock_run_cmd = mocker.patch(
            "muc_one_up.read_simulator.wrappers.minimap2_wrapper.run_command"
        )
        mock_convert = mocker.patch(
            "muc_one_up.read_simulator.wrappers.minimap2_wrapper.convert_sam_to_bam"
        )
        mock_sort = mocker.patch(
            "muc_one_up.read_simulator.wrappers.minimap2_wrapper.sort_and_index_bam"
        )

        # Create SAM
        def create_sam(*args, **kwargs):
            output_sam.write_text("SAM data")

        mock_run_cmd.side_effect = create_sam

        def create_bam(*args, **kwargs):
            output_bam.write_text("BAM data")
            return str(output_bam)

        mock_convert.side_effect = create_bam
        mock_sort.return_value = str(output_bam)

        # Act
        align_reads_with_minimap2(
            minimap2_cmd="minimap2",
            samtools_cmd="samtools",
            reference=str(ref_fa),
            reads_fastq=str(reads_fastq),
            output_bam=str(output_bam),
        )

        # Assert: SAM was deleted
        assert not output_sam.exists()


class TestValidatePreset:
    """Test validate_preset validation logic."""

    def test_validates_ont_preset(self):
        """Test ONT preset is valid."""
        # Should not raise
        validate_preset(MINIMAP2_PRESET_ONT)

    def test_validates_pacbio_clr_preset(self):
        """Test PacBio CLR preset is valid."""
        # Should not raise
        validate_preset(MINIMAP2_PRESET_PACBIO_CLR)

    def test_validates_pacbio_hifi_preset(self):
        """Test PacBio HiFi preset is valid."""
        # Should not raise
        validate_preset(MINIMAP2_PRESET_PACBIO_HIFI)

    def test_rejects_invalid_preset(self):
        """Test invalid preset raises ValueError."""
        with pytest.raises(ValueError, match="Invalid minimap2 preset"):
            validate_preset("invalid-preset")


class TestPrebuiltIndexReuse:
    """#106: reuse a prebuilt minimap2 index instead of re-indexing the FASTA each run."""

    @staticmethod
    def _align(mocker, tmp_path, index_names):
        ref_fa = tmp_path / "ref.fa"
        ref_fa.write_text(">chr1\nACGT\n")
        for name in index_names:
            (tmp_path / name).write_bytes(b"MMI")
        reads = tmp_path / "r.fastq"
        reads.write_text("@r\nACGT\n+\nIIII\n")
        out = tmp_path / "a.bam"
        out.write_bytes(b"BAM")
        mod = "muc_one_up.read_simulator.wrappers.minimap2_wrapper"
        run = mocker.patch(
            f"{mod}.run_command", side_effect=lambda *a, **k: (tmp_path / "a.sam").write_text("SAM")
        )
        mocker.patch(f"{mod}.convert_sam_to_bam", return_value=str(out))
        mocker.patch(f"{mod}.sort_and_index_bam", return_value=str(out))
        align_reads_with_minimap2(
            "minimap2", "samtools", str(ref_fa), str(reads), str(out), MINIMAP2_PRESET_ONT, 2
        )
        return [str(x) for x in run.call_args[0][0]], ref_fa

    def test_uses_fasta_without_index(self, mocker, tmp_path):
        cmd, ref = self._align(mocker, tmp_path, [])
        assert str(ref) in cmd

    def test_prefers_preset_specific_index(self, mocker, tmp_path):
        cmd, ref = self._align(
            mocker, tmp_path, ["ref.fa.mmi", f"ref.fa.{MINIMAP2_PRESET_ONT}.mmi"]
        )
        assert f"{ref}.{MINIMAP2_PRESET_ONT}.mmi" in cmd and str(ref) not in cmd

    def test_ignores_generic_index(self, mocker, tmp_path):
        """A generic .mmi may have been built with another preset's -k/-w (#117)."""
        cmd, ref = self._align(mocker, tmp_path, ["ref.fa.mmi"])
        assert str(ref) in cmd and f"{ref}.mmi" not in cmd

    def test_ignores_index_older_than_fasta(self, mocker, tmp_path):
        """A stale index may hold an older reference sequence (#117)."""
        import os

        index = tmp_path / f"ref.fa.{MINIMAP2_PRESET_ONT}.mmi"
        index.write_bytes(b"MMI")
        os.utime(index, (1_000_000, 1_000_000))
        cmd, ref = self._align(mocker, tmp_path, [])
        assert str(ref) in cmd and str(index) not in cmd
