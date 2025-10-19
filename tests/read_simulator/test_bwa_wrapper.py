"""Tests for BWA wrapper - focuses on OUR command construction logic.

Following Phase 2 testing principles:
- Mock ONLY at system boundary (subprocess.Popen, subprocess.run)
- Test OUR code's logic (command construction, error handling, file validation)
- NOT testing external tools (BWA, samtools) - that's integration testing

These tests verify that our wrapper correctly:
1. Constructs BWA mem commands
2. Chains BWA → samtools view using subprocess.Popen
3. Sorts and indexes BAM files
4. Validates input files
5. Handles tool failures appropriately
"""

from pathlib import Path
from unittest.mock import Mock

import pytest

from muc_one_up.exceptions import ExternalToolError, FileOperationError
from muc_one_up.read_simulator.wrappers.bwa_wrapper import align_reads


class TestBWAWrapperCommandConstruction:
    """Test BWA wrapper command construction (tests OUR code)."""

    def test_constructs_correct_bwa_mem_command(self, mocker, tmp_path, tools_dict):
        """Test that BWA mem command is constructed with correct arguments."""
        # Arrange: Create real test files
        ref = tmp_path / "ref.fa"
        r1 = tmp_path / "R1.fq"
        r2 = tmp_path / "R2.fq"
        out = tmp_path / "out.bam"

        ref.write_text(">chr1\nACGT")
        r1.write_text("@read1\nACGT\n+\nIIII")
        r2.write_text("@read1\nACGT\n+\nIIII")

        # Mock subprocess.Popen for BWA | samtools pipeline
        mock_bwa = Mock()
        mock_bwa.returncode = 0
        mock_bwa.stdout = Mock()
        mock_bwa.communicate = Mock(return_value=(b"", b""))

        mock_samtools = Mock()
        mock_samtools.returncode = 0
        mock_samtools.communicate = Mock(return_value=(b"", b""))

        mock_popen = mocker.patch("subprocess.Popen", side_effect=[mock_bwa, mock_samtools])

        # Mock subprocess.run for sort and index - create output files
        def mock_run_side_effect(cmd, **kwargs):
            if "sort" in cmd:
                out.write_bytes(b"SORTED_BAM_DATA")
            elif "index" in cmd:
                Path(str(out) + ".bai").write_bytes(b"INDEX_DATA")
            return Mock(returncode=0, stderr=b"")

        mocker.patch("subprocess.run", side_effect=mock_run_side_effect)

        # Act: Call our wrapper
        align_reads(
            read1=str(r1),
            read2=str(r2),
            human_reference=str(ref),
            output_bam=str(out),
            tools=tools_dict,
            threads=4,
        )

        # Assert: Verify BWA mem command construction (tests OUR logic)
        assert mock_popen.call_count == 2  # BWA + samtools

        bwa_call = mock_popen.call_args_list[0]
        bwa_cmd = bwa_call[0][0]  # First positional argument

        # Our code should construct this exact command
        assert bwa_cmd == ["bwa", "mem", "-t", "4", str(ref), str(r1), str(r2)]

    def test_constructs_correct_samtools_view_command(self, mocker, tmp_path, tools_dict):
        """Test that samtools view command is constructed correctly."""
        # Arrange
        ref = tmp_path / "ref.fa"
        r1 = tmp_path / "R1.fq"
        r2 = tmp_path / "R2.fq"
        out = tmp_path / "out.bam"

        ref.write_text(">chr1\nACGT")
        r1.write_text("@read1\nACGT\n+\nIIII")
        r2.write_text("@read1\nACGT\n+\nIIII")

        mock_bwa = Mock(returncode=0, stdout=Mock(), communicate=Mock(return_value=(b"", b"")))
        mock_samtools = Mock(returncode=0, communicate=Mock(return_value=(b"", b"")))

        mock_popen = mocker.patch("subprocess.Popen", side_effect=[mock_bwa, mock_samtools])

        # Mock subprocess.run to create output files
        def mock_run_side_effect(cmd, **kwargs):
            if "sort" in cmd:
                out.write_bytes(b"SORTED_BAM_DATA")
            elif "index" in cmd:
                Path(str(out) + ".bai").write_bytes(b"INDEX_DATA")
            return Mock(returncode=0, stderr=b"")

        mocker.patch("subprocess.run", side_effect=mock_run_side_effect)

        # Act
        align_reads(str(r1), str(r2), str(ref), str(out), tools_dict, threads=4)

        # Assert: Verify samtools view command (tests OUR logic)
        samtools_call = mock_popen.call_args_list[1]
        samtools_cmd = samtools_call[0][0]

        assert samtools_cmd == ["samtools", "view", "-bS", "-"]

    def test_uses_correct_thread_count(self, mocker, tmp_path, tools_dict):
        """Test that thread count is correctly passed to BWA."""
        # Arrange
        ref = tmp_path / "ref.fa"
        r1 = tmp_path / "R1.fq"
        r2 = tmp_path / "R2.fq"
        out = tmp_path / "out.bam"

        ref.write_text(">chr1\nACGT")
        r1.write_text("@read1\nACGT\n+\nIIII")
        r2.write_text("@read1\nACGT\n+\nIIII")

        mock_bwa = Mock(returncode=0, stdout=Mock(), communicate=Mock(return_value=(b"", b"")))
        mock_samtools = Mock(returncode=0, communicate=Mock(return_value=(b"", b"")))

        mock_popen = mocker.patch("subprocess.Popen", side_effect=[mock_bwa, mock_samtools])

        # Mock subprocess.run to create output files
        def mock_run_side_effect(cmd, **kwargs):
            if "sort" in cmd:
                out.write_bytes(b"SORTED_BAM_DATA")
            elif "index" in cmd:
                Path(str(out) + ".bai").write_bytes(b"INDEX_DATA")
            return Mock(returncode=0, stderr=b"")

        mock_run = mocker.patch("subprocess.run", side_effect=mock_run_side_effect)

        # Act: Use 8 threads
        align_reads(
            str(r1),
            str(r2),
            str(ref),
            str(out),
            tools_dict,
            threads=8,
        )

        # Assert: BWA command should have -t 8
        bwa_call = mock_popen.call_args_list[0]
        bwa_cmd = bwa_call[0][0]

        assert "-t" in bwa_cmd
        t_index = bwa_cmd.index("-t")
        assert bwa_cmd[t_index + 1] == "8"

        # Assert: samtools sort should also have 8 threads
        sort_call = mock_run.call_args_list[0]
        sort_cmd = sort_call[0][0]

        assert "-@" in sort_cmd
        at_index = sort_cmd.index("-@")
        assert sort_cmd[at_index + 1] == "8"


class TestBWAWrapperInputValidation:
    """Test that BWA wrapper validates inputs (tests OUR code)."""

    def test_raises_error_when_read1_missing(self, tmp_path, tools_dict):
        """Test that missing read1 file raises FileOperationError."""
        ref = tmp_path / "ref.fa"
        r2 = tmp_path / "R2.fq"

        ref.write_text(">chr1\nACGT")
        r2.write_text("@read1\nACGT\n+\nIIII")

        with pytest.raises(FileOperationError, match="Input file not found"):
            align_reads(
                read1="nonexistent.fq",
                read2=str(r2),
                human_reference=str(ref),
                output_bam=str(tmp_path / "out.bam"),
                tools=tools_dict,
                threads=4,
            )

    def test_raises_error_when_read2_missing(self, tmp_path, tools_dict):
        """Test that missing read2 file raises FileOperationError."""
        ref = tmp_path / "ref.fa"
        r1 = tmp_path / "R1.fq"

        ref.write_text(">chr1\nACGT")
        r1.write_text("@read1\nACGT\n+\nIIII")

        with pytest.raises(FileOperationError, match="Input file not found"):
            align_reads(
                read1=str(r1),
                read2="nonexistent.fq",
                human_reference=str(ref),
                output_bam=str(tmp_path / "out.bam"),
                tools=tools_dict,
                threads=4,
            )

    def test_raises_error_when_reference_missing(self, tmp_path, tools_dict):
        """Test that missing reference file raises FileOperationError."""
        r1 = tmp_path / "R1.fq"
        r2 = tmp_path / "R2.fq"

        r1.write_text("@read1\nACGT\n+\nIIII")
        r2.write_text("@read1\nACGT\n+\nIIII")

        with pytest.raises(FileOperationError, match="Input file not found"):
            align_reads(
                read1=str(r1),
                read2=str(r2),
                human_reference="nonexistent.fa",
                output_bam=str(tmp_path / "out.bam"),
                tools=tools_dict,
                threads=4,
            )


class TestBWAWrapperErrorHandling:
    """Test that BWA wrapper handles tool failures correctly (tests OUR code)."""

    def test_raises_error_when_bwa_fails(self, mocker, tmp_path, tools_dict):
        """Test that BWA failure is properly caught and reported."""
        # Arrange
        ref = tmp_path / "ref.fa"
        r1 = tmp_path / "R1.fq"
        r2 = tmp_path / "R2.fq"

        ref.write_text(">chr1\nACGT")
        r1.write_text("@read1\nACGT\n+\nIIII")
        r2.write_text("@read1\nACGT\n+\nIIII")

        # Mock BWA failure
        mock_bwa = Mock()
        mock_bwa.returncode = 1  # Non-zero = failure
        mock_bwa.stdout = Mock()
        mock_bwa.communicate = Mock(return_value=(b"", b"BWA alignment failed"))

        mocker.patch("subprocess.Popen", return_value=mock_bwa)

        # Act & Assert: Should raise ExternalToolError with details
        with pytest.raises(ExternalToolError, match=r"bwa mem.*failed"):
            align_reads(
                str(r1),
                str(r2),
                str(ref),
                str(tmp_path / "out.bam"),
                tools_dict,
                threads=4,
            )

    def test_raises_error_when_samtools_view_fails(self, mocker, tmp_path, tools_dict):
        """Test that samtools view failure is properly caught and reported."""
        # Arrange
        ref = tmp_path / "ref.fa"
        r1 = tmp_path / "R1.fq"
        r2 = tmp_path / "R2.fq"

        ref.write_text(">chr1\nACGT")
        r1.write_text("@read1\nACGT\n+\nIIII")
        r2.write_text("@read1\nACGT\n+\nIIII")

        # BWA succeeds
        mock_bwa = Mock(returncode=0, stdout=Mock(), communicate=Mock(return_value=(b"", b"")))

        # samtools view fails
        mock_samtools = Mock()
        mock_samtools.returncode = 1
        mock_samtools.communicate = Mock(return_value=(b"", b"samtools error"))

        mocker.patch("subprocess.Popen", side_effect=[mock_bwa, mock_samtools])

        # Act & Assert
        with pytest.raises(ExternalToolError, match=r"samtools view.*failed"):
            align_reads(
                str(r1),
                str(r2),
                str(ref),
                str(tmp_path / "out.bam"),
                tools_dict,
                threads=4,
            )

    def test_raises_error_when_sort_fails(self, mocker, tmp_path, tools_dict):
        """Test that samtools sort failure is properly caught and reported."""
        # Arrange
        ref = tmp_path / "ref.fa"
        r1 = tmp_path / "R1.fq"
        r2 = tmp_path / "R2.fq"

        ref.write_text(">chr1\nACGT")
        r1.write_text("@read1\nACGT\n+\nIIII")
        r2.write_text("@read1\nACGT\n+\nIIII")

        # BWA and samtools view succeed
        mock_bwa = Mock(returncode=0, stdout=Mock(), communicate=Mock(return_value=(b"", b"")))
        mock_samtools = Mock(returncode=0, communicate=Mock(return_value=(b"", b"")))

        mocker.patch("subprocess.Popen", side_effect=[mock_bwa, mock_samtools])

        # samtools sort fails
        from subprocess import CalledProcessError

        mocker.patch(
            "subprocess.run",
            side_effect=CalledProcessError(1, "samtools sort", stderr=b"sort failed"),
        )

        # Act & Assert
        with pytest.raises(ExternalToolError, match=r"samtools sort.*failed"):
            align_reads(
                str(r1),
                str(r2),
                str(ref),
                str(tmp_path / "out.bam"),
                tools_dict,
                threads=4,
            )


class TestBWAWrapperOutputValidation:
    """Test that BWA wrapper validates outputs (tests OUR code)."""

    def test_creates_output_bam_file(self, mocker, tmp_path, tools_dict):
        """Test that output BAM file is created."""
        # Arrange
        ref = tmp_path / "ref.fa"
        r1 = tmp_path / "R1.fq"
        r2 = tmp_path / "R2.fq"
        out = tmp_path / "out.bam"

        ref.write_text(">chr1\nACGT")
        r1.write_text("@read1\nACGT\n+\nIIII")
        r2.write_text("@read1\nACGT\n+\nIIII")

        mock_bwa = Mock(returncode=0, stdout=Mock(), communicate=Mock(return_value=(b"", b"")))
        mock_samtools = Mock(returncode=0, communicate=Mock(return_value=(b"", b"")))

        mocker.patch("subprocess.Popen", side_effect=[mock_bwa, mock_samtools])

        # Mock subprocess.run to create output files
        def mock_run_side_effect(cmd, **kwargs):
            # Create the output files that tools would create
            if "sort" in cmd:
                out.write_bytes(b"SORTED_BAM_DATA")
            elif "index" in cmd:
                index_file = Path(str(out) + ".bai")
                index_file.write_bytes(b"INDEX_DATA")
            return Mock(returncode=0, stderr=b"")

        mocker.patch("subprocess.run", side_effect=mock_run_side_effect)

        # Act
        align_reads(str(r1), str(r2), str(ref), str(out), tools_dict, threads=4)

        # Assert: Output files should exist
        assert out.exists()
        assert Path(str(out) + ".bai").exists()

    def test_raises_error_when_output_bam_missing(self, mocker, tmp_path, tools_dict):
        """Test that missing output BAM raises FileOperationError."""
        # Arrange
        ref = tmp_path / "ref.fa"
        r1 = tmp_path / "R1.fq"
        r2 = tmp_path / "R2.fq"
        out = tmp_path / "out.bam"

        ref.write_text(">chr1\nACGT")
        r1.write_text("@read1\nACGT\n+\nIIII")
        r2.write_text("@read1\nACGT\n+\nIIII")

        mock_bwa = Mock(returncode=0, stdout=Mock(), communicate=Mock(return_value=(b"", b"")))
        mock_samtools = Mock(returncode=0, communicate=Mock(return_value=(b"", b"")))

        mocker.patch("subprocess.Popen", side_effect=[mock_bwa, mock_samtools])

        # Mock subprocess.run but DON'T create output files
        mocker.patch("subprocess.run", return_value=Mock(returncode=0, stderr=b""))

        # Act & Assert: Should raise error because BAM file wasn't created
        with pytest.raises(FileOperationError, match="output file missing or empty"):
            align_reads(str(r1), str(r2), str(ref), str(out), tools_dict, threads=4)
