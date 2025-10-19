"""Tests for reseq wrapper - focuses on OUR command construction logic.

Following Phase 2 testing principles:
- Mock ONLY at system boundary (subprocess.Popen for run_command)
- Test OUR code's logic (command construction, error handling, file I/O)
- NOT testing reseq itself

These tests verify that our wrapper correctly:
1. Constructs reseq commands with proper arguments
2. Validates output files are created
3. Handles timeout scenarios for seqToIllumina
4. Splits interleaved FASTQ files correctly
"""

import gzip
from unittest.mock import Mock

import pytest

from muc_one_up.exceptions import FileOperationError
from muc_one_up.read_simulator.wrappers.reseq_wrapper import (
    create_reads,
    generate_systematic_errors,
    replace_Ns,
    split_reads,
)


class TestReplaceNs:
    """Test replace_Ns command construction."""

    def test_constructs_correct_command(self, mocker, tmp_path, tools_dict):
        """Test that replaceN command is constructed correctly."""
        # Arrange
        input_fa = tmp_path / "input.fa"
        output_fa = tmp_path / "output.fa"
        input_fa.write_text(">seq1\nACGTNNACGT\n")

        # Mock subprocess.Popen (used by run_command)
        mock_proc = Mock()
        mock_proc.returncode = 0
        mock_proc.stdout = Mock()
        mock_proc.stderr = Mock()
        mock_proc.stdout.readline = Mock(return_value=b"")
        mock_proc.stderr.readline = Mock(return_value=b"")
        mock_proc.wait = Mock(return_value=0)
        mock_proc.pid = 12345

        mock_popen = mocker.patch("subprocess.Popen", return_value=mock_proc)

        # Act
        replace_Ns(str(input_fa), str(output_fa), tools_dict)

        # Assert: Verify command construction
        assert mock_popen.call_count == 1
        cmd = mock_popen.call_args[0][0]

        assert cmd[0] == "reseq"
        assert "replaceN" in cmd
        assert "-r" in cmd
        assert "-R" in cmd
        # Verify input/output paths are in command
        assert str(input_fa) in cmd
        assert str(output_fa) in cmd


class TestGenerateSystematicErrors:
    """Test generate_systematic_errors command construction."""

    def test_constructs_correct_illuminaPE_command(self, mocker, tmp_path, tools_dict):
        """Test that illuminaPE command is constructed with correct arguments."""
        # Arrange
        input_fa = tmp_path / "input.fa"
        model_file = tmp_path / "model.txt"
        output_fq = tmp_path / "output.fq"

        input_fa.write_text(">seq1\nACGT\n")
        model_file.write_text("MODEL_DATA")

        # Mock subprocess.Popen
        def mock_popen_side_effect(cmd, **kwargs):
            # Create output when illuminaPE is called
            if "illuminaPE" in cmd:
                output_fq.write_text("@read1\nACGT\n+\nIIII\n")

            mock_proc = Mock()
            mock_proc.returncode = 0
            mock_proc.stdout = Mock()
            mock_proc.stderr = Mock()
            mock_proc.stdout.readline = Mock(return_value=b"")
            mock_proc.stderr.readline = Mock(return_value=b"")
            mock_proc.wait = Mock(return_value=0)
            mock_proc.pid = 12345
            return mock_proc

        mock_popen = mocker.patch("subprocess.Popen", side_effect=mock_popen_side_effect)

        # Act
        generate_systematic_errors(str(input_fa), str(model_file), str(output_fq), tools_dict)

        # Assert: Verify command construction
        assert mock_popen.call_count == 1
        cmd = mock_popen.call_args[0][0]

        assert cmd[0] == "reseq"
        assert "illuminaPE" in cmd
        assert "-r" in cmd
        assert "-s" in cmd
        assert "--stopAfterEstimation" in cmd
        assert "--writeSysError" in cmd
        assert str(output_fq) in cmd

    def test_raises_error_when_output_missing(self, mocker, tmp_path, tools_dict):
        """Test that missing output raises FileOperationError."""
        # Arrange
        input_fa = tmp_path / "input.fa"
        model_file = tmp_path / "model.txt"
        output_fq = tmp_path / "output.fq"

        input_fa.write_text(">seq1\nACGT\n")
        model_file.write_text("MODEL_DATA")

        # Mock Popen but DON'T create output file
        mock_proc = Mock()
        mock_proc.returncode = 0
        mock_proc.stdout = Mock()
        mock_proc.stderr = Mock()
        mock_proc.stdout.readline = Mock(return_value=b"")
        mock_proc.stderr.readline = Mock(return_value=b"")
        mock_proc.wait = Mock(return_value=0)
        mock_proc.pid = 12345

        mocker.patch("subprocess.Popen", return_value=mock_proc)

        # Act & Assert
        with pytest.raises(FileOperationError, match="missing or empty"):
            generate_systematic_errors(str(input_fa), str(model_file), str(output_fq), tools_dict)


class TestCreateReads:
    """Test create_reads command construction and timeout handling."""

    def test_constructs_correct_seqToIllumina_command(self, mocker, tmp_path, tools_dict):
        """Test that seqToIllumina command includes threads and paths."""
        # Arrange
        fragments = tmp_path / "fragments.fa"
        model_file = tmp_path / "model.txt"
        output_fq = tmp_path / "output.fq"

        fragments.write_text(">frag1\nACGT\n")
        model_file.write_text("MODEL_DATA")

        # Mock subprocess.Popen
        def mock_popen_side_effect(cmd, **kwargs):
            # Create output when seqToIllumina is called
            if "seqToIllumina" in cmd:
                output_fq.write_text("@read1\nACGT\n+\nIIII\n")

            mock_proc = Mock()
            mock_proc.returncode = 0
            mock_proc.stdout = Mock()
            mock_proc.stderr = Mock()
            mock_proc.stdout.readline = Mock(return_value=b"")
            mock_proc.stderr.readline = Mock(return_value=b"")
            mock_proc.wait = Mock(return_value=0)
            mock_proc.pid = 12345
            return mock_proc

        mock_popen = mocker.patch("subprocess.Popen", side_effect=mock_popen_side_effect)

        # Act: Use 8 threads
        create_reads(str(fragments), str(model_file), str(output_fq), 8, tools_dict)

        # Assert: Verify command construction
        assert mock_popen.call_count == 1
        cmd = mock_popen.call_args[0][0]

        assert cmd[0] == "reseq"
        assert "seqToIllumina" in cmd
        assert "-j" in cmd
        j_index = cmd.index("-j")
        assert cmd[j_index + 1] == "8"  # Thread count
        assert "-s" in cmd
        assert "-i" in cmd
        assert "-o" in cmd

    def test_handles_timeout_with_valid_output(self, mocker, tmp_path, tools_dict):
        """Test that timeout with existing output logs warning and continues."""
        # Arrange
        fragments = tmp_path / "fragments.fa"
        model_file = tmp_path / "model.txt"
        output_fq = tmp_path / "output.fq"

        fragments.write_text(">frag1\nACGT\n")
        model_file.write_text("MODEL_DATA")
        # Pre-create output file (simulates timeout after file was created)
        output_fq.write_text("@read1\nACGT\n+\nIIII\n")

        # Mock subprocess.Popen to simulate timeout
        def mock_popen_side_effect(cmd, **kwargs):
            mock_proc = Mock()
            mock_proc.returncode = -15  # SIGTERM (timeout)
            mock_proc.stdout = Mock()
            mock_proc.stderr = Mock()
            mock_proc.stdout.readline = Mock(return_value=b"")
            mock_proc.stderr.readline = Mock(return_value=b"")
            mock_proc.wait = Mock(return_value=-15)
            mock_proc.pid = 12345
            return mock_proc

        mocker.patch("subprocess.Popen", side_effect=mock_popen_side_effect)

        # Act: Should log warning but NOT raise error
        create_reads(str(fragments), str(model_file), str(output_fq), 4, tools_dict, timeout=1)

        # Assert: No exception raised, output exists
        assert output_fq.exists()

    def test_raises_error_on_timeout_without_output(self, mocker, tmp_path, tools_dict):
        """Test that timeout without output raises FileOperationError."""
        # Arrange
        fragments = tmp_path / "fragments.fa"
        model_file = tmp_path / "model.txt"
        output_fq = tmp_path / "output.fq"

        fragments.write_text(">frag1\nACGT\n")
        model_file.write_text("MODEL_DATA")
        # DON'T create output file

        # Mock subprocess.Popen to simulate timeout
        def mock_popen_side_effect(cmd, **kwargs):
            mock_proc = Mock()
            mock_proc.returncode = -15  # SIGTERM (timeout)
            mock_proc.stdout = Mock()
            mock_proc.stderr = Mock()
            mock_proc.stdout.readline = Mock(return_value=b"")
            mock_proc.stderr.readline = Mock(return_value=b"")
            mock_proc.wait = Mock(return_value=-15)
            mock_proc.pid = 12345
            return mock_proc

        mocker.patch("subprocess.Popen", side_effect=mock_popen_side_effect)

        # Act & Assert
        with pytest.raises(FileOperationError, match="missing or empty"):
            create_reads(str(fragments), str(model_file), str(output_fq), 4, tools_dict, timeout=1)


class TestSplitReads:
    """Test split_reads FASTQ file splitting logic."""

    def test_splits_interleaved_fastq_correctly(self, tmp_path):
        """Test that interleaved FASTQ is split into two gzipped files."""
        # Arrange: Create interleaved FASTQ (8 lines: 2 pairs of 4 lines each)
        interleaved = tmp_path / "interleaved.fq"
        out_fq1 = tmp_path / "R1.fq.gz"
        out_fq2 = tmp_path / "R2.fq.gz"

        # Write interleaved FASTQ: read1, read2 (4 lines each)
        interleaved.write_text(
            "@read1/1\nACGTACGT\n+\nIIIIIIII\n"  # Read 1
            "@read1/2\nTGCATGCA\n+\nIIIIIIII\n"  # Read 2
        )

        # Act
        split_reads(str(interleaved), str(out_fq1), str(out_fq2))

        # Assert: Both output files exist
        assert out_fq1.exists()
        assert out_fq2.exists()

        # Verify R1 content
        with gzip.open(out_fq1, "rt") as f1:
            r1_lines = f1.readlines()
            assert len(r1_lines) == 4
            assert r1_lines[0].strip() == "@read1/1"
            assert r1_lines[1].strip() == "ACGTACGT"

        # Verify R2 content
        with gzip.open(out_fq2, "rt") as f2:
            r2_lines = f2.readlines()
            assert len(r2_lines) == 4
            assert r2_lines[0].strip() == "@read1/2"
            assert r2_lines[1].strip() == "TGCATGCA"

    def test_creates_gzipped_output(self, tmp_path):
        """Test that output files are gzip compressed."""
        # Arrange
        interleaved = tmp_path / "interleaved.fq"
        out_fq1 = tmp_path / "R1.fq.gz"
        out_fq2 = tmp_path / "R2.fq.gz"

        interleaved.write_text("@read1/1\nACGT\n+\nIIII\n@read1/2\nTGCA\n+\nIIII\n")

        # Act
        split_reads(str(interleaved), str(out_fq1), str(out_fq2))

        # Assert: Files should be gzip compressed (check magic bytes)
        with out_fq1.open("rb") as f1:
            magic = f1.read(2)
            assert magic == b"\x1f\x8b"  # Gzip magic number

        with out_fq2.open("rb") as f2:
            magic = f2.read(2)
            assert magic == b"\x1f\x8b"  # Gzip magic number
