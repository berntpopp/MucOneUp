"""Tests for UCSC tools wrapper - focuses on OUR command construction logic.

Following Phase 2 testing principles:
- Mock ONLY at system boundary (subprocess.Popen for run_command)
- Test OUR code's logic (command construction, error handling, file validation)
- NOT testing faToTwoBit or pblat themselves

These tests verify that our wrapper correctly:
1. Constructs UCSC tool commands with proper arguments
2. Validates input and output files
3. Handles pblat error cases (continues if output exists)
"""

from unittest.mock import Mock

import pytest

from muc_one_up.exceptions import FileOperationError
from muc_one_up.read_simulator.wrappers.ucsc_tools_wrapper import fa_to_twobit, run_pblat


class TestFaToTwoBit:
    """Test fa_to_twobit command construction."""

    def test_constructs_correct_command(self, mocker, tmp_path, tools_dict):
        """Test that faToTwoBit command is constructed correctly."""
        # Arrange
        input_fa = tmp_path / "input.fa"
        output_2bit = tmp_path / "output.2bit"
        input_fa.write_text(">chr1\nACGT\n")

        # Mock subprocess.Popen
        def mock_popen_side_effect(cmd, **kwargs):
            # Create output when faToTwoBit is called
            if "faToTwoBit" in cmd:
                output_2bit.write_bytes(b"TWOBIT_DATA")

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
        fa_to_twobit(str(input_fa), str(output_2bit), tools_dict)

        # Assert: Verify command construction
        assert mock_popen.call_count == 1
        cmd = mock_popen.call_args[0][0]

        assert "faToTwoBit" in cmd
        assert str(input_fa) in cmd
        assert str(output_2bit) in cmd

    def test_validates_output_exists(self, mocker, tmp_path, tools_dict):
        """Test that output file is validated after creation."""
        # Arrange
        input_fa = tmp_path / "input.fa"
        output_2bit = tmp_path / "output.2bit"
        input_fa.write_text(">chr1\nACGT\n")

        def mock_popen_side_effect(cmd, **kwargs):
            output_2bit.write_bytes(b"TWOBIT_DATA")

            mock_proc = Mock()
            mock_proc.returncode = 0
            mock_proc.stdout = Mock()
            mock_proc.stderr = Mock()
            mock_proc.stdout.readline = Mock(return_value=b"")
            mock_proc.stderr.readline = Mock(return_value=b"")
            mock_proc.wait = Mock(return_value=0)
            mock_proc.pid = 12345
            return mock_proc

        mocker.patch("subprocess.Popen", side_effect=mock_popen_side_effect)

        # Act
        fa_to_twobit(str(input_fa), str(output_2bit), tools_dict)

        # Assert: Output should exist
        assert output_2bit.exists()
        assert output_2bit.stat().st_size > 0

    def test_raises_error_when_output_missing(self, mocker, tmp_path, tools_dict):
        """Test that missing output raises FileOperationError."""
        # Arrange
        input_fa = tmp_path / "input.fa"
        output_2bit = tmp_path / "output.2bit"
        input_fa.write_text(">chr1\nACGT\n")

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
            fa_to_twobit(str(input_fa), str(output_2bit), tools_dict)


class TestRunPblat:
    """Test run_pblat command construction and error handling."""

    def test_constructs_pblat_command(self, mocker, tmp_path, tools_dict):
        """Test that pblat command is constructed with correct parameters."""
        # Arrange
        twobit = tmp_path / "input.2bit"
        ref = tmp_path / "ref.fa"
        output_psl = tmp_path / "output.psl"

        twobit.write_bytes(b"TWOBIT_DATA")
        ref.write_text(">chr1\nACGT\n")

        def mock_popen_side_effect(cmd, **kwargs):
            # Create output when pblat is called
            if "pblat" in cmd:
                output_psl.write_text("PSL_DATA\n")

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

        # Act: Run with custom parameters
        run_pblat(
            str(twobit),
            str(ref),
            str(output_psl),
            tools_dict,
            threads=8,
            min_score=90,
            min_identity=90,
        )

        # Assert: Verify command construction
        assert mock_popen.call_count == 1
        cmd = mock_popen.call_args[0][0]

        assert "pblat" in cmd
        assert "-threads=8" in cmd
        assert "-minScore=90" in cmd
        assert "-minIdentity=90" in cmd
        assert str(twobit) in cmd
        assert str(ref) in cmd
        assert str(output_psl) in cmd

    def test_validates_input_files(self, tmp_path, tools_dict):
        """Test that missing input files raise FileOperationError."""
        # Arrange: Create only one input file
        twobit = tmp_path / "input.2bit"
        ref = tmp_path / "ref.fa"
        output_psl = tmp_path / "output.psl"

        twobit.write_bytes(b"TWOBIT_DATA")
        # Don't create ref.fa

        # Act & Assert: Should raise error for missing ref
        with pytest.raises(FileOperationError, match="missing or empty"):
            run_pblat(str(twobit), str(ref), str(output_psl), tools_dict)

    def test_handles_pblat_error_with_output(self, mocker, tmp_path, tools_dict):
        """Test that pblat error is tolerated if output exists."""
        # Arrange
        twobit = tmp_path / "input.2bit"
        ref = tmp_path / "ref.fa"
        output_psl = tmp_path / "output.psl"

        twobit.write_bytes(b"TWOBIT_DATA")
        ref.write_text(">chr1\nACGT\n")
        # Pre-create output to simulate error after file was created
        output_psl.write_text("PSL_DATA\n")

        # Mock Popen to simulate error
        def mock_popen_side_effect(cmd, **kwargs):
            mock_proc = Mock()
            mock_proc.returncode = 1  # Non-zero = error
            mock_proc.stdout = Mock()
            mock_proc.stderr = Mock()
            mock_proc.stdout.readline = Mock(return_value=b"")
            mock_proc.stderr.readline = Mock(return_value=b"")
            mock_proc.wait = Mock(return_value=1)
            mock_proc.pid = 12345
            return mock_proc

        mocker.patch("subprocess.Popen", side_effect=mock_popen_side_effect)

        # Act: Should log warning but NOT raise error
        run_pblat(str(twobit), str(ref), str(output_psl), tools_dict)

        # Assert: Output exists, no exception raised
        assert output_psl.exists()

    def test_raises_error_when_no_output(self, mocker, tmp_path, tools_dict):
        """Test that pblat error without output raises FileOperationError."""
        # Arrange
        twobit = tmp_path / "input.2bit"
        ref = tmp_path / "ref.fa"
        output_psl = tmp_path / "output.psl"

        twobit.write_bytes(b"TWOBIT_DATA")
        ref.write_text(">chr1\nACGT\n")
        # DON'T create output file

        # Mock Popen to simulate error
        def mock_popen_side_effect(cmd, **kwargs):
            mock_proc = Mock()
            mock_proc.returncode = 1
            mock_proc.stdout = Mock()
            mock_proc.stderr = Mock()
            mock_proc.stdout.readline = Mock(return_value=b"")
            mock_proc.stderr.readline = Mock(return_value=b"")
            mock_proc.wait = Mock(return_value=1)
            mock_proc.pid = 12345
            return mock_proc

        mocker.patch("subprocess.Popen", side_effect=mock_popen_side_effect)

        # Act & Assert: Should raise error
        with pytest.raises(FileOperationError, match="pblat alignment failed"):
            run_pblat(str(twobit), str(ref), str(output_psl), tools_dict)
