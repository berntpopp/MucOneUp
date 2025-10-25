"""Tests for samtools wrapper - focuses on OUR command construction logic.

Following Phase 2 testing principles:
- Mock ONLY at system boundary (subprocess.run, subprocess.Popen)
- Test OUR code's logic (command construction, error handling)
- NOT testing samtools itself

These tests verify that our wrapper correctly:
1. Constructs samtools commands with proper arguments
2. Handles output redirection securely
3. Calculates coverage statistics correctly
4. Validates inputs and outputs appropriately
"""

from pathlib import Path
from unittest.mock import Mock

import pytest

from muc_one_up.exceptions import FileOperationError
from muc_one_up.read_simulator.wrappers.samtools_wrapper import (
    calculate_target_coverage,
    calculate_vntr_coverage,
    convert_bam_to_fastq,
    convert_sam_to_bam,
    downsample_entire_bam,
    extract_subset_reference,
    sort_and_index_bam,
)


class TestExtractSubsetReference:
    """Test extract_subset_reference command construction."""

    def test_constructs_correct_collate_command(self, mocker, tmp_path, tools_dict):
        """Test that samtools collate command is constructed correctly."""
        # Arrange
        input_bam = tmp_path / "input.bam"
        output_fa = tmp_path / "output.fa"
        input_bam.write_bytes(b"FAKE_BAM_DATA")

        # Mock run_command for collate
        mock_run_cmd = mocker.patch(
            "muc_one_up.read_simulator.wrappers.samtools_wrapper.run_command"
        )

        # Mock subprocess.run for fasta extraction - create output
        def mock_run_side_effect(cmd, **kwargs):
            if "fasta" in cmd:
                output_fa.write_text(">seq1\nACGT\n")
            return Mock(returncode=0, stderr=b"")

        mocker.patch("subprocess.run", side_effect=mock_run_side_effect)

        # Act
        extract_subset_reference(str(input_bam), str(output_fa), tools_dict)

        # Assert: Check collate command (tests OUR logic)
        assert mock_run_cmd.call_count == 1
        collate_call = mock_run_cmd.call_args_list[0]
        collate_cmd = collate_call[0][0]

        assert collate_cmd[0] == "samtools"
        assert "collate" in collate_cmd
        assert "-o" in collate_cmd

    def test_creates_output_fasta_file(self, mocker, tmp_path, tools_dict):
        """Test that output FASTA file is created."""
        # Arrange
        input_bam = tmp_path / "input.bam"
        output_fa = tmp_path / "output.fa"
        input_bam.write_bytes(b"FAKE_BAM_DATA")

        mocker.patch("muc_one_up.read_simulator.wrappers.samtools_wrapper.run_command")

        def mock_run_side_effect(cmd, **kwargs):
            output_fa.write_text(">seq1\nACGT\n")
            return Mock(returncode=0, stderr=b"")

        mocker.patch("subprocess.run", side_effect=mock_run_side_effect)

        # Act
        extract_subset_reference(str(input_bam), str(output_fa), tools_dict)

        # Assert
        assert output_fa.exists()
        assert output_fa.stat().st_size > 0

    def test_raises_error_when_output_missing(self, mocker, tmp_path, tools_dict):
        """Test that missing output FASTA raises FileOperationError."""
        # Arrange
        input_bam = tmp_path / "input.bam"
        output_fa = tmp_path / "output.fa"
        input_bam.write_bytes(b"FAKE_BAM_DATA")

        mocker.patch("muc_one_up.read_simulator.wrappers.samtools_wrapper.run_command")
        # Don't create output file
        mocker.patch("subprocess.run", return_value=Mock(returncode=0, stderr=b""))

        # Act & Assert
        with pytest.raises(FileOperationError, match="missing or empty"):
            extract_subset_reference(str(input_bam), str(output_fa), tools_dict)


class TestCalculateVNTRCoverage:
    """Test calculate_vntr_coverage statistics calculation."""

    def test_constructs_correct_depth_command(self, mocker, tmp_path):
        """Test that samtools depth command is constructed with correct region."""
        # Arrange
        bam = tmp_path / "input.bam"
        bam.write_bytes(b"BAM_DATA")

        depth_file = tmp_path / "test_vntr_depth.txt"
        depth_file.write_text("chr1\t100\t50\nchr1\t101\t50\n")

        # Mock subprocess.run
        def mock_run_side_effect(cmd, **kwargs):
            # Create depth file when depth command is run
            if "depth" in cmd:
                depth_file.write_text("chr1\t100\t50\nchr1\t101\t50\n")
            return Mock(returncode=0, stderr=b"")

        mock_run = mocker.patch("subprocess.run", side_effect=mock_run_side_effect)

        # Act
        calculate_vntr_coverage(
            samtools_exe="samtools",
            bam_file=str(bam),
            region="chr1:100-200",
            threads=4,
            output_dir=str(tmp_path),
            output_name="test",
        )

        # Assert: Verify depth command (tests OUR logic)
        assert mock_run.call_count >= 1
        depth_call = mock_run.call_args_list[0]
        depth_cmd = depth_call[0][0]

        assert depth_cmd[0] == "samtools"
        assert "depth" in depth_cmd
        assert "-r" in depth_cmd
        # Region should be in command
        r_index = depth_cmd.index("-r")
        assert depth_cmd[r_index + 1] == "chr1:100-200"

    def test_calculates_mean_coverage_correctly(self, mocker, tmp_path):
        """Test that mean coverage is calculated correctly from depth file."""
        # Arrange
        bam = tmp_path / "input.bam"
        bam.write_bytes(b"BAM_DATA")

        depth_file = tmp_path / "test_vntr_depth.txt"

        def mock_run_side_effect(cmd, **kwargs):
            # Create depth file with known values: depths of 10, 20, 30 = mean 20
            depth_file.write_text("chr1\t100\t10\nchr1\t101\t20\nchr1\t102\t30\n")
            return Mock(returncode=0, stderr=b"")

        mocker.patch("subprocess.run", side_effect=mock_run_side_effect)

        # Act
        coverage, depth_file = calculate_vntr_coverage(
            "samtools", str(bam), "chr1:100-200", 4, str(tmp_path), "test"
        )

        # Assert: Mean of 10, 20, 30 = 20.0
        assert coverage == 20.0
        # Verify depth file was created
        assert Path(depth_file).exists()
        assert Path(depth_file).stat().st_size > 0


class TestCalculateTargetCoverage:
    """Test calculate_target_coverage with BED file."""

    def test_constructs_depth_command_with_bed_file(self, mocker, tmp_path):
        """Test that depth command includes BED file."""
        # Arrange
        bam = tmp_path / "input.bam"
        bed = tmp_path / "targets.bed"
        bam.write_bytes(b"BAM_DATA")
        bed.write_text("chr1\t100\t200\n")

        depth_file = tmp_path / "test_target_depth.txt"

        def mock_run_side_effect(cmd, **kwargs):
            depth_file.write_text("chr1\t150\t25\n")
            return Mock(returncode=0, stderr=b"")

        mock_run = mocker.patch("subprocess.run", side_effect=mock_run_side_effect)

        # Act
        calculate_target_coverage("samtools", str(bam), str(bed), 4, str(tmp_path), "test")

        # Assert
        depth_call = mock_run.call_args_list[0]
        depth_cmd = depth_call[0][0]

        assert "-b" in depth_cmd
        b_index = depth_cmd.index("-b")
        assert depth_cmd[b_index + 1] == str(bed)

    def test_raises_error_when_bed_file_missing(self, tmp_path):
        """Test that missing BED file raises FileOperationError."""
        bam = tmp_path / "input.bam"
        bam.write_bytes(b"BAM_DATA")

        with pytest.raises(FileOperationError, match="BED file not found"):
            calculate_target_coverage(
                "samtools", str(bam), "nonexistent.bed", 4, str(tmp_path), "test"
            )


class TestDownsampleBAM:
    """Test downsample_entire_bam fraction calculation."""

    def test_constructs_fraction_string_correctly(self, mocker, tmp_path):
        """Test that seed.fraction string format is correct."""
        # Arrange
        input_bam = tmp_path / "input.bam"
        output_bam = tmp_path / "output.bam"
        input_bam.write_bytes(b"BAM_DATA")

        # Track commands called
        commands_called = []

        # Mock subprocess.Popen (used by run_command) to create output
        def mock_popen_side_effect(cmd, **kwargs):
            commands_called.append(cmd)
            mock_proc = Mock()
            mock_proc.returncode = 0
            mock_proc.stdout = Mock()
            mock_proc.stderr = Mock()
            mock_proc.stdout.readline = Mock(return_value=b"")
            mock_proc.stderr.readline = Mock(return_value=b"")
            mock_proc.wait = Mock(return_value=0)
            mock_proc.pid = 12345

            # Create output files when appropriate command is run
            if "view" in cmd:
                output_bam.write_bytes(b"DOWNSAMPLED_BAM")
            elif "index" in cmd:
                Path(str(output_bam) + ".bai").write_bytes(b"INDEX")

            return mock_proc

        mocker.patch("subprocess.Popen", side_effect=mock_popen_side_effect)

        # Act: Downsample to 50% with seed 42
        downsample_entire_bam("samtools", str(input_bam), str(output_bam), 0.5, 42, 4)

        # Assert: Check fraction string format (tests OUR logic)
        view_call = [c for c in commands_called if "view" in c]
        assert len(view_call) >= 1
        view_cmd = view_call[0]

        assert "-s" in view_cmd
        s_index = view_cmd.index("-s")
        fraction_str = view_cmd[s_index + 1]

        # Format should be: seed.fraction where fraction is 4 digits
        # 0.5 = 5000/10000 = .5000, so seed.fraction = 42.5000
        assert fraction_str == "42.5000"

    def test_uses_correct_thread_count(self, mocker, tmp_path):
        """Test that thread count is passed correctly."""
        # Arrange
        input_bam = tmp_path / "input.bam"
        output_bam = tmp_path / "output.bam"
        input_bam.write_bytes(b"BAM_DATA")

        commands_called = []

        def mock_popen_side_effect(cmd, **kwargs):
            commands_called.append(cmd)
            mock_proc = Mock()
            mock_proc.returncode = 0
            mock_proc.stdout = Mock()
            mock_proc.stderr = Mock()
            mock_proc.stdout.readline = Mock(return_value=b"")
            mock_proc.stderr.readline = Mock(return_value=b"")
            mock_proc.wait = Mock(return_value=0)
            mock_proc.pid = 12345

            # Create output files
            if "view" in cmd:
                output_bam.write_bytes(b"DOWNSAMPLED_BAM")
            elif "index" in cmd:
                Path(str(output_bam) + ".bai").write_bytes(b"INDEX")

            return mock_proc

        mocker.patch("subprocess.Popen", side_effect=mock_popen_side_effect)

        # Act: Use 8 threads
        downsample_entire_bam("samtools", str(input_bam), str(output_bam), 0.3, 123, 8)

        # Assert
        view_call = [c for c in commands_called if "view" in c]
        assert len(view_call) >= 1
        view_cmd = view_call[0]

        assert "-@" in view_cmd
        at_index = view_cmd.index("-@")
        assert view_cmd[at_index + 1] == "8"


class TestConvertSamToBam:
    """Test convert_sam_to_bam SAM→BAM conversion."""

    def test_constructs_correct_view_command(self, mocker, tmp_path):
        """Test that samtools view command is constructed correctly."""
        # Arrange
        input_sam = tmp_path / "input.sam"
        output_bam = tmp_path / "output.bam"
        input_sam.write_text("SAM data")

        # Mock run_command to create output
        def mock_run_cmd(cmd, **kwargs):
            output_bam.write_bytes(b"BAM data")

        mocker.patch(
            "muc_one_up.read_simulator.wrappers.samtools_wrapper.run_command",
            side_effect=mock_run_cmd,
        )

        # Act
        result = convert_sam_to_bam(
            samtools_cmd="samtools", input_sam=str(input_sam), output_bam=str(output_bam), threads=8
        )

        # Assert: Verify command construction
        from muc_one_up.read_simulator.wrappers.samtools_wrapper import run_command

        cmd = run_command.call_args[0][0]
        assert "samtools" in cmd
        assert "view" in cmd
        assert "-b" in cmd  # BAM output
        assert "-@" in cmd
        assert "8" in cmd  # threads
        assert str(input_sam) in cmd
        assert str(output_bam) in cmd or "-o" in cmd

        # Verify return value
        assert result == str(output_bam)

    def test_validates_input_sam_exists(self, tmp_path):
        """Test that missing input SAM is detected."""
        # Arrange
        missing_sam = tmp_path / "nonexistent.sam"
        output_bam = tmp_path / "output.bam"

        # Act & Assert
        with pytest.raises(FileOperationError, match="Input SAM file not found"):
            convert_sam_to_bam("samtools", str(missing_sam), str(output_bam))

    def test_validates_output_bam_created(self, mocker, tmp_path):
        """Test that missing/empty output BAM raises error."""
        # Arrange
        input_sam = tmp_path / "input.sam"
        output_bam = tmp_path / "output.bam"
        input_sam.write_text("SAM data")

        # Mock run_command but don't create output
        mocker.patch("muc_one_up.read_simulator.wrappers.samtools_wrapper.run_command")

        # Act & Assert
        with pytest.raises(FileOperationError, match=r"Output BAM.*missing or empty"):
            convert_sam_to_bam("samtools", str(input_sam), str(output_bam))

    def test_validates_output_bam_not_empty(self, mocker, tmp_path):
        """Test that empty output BAM is detected."""
        # Arrange
        input_sam = tmp_path / "input.sam"
        output_bam = tmp_path / "output.bam"
        input_sam.write_text("SAM data")

        # Mock run_command to create empty output
        def mock_run_cmd(cmd, **kwargs):
            output_bam.write_text("")  # Empty!

        mocker.patch(
            "muc_one_up.read_simulator.wrappers.samtools_wrapper.run_command",
            side_effect=mock_run_cmd,
        )

        # Act & Assert
        with pytest.raises(FileOperationError, match=r"Output BAM.*missing or empty"):
            convert_sam_to_bam("samtools", str(input_sam), str(output_bam))


class TestConvertBamToFastq:
    """Test convert_bam_to_fastq BAM→FASTQ conversion."""

    def test_constructs_correct_fastq_command(self, mocker, tmp_path):
        """Test that samtools fastq command is constructed correctly."""
        # Arrange
        input_bam = tmp_path / "input.bam"
        output_fastq = tmp_path / "output.fastq"
        input_bam.write_bytes(b"BAM data")

        # Mock run_command to create output
        def mock_run_cmd(cmd, **kwargs):
            output_fastq.write_text("@read1\nACGT\n+\nIIII\n")

        mocker.patch(
            "muc_one_up.read_simulator.wrappers.samtools_wrapper.run_command",
            side_effect=mock_run_cmd,
        )

        # Act
        result = convert_bam_to_fastq(
            samtools_cmd="samtools",
            input_bam=str(input_bam),
            output_fastq=str(output_fastq),
            threads=8,
        )

        # Assert: Verify command construction
        from muc_one_up.read_simulator.wrappers.samtools_wrapper import run_command

        cmd = run_command.call_args[0][0]
        assert "samtools" in cmd
        assert "fastq" in cmd
        assert "-@" in cmd
        assert "8" in cmd  # threads
        assert "-0" in cmd  # Output unpaired reads
        assert str(input_bam) in cmd
        assert str(output_fastq) in cmd

        # Verify return value
        assert result == str(output_fastq)

    def test_validates_input_bam_exists(self, tmp_path):
        """Test that missing input BAM is detected."""
        # Arrange
        missing_bam = tmp_path / "nonexistent.bam"
        output_fastq = tmp_path / "output.fastq"

        # Act & Assert
        with pytest.raises(FileOperationError, match="Input BAM file not found"):
            convert_bam_to_fastq("samtools", str(missing_bam), str(output_fastq))

    def test_validates_output_fastq_created(self, mocker, tmp_path):
        """Test that missing/empty output FASTQ raises error."""
        # Arrange
        input_bam = tmp_path / "input.bam"
        output_fastq = tmp_path / "output.fastq"
        input_bam.write_bytes(b"BAM data")

        # Mock run_command but don't create output
        mocker.patch("muc_one_up.read_simulator.wrappers.samtools_wrapper.run_command")

        # Act & Assert
        with pytest.raises(FileOperationError, match=r"Output FASTQ.*missing or empty"):
            convert_bam_to_fastq("samtools", str(input_bam), str(output_fastq))

    def test_validates_output_fastq_not_empty(self, mocker, tmp_path):
        """Test that empty output FASTQ is detected."""
        # Arrange
        input_bam = tmp_path / "input.bam"
        output_fastq = tmp_path / "output.fastq"
        input_bam.write_bytes(b"BAM data")

        # Mock run_command to create empty output
        def mock_run_cmd(cmd, **kwargs):
            output_fastq.write_text("")  # Empty!

        mocker.patch(
            "muc_one_up.read_simulator.wrappers.samtools_wrapper.run_command",
            side_effect=mock_run_cmd,
        )

        # Act & Assert
        with pytest.raises(FileOperationError, match=r"Output FASTQ.*missing or empty"):
            convert_bam_to_fastq("samtools", str(input_bam), str(output_fastq))


class TestSortAndIndexBam:
    """Test sort_and_index_bam sorting and indexing."""

    def test_constructs_correct_sort_command(self, mocker, tmp_path):
        """Test that samtools sort command is constructed correctly."""
        # Arrange
        input_bam = tmp_path / "input.bam"
        output_bam = tmp_path / "output.bam"
        input_bam.write_bytes(b"BAM data")

        commands_called = []

        # Mock run_command to track calls and create outputs
        def mock_run_cmd(cmd, **kwargs):
            commands_called.append(cmd)
            if "sort" in cmd:
                # Create temporary sorting file
                temp_bam = tmp_path / "output.bam.sorting.bam"
                temp_bam.write_bytes(b"SORTED BAM")
            elif "index" in cmd:
                # Create index file
                (tmp_path / "output.bam.bai").write_bytes(b"INDEX")

        mocker.patch(
            "muc_one_up.read_simulator.wrappers.samtools_wrapper.run_command",
            side_effect=mock_run_cmd,
        )

        # Act
        result = sort_and_index_bam(
            samtools_exe="samtools", input_bam=str(input_bam), output_bam=str(output_bam), threads=8
        )

        # Assert: Verify sort command
        sort_calls = [c for c in commands_called if "sort" in c]
        assert len(sort_calls) >= 1
        sort_cmd = sort_calls[0]

        assert "samtools" in sort_cmd
        assert "sort" in sort_cmd
        assert "-@" in sort_cmd
        assert "8" in sort_cmd
        assert str(input_bam) in sort_cmd

        # Verify return value
        assert result == str(output_bam)

    def test_sorts_by_name_when_requested(self, mocker, tmp_path):
        """Test that sort by name uses -n flag."""
        # Arrange
        input_bam = tmp_path / "input.bam"
        output_bam = tmp_path / "output.bam"
        input_bam.write_bytes(b"BAM data")

        commands_called = []

        def mock_run_cmd(cmd, **kwargs):
            commands_called.append(cmd)
            if "sort" in cmd:
                temp_bam = tmp_path / "output.bam.sorting.bam"
                temp_bam.write_bytes(b"SORTED BAM")

        mocker.patch(
            "muc_one_up.read_simulator.wrappers.samtools_wrapper.run_command",
            side_effect=mock_run_cmd,
        )

        # Act
        sort_and_index_bam(
            samtools_exe="samtools",
            input_bam=str(input_bam),
            output_bam=str(output_bam),
            by_name=True,
        )

        # Assert: Verify -n flag is present
        sort_calls = [c for c in commands_called if "sort" in c]
        assert len(sort_calls) >= 1
        sort_cmd = sort_calls[0]
        assert "-n" in sort_cmd

        # Verify index command was NOT called (name-sorted BAMs aren't indexed)
        index_calls = [c for c in commands_called if "index" in c]
        assert len(index_calls) == 0

    def test_creates_index_for_coordinate_sorted(self, mocker, tmp_path):
        """Test that coordinate-sorted BAMs are indexed."""
        # Arrange
        input_bam = tmp_path / "input.bam"
        output_bam = tmp_path / "output.bam"
        input_bam.write_bytes(b"BAM data")

        commands_called = []

        def mock_run_cmd(cmd, **kwargs):
            commands_called.append(cmd)
            if "sort" in cmd:
                temp_bam = tmp_path / "output.bam.sorting.bam"
                temp_bam.write_bytes(b"SORTED BAM")
            elif "index" in cmd:
                (tmp_path / "output.bam.bai").write_bytes(b"INDEX")

        mocker.patch(
            "muc_one_up.read_simulator.wrappers.samtools_wrapper.run_command",
            side_effect=mock_run_cmd,
        )

        # Act
        sort_and_index_bam(
            samtools_exe="samtools",
            input_bam=str(input_bam),
            output_bam=str(output_bam),
            by_name=False,  # Coordinate sort
        )

        # Assert: Verify index command was called
        index_calls = [c for c in commands_called if "index" in c]
        assert len(index_calls) >= 1
        index_cmd = index_calls[0]
        assert "samtools" in index_cmd
        assert "index" in index_cmd
        assert str(output_bam) in index_cmd

    def test_replaces_input_when_output_is_none(self, mocker, tmp_path):
        """Test that output_bam=None replaces input file."""
        # Arrange
        input_bam = tmp_path / "input.bam"
        input_bam.write_bytes(b"BAM data")

        def mock_run_cmd(cmd, **kwargs):
            if "sort" in cmd:
                temp_bam = tmp_path / "input.bam.sorting.bam"
                temp_bam.write_bytes(b"SORTED BAM")
            elif "index" in cmd:
                (tmp_path / "input.bam.bai").write_bytes(b"INDEX")

        mocker.patch(
            "muc_one_up.read_simulator.wrappers.samtools_wrapper.run_command",
            side_effect=mock_run_cmd,
        )

        # Act
        result = sort_and_index_bam(
            samtools_exe="samtools", input_bam=str(input_bam), output_bam=None
        )

        # Assert: Output should be same as input (in-place sort)
        assert result == str(input_bam)
