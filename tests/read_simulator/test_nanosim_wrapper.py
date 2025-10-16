"""Tests for NanoSim wrapper - focuses on OUR command construction logic.

Following Phase 2 testing principles:
- Mock ONLY at system boundary (subprocess.Popen, subprocess.run)
- Test OUR code's logic (command construction, error handling, file management)
- NOT testing NanoSim or minimap2 themselves

These tests verify that our wrapper correctly:
1. Constructs NanoSim commands with proper arguments
2. Handles conda/mamba commands using shlex.split
3. Chains minimap2 → samtools commands correctly
4. Manages temporary files and cleanup
5. Validates output files are created
"""

from pathlib import Path
from unittest.mock import Mock

import pytest

from muc_one_up.read_simulator.wrappers.nanosim_wrapper import (
    align_ont_reads_with_minimap2,
    run_nanosim_simulation,
)


class TestRunNanoSimSimulation:
    """Test run_nanosim_simulation command construction."""

    def test_constructs_basic_command(self, mocker, tmp_path):
        """Test that NanoSim genome command is constructed correctly."""
        # Arrange
        ref_fa = tmp_path / "ref.fa"
        model = tmp_path / "model"
        output_prefix = tmp_path / "output"
        output_fastq = tmp_path / "output_aligned_reads.fastq"

        ref_fa.write_text(">chr1\nACGT\n")
        model.mkdir()

        # Mock subprocess.Popen
        def mock_popen_side_effect(cmd, **kwargs):
            # Create output when NanoSim is called
            if "genome" in cmd:
                output_fastq.write_text("@read1\nACGT\n+\nIIII\n")

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
        result = run_nanosim_simulation(
            nanosim_cmd="nanosim",
            reference_fasta=str(ref_fa),
            output_prefix=str(output_prefix),
            training_model=str(model),
            coverage=30.0,
            threads=4,
        )

        # Assert: Verify command construction
        assert mock_popen.call_count == 1
        cmd = mock_popen.call_args[0][0]

        assert "nanosim" in cmd
        assert "genome" in cmd
        assert "-rg" in cmd
        assert str(ref_fa) in cmd
        assert "-c" in cmd
        assert str(model) in cmd
        assert "-o" in cmd
        assert str(output_prefix) in cmd
        assert "-t" in cmd
        assert "4" in cmd
        assert "-x" in cmd
        assert "30.0" in cmd
        assert "--fastq" in cmd

        # Verify output path returned
        assert result == str(output_fastq)

    def test_handles_conda_command(self, mocker, tmp_path):
        """Test that conda/mamba commands are parsed with shlex.split."""
        # Arrange
        ref_fa = tmp_path / "ref.fa"
        model = tmp_path / "model"
        output_prefix = tmp_path / "output"
        output_fastq = tmp_path / "output_aligned_reads.fastq"

        ref_fa.write_text(">chr1\nACGT\n")
        model.mkdir()

        def mock_popen_side_effect(cmd, **kwargs):
            if "genome" in cmd:
                output_fastq.write_text("@read1\nACGT\n+\nIIII\n")

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

        # Act: Use conda command with spaces
        run_nanosim_simulation(
            nanosim_cmd="conda run -n nanosim simulator.py",
            reference_fasta=str(ref_fa),
            output_prefix=str(output_prefix),
            training_model=str(model),
            coverage=30.0,
        )

        # Assert: Command should be parsed correctly
        cmd = mock_popen.call_args[0][0]
        assert "conda" in cmd
        assert "run" in cmd
        assert "-n" in cmd
        assert "nanosim" in cmd
        assert "simulator.py" in cmd

    def test_includes_optional_parameters(self, mocker, tmp_path):
        """Test that min/max read lengths are included when specified."""
        # Arrange
        ref_fa = tmp_path / "ref.fa"
        model = tmp_path / "model"
        output_prefix = tmp_path / "output"
        output_fastq = tmp_path / "output_aligned_reads.fastq"

        ref_fa.write_text(">chr1\nACGT\n")
        model.mkdir()

        def mock_popen_side_effect(cmd, **kwargs):
            if "genome" in cmd:
                output_fastq.write_text("@read1\nACGT\n+\nIIII\n")

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

        # Act: Include min/max read lengths
        run_nanosim_simulation(
            nanosim_cmd="nanosim",
            reference_fasta=str(ref_fa),
            output_prefix=str(output_prefix),
            training_model=str(model),
            coverage=30.0,
            min_read_length=100,
            max_read_length=50000,
        )

        # Assert
        cmd = mock_popen.call_args[0][0]
        assert "--min_len" in cmd
        assert "100" in cmd
        assert "--max_len" in cmd
        assert "50000" in cmd

    def test_finds_aligned_reads_output(self, mocker, tmp_path):
        """Test that primary output pattern is found and returned."""
        # Arrange
        ref_fa = tmp_path / "ref.fa"
        model = tmp_path / "model"
        output_prefix = tmp_path / "output"
        output_fastq = tmp_path / "output_aligned_reads.fastq"

        ref_fa.write_text(">chr1\nACGT\n")
        model.mkdir()

        def mock_popen_side_effect(cmd, **kwargs):
            # Create primary output pattern
            output_fastq.write_text("@read1\nACGT\n+\nIIII\n")

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
        result = run_nanosim_simulation(
            "nanosim", str(ref_fa), str(output_prefix), str(model), 30.0
        )

        # Assert: Returns primary pattern
        assert result == str(output_fastq)

    def test_finds_alternative_output(self, mocker, tmp_path):
        """Test that alternative output pattern is found when primary is missing."""
        # Arrange
        ref_fa = tmp_path / "ref.fa"
        model = tmp_path / "model"
        output_prefix = tmp_path / "output"
        # Create alternative output pattern (without _aligned_reads)
        alt_output_fastq = tmp_path / "output.fastq"

        ref_fa.write_text(">chr1\nACGT\n")
        model.mkdir()

        def mock_popen_side_effect(cmd, **kwargs):
            # Create alternative output pattern
            alt_output_fastq.write_text("@read1\nACGT\n+\nIIII\n")

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
        result = run_nanosim_simulation(
            "nanosim", str(ref_fa), str(output_prefix), str(model), 30.0
        )

        # Assert: Returns alternative pattern
        assert result == str(alt_output_fastq)

    def test_raises_error_when_no_output(self, mocker, tmp_path):
        """Test that missing output raises RuntimeError."""
        # Arrange
        ref_fa = tmp_path / "ref.fa"
        model = tmp_path / "model"
        output_prefix = tmp_path / "output"

        ref_fa.write_text(">chr1\nACGT\n")
        model.mkdir()

        # Mock Popen but DON'T create any output file
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
        with pytest.raises(RuntimeError, match="Expected output file not found"):
            run_nanosim_simulation("nanosim", str(ref_fa), str(output_prefix), str(model), 30.0)


class TestAlignONTReadsWithMinimap2:
    """Test align_ont_reads_with_minimap2 command chaining."""

    def test_constructs_minimap2_command(self, mocker, tmp_path):
        """Test that minimap2 command is constructed correctly."""
        # Arrange
        ref = tmp_path / "ref.fa"
        reads = tmp_path / "reads.fq"
        output_bam = tmp_path / "output.bam"

        ref.write_text(">chr1\nACGT\n")
        reads.write_text("@read1\nACGT\n+\nIIII\n")

        # Track subprocess.run calls (for minimap2)
        run_calls = []

        def mock_run_side_effect(cmd, **kwargs):
            run_calls.append(cmd)
            # Create SAM output if stdout is redirected
            if hasattr(kwargs.get("stdout"), "write"):
                kwargs["stdout"].write("@HD\tVN:1.0\n")
            return Mock(returncode=0, stderr="")

        mocker.patch("subprocess.run", side_effect=mock_run_side_effect)

        # Mock subprocess.Popen for samtools commands
        def mock_popen_side_effect(cmd, **kwargs):
            # Create outputs for samtools commands
            if "view" in cmd:
                Path(str(output_bam) + ".unsorted").write_bytes(b"BAM_DATA")
            elif "sort" in cmd:
                output_bam.write_bytes(b"SORTED_BAM")
            elif "index" in cmd:
                Path(str(output_bam) + ".bai").write_bytes(b"INDEX")

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
        align_ont_reads_with_minimap2(
            minimap2_cmd="minimap2",
            samtools_cmd="samtools",
            human_reference=str(ref),
            reads_fastq=str(reads),
            output_bam=str(output_bam),
            threads=8,
        )

        # Assert: Verify minimap2 command
        assert len(run_calls) >= 1
        minimap2_cmd = run_calls[0]

        assert "minimap2" in minimap2_cmd
        assert "-t" in minimap2_cmd
        assert "8" in minimap2_cmd
        assert "-ax" in minimap2_cmd
        assert "map-ont" in minimap2_cmd
        assert str(ref) in minimap2_cmd
        assert str(reads) in minimap2_cmd

    def test_chains_samtools_commands(self, mocker, tmp_path):
        """Test that samtools view, sort, and index are called in sequence."""
        # Arrange
        ref = tmp_path / "ref.fa"
        reads = tmp_path / "reads.fq"
        output_bam = tmp_path / "output.bam"

        ref.write_text(">chr1\nACGT\n")
        reads.write_text("@read1\nACGT\n+\nIIII\n")

        # Mock subprocess.run (minimap2)
        def mock_run_side_effect(cmd, **kwargs):
            if hasattr(kwargs.get("stdout"), "write"):
                kwargs["stdout"].write("@HD\tVN:1.0\n")
            return Mock(returncode=0, stderr="")

        mocker.patch("subprocess.run", side_effect=mock_run_side_effect)

        # Track Popen calls (samtools)
        popen_calls = []

        def mock_popen_side_effect(cmd, **kwargs):
            popen_calls.append(cmd)

            if "view" in cmd:
                Path(str(output_bam) + ".unsorted").write_bytes(b"BAM_DATA")
            elif "sort" in cmd:
                output_bam.write_bytes(b"SORTED_BAM")
            elif "index" in cmd:
                Path(str(output_bam) + ".bai").write_bytes(b"INDEX")

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
        align_ont_reads_with_minimap2(
            "minimap2", "samtools", str(ref), str(reads), str(output_bam), threads=4
        )

        # Assert: Verify samtools commands
        assert len(popen_calls) == 3

        # Check view command
        view_cmd = popen_calls[0]
        assert "samtools" in view_cmd
        assert "view" in view_cmd
        assert "-F" in view_cmd
        assert "4" in view_cmd  # Filter unmapped

        # Check sort command
        sort_cmd = popen_calls[1]
        assert "samtools" in sort_cmd
        assert "sort" in sort_cmd

        # Check index command
        index_cmd = popen_calls[2]
        assert "samtools" in index_cmd
        assert "index" in index_cmd

    def test_cleans_up_temporary_files(self, mocker, tmp_path):
        """Test that temporary SAM and unsorted BAM are cleaned up."""
        # Arrange
        ref = tmp_path / "ref.fa"
        reads = tmp_path / "reads.fq"
        output_bam = tmp_path / "output.bam"

        ref.write_text(">chr1\nACGT\n")
        reads.write_text("@read1\nACGT\n+\nIIII\n")

        # Track temp files created
        temp_files_created = []

        def mock_run_side_effect(cmd, **kwargs):
            stdout_file = kwargs.get("stdout")
            if hasattr(stdout_file, "name"):
                temp_files_created.append(stdout_file.name)
            if hasattr(stdout_file, "write"):
                stdout_file.write("@HD\tVN:1.0\n")
            return Mock(returncode=0, stderr="")

        mocker.patch("subprocess.run", side_effect=mock_run_side_effect)

        def mock_popen_side_effect(cmd, **kwargs):
            if "view" in cmd:
                unsorted_path = Path(str(output_bam) + ".unsorted")
                unsorted_path.write_bytes(b"BAM_DATA")
                temp_files_created.append(str(unsorted_path))
            elif "sort" in cmd:
                output_bam.write_bytes(b"SORTED_BAM")
            elif "index" in cmd:
                Path(str(output_bam) + ".bai").write_bytes(b"INDEX")

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
        align_ont_reads_with_minimap2("minimap2", "samtools", str(ref), str(reads), str(output_bam))

        # Assert: Temp files should be cleaned up
        for temp_file in temp_files_created:
            assert not Path(temp_file).exists(), f"Temp file not cleaned up: {temp_file}"

    def test_handles_minimap2_failure(self, mocker, tmp_path):
        """Test that minimap2 failure is properly caught and reported."""
        # Arrange
        ref = tmp_path / "ref.fa"
        reads = tmp_path / "reads.fq"
        output_bam = tmp_path / "output.bam"

        ref.write_text(">chr1\nACGT\n")
        reads.write_text("@read1\nACGT\n+\nIIII\n")

        # Mock minimap2 failure
        import subprocess

        mocker.patch(
            "subprocess.run",
            side_effect=subprocess.CalledProcessError(1, "minimap2", stderr="Alignment failed"),
        )

        # Act & Assert: Function wraps ExternalToolError in RuntimeError
        with pytest.raises(RuntimeError, match="ONT read alignment failed"):
            align_ont_reads_with_minimap2(
                "minimap2", "samtools", str(ref), str(reads), str(output_bam)
            )

    def test_handles_conda_commands(self, mocker, tmp_path):
        """Test that conda/mamba commands are parsed correctly."""
        # Arrange
        ref = tmp_path / "ref.fa"
        reads = tmp_path / "reads.fq"
        output_bam = tmp_path / "output.bam"

        ref.write_text(">chr1\nACGT\n")
        reads.write_text("@read1\nACGT\n+\nIIII\n")

        run_calls = []

        def mock_run_side_effect(cmd, **kwargs):
            run_calls.append(cmd)
            if hasattr(kwargs.get("stdout"), "write"):
                kwargs["stdout"].write("@HD\tVN:1.0\n")
            return Mock(returncode=0, stderr="")

        mocker.patch("subprocess.run", side_effect=mock_run_side_effect)

        def mock_popen_side_effect(cmd, **kwargs):
            if "view" in cmd:
                Path(str(output_bam) + ".unsorted").write_bytes(b"BAM_DATA")
            elif "sort" in cmd:
                output_bam.write_bytes(b"SORTED_BAM")
            elif "index" in cmd:
                Path(str(output_bam) + ".bai").write_bytes(b"INDEX")

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

        # Act: Use conda command with spaces
        align_ont_reads_with_minimap2(
            minimap2_cmd="conda run -n ont minimap2",
            samtools_cmd="conda run -n ont samtools",
            human_reference=str(ref),
            reads_fastq=str(reads),
            output_bam=str(output_bam),
        )

        # Assert: Conda command should be parsed
        minimap2_cmd = run_calls[0]
        assert "conda" in minimap2_cmd
        assert "run" in minimap2_cmd
        assert "-n" in minimap2_cmd
        assert "ont" in minimap2_cmd
        assert "minimap2" in minimap2_cmd
