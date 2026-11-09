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
from pathlib import Path
from unittest.mock import Mock

import pytest

from muc_one_up.exceptions import ExternalToolError, FileOperationError
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

    def test_constructs_command_with_seed(self, mocker, tmp_path, tools_dict):
        """Test that replaceN command includes --seed when provided."""
        input_fa = tmp_path / "input.fa"
        output_fa = tmp_path / "output.fa"
        input_fa.write_text(">seq1\nACGTNNACGT\n")

        mock_proc = Mock()
        mock_proc.returncode = 0
        mock_proc.stdout = Mock()
        mock_proc.stderr = Mock()
        mock_proc.stdout.readline = Mock(return_value=b"")
        mock_proc.stderr.readline = Mock(return_value=b"")
        mock_proc.wait = Mock(return_value=0)
        mock_proc.pid = 12345

        mock_popen = mocker.patch("subprocess.Popen", return_value=mock_proc)

        replace_Ns(str(input_fa), str(output_fa), tools_dict, seed=42)

        cmd = mock_popen.call_args[0][0]
        assert "--seed" in cmd
        assert cmd[cmd.index("--seed") + 1] == "42"


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

    def test_constructs_illuminaPE_command_with_seed(self, mocker, tmp_path, tools_dict):
        """Test that illuminaPE command includes --seed when provided."""
        input_fa = tmp_path / "input.fa"
        model_file = tmp_path / "model.txt"
        output_fq = tmp_path / "output.fq"
        input_fa.write_text(">seq1\nACGT\n")
        model_file.write_text("MODEL_DATA")

        def mock_popen_side_effect(cmd, **kwargs):
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

        generate_systematic_errors(
            str(input_fa), str(model_file), str(output_fq), tools_dict, seed=123
        )

        cmd = mock_popen.call_args[0][0]
        assert "--seed" in cmd
        assert cmd[cmd.index("--seed") + 1] == "123"

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
    """Test create_reads command construction and output validation."""

    @pytest.fixture(autouse=True)
    def _reseq2(self, mocker):
        """Report ReSeq2 so the version probe neither runs nor triggers chunking."""
        mocker.patch(
            "muc_one_up.read_simulator.wrappers.reseq_wrapper.get_tool_version",
            return_value="ReSeq2 version 2.0.3",
        )

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

    def test_constructs_seqToIllumina_command_with_seed(self, mocker, tmp_path, tools_dict):
        """Test that seqToIllumina command includes --seed when provided."""
        fragments = tmp_path / "fragments.fa"
        model_file = tmp_path / "model.txt"
        output_fq = tmp_path / "output.fq"
        fragments.write_text(">frag1\nACGT\n")
        model_file.write_text("MODEL_DATA")

        def mock_popen_side_effect(cmd, **kwargs):
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

        create_reads(str(fragments), str(model_file), str(output_fq), 4, tools_dict, seed=456)

        cmd = mock_popen.call_args[0][0]
        assert "--seed" in cmd
        assert cmd[cmd.index("--seed") + 1] == "456"

    def test_raises_error_on_empty_output(self, mocker, tmp_path, tools_dict):
        """Test that empty output after successful run raises FileOperationError."""
        fragments = tmp_path / "fragments.fa"
        model_file = tmp_path / "model.txt"
        output_fq = tmp_path / "output.fq"

        fragments.write_text(">frag1\nACGT\n")
        model_file.write_text("MODEL_DATA")

        def mock_popen_side_effect(cmd, **kwargs):
            # Command succeeds but produces empty output
            output_fq.write_text("")
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

        with pytest.raises(FileOperationError, match="missing or empty"):
            create_reads(str(fragments), str(model_file), str(output_fq), 4, tools_dict)

    def test_raises_error_on_command_failure(self, mocker, tmp_path, tools_dict):
        """Test that command failure raises ExternalToolError."""
        fragments = tmp_path / "fragments.fa"
        model_file = tmp_path / "model.txt"
        output_fq = tmp_path / "output.fq"

        fragments.write_text(">frag1\nACGT\n")
        model_file.write_text("MODEL_DATA")

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

        with pytest.raises(ExternalToolError):
            create_reads(str(fragments), str(model_file), str(output_fq), 4, tools_dict)


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


class TestSplitReadsValidation:
    """split_reads must not silently put R1/R2 out of step (audit a)."""

    @pytest.mark.parametrize(
        "content",
        [
            "@r1/1\nACGT\n+\nIIII\n",  # odd number of records
            "@r1/1\nACGT\n+\nIIII\n@r1/2\nTGCA\n",  # truncated R2
            "@r1/1\nACGT\n+\nIIII\n@r2/2\nTGCA\n+\nIIII\n",  # mismatched mate
        ],
        ids=["odd", "truncated", "mismatch"],
    )
    def test_rejects_broken_interleaved_input(self, tmp_path, content):
        interleaved = tmp_path / "in.fq"
        interleaved.write_text(content)

        with pytest.raises(FileOperationError):
            split_reads(str(interleaved), str(tmp_path / "R1.fq.gz"), str(tmp_path / "R2.fq.gz"))

    def test_accepts_identical_names_with_comments(self, tmp_path):
        interleaved = tmp_path / "in.fq"
        interleaved.write_text("@r1 a\nACGT\n+\nIIII\n@r1 b\nTGCA\n+\nIIII")

        split_reads(str(interleaved), str(tmp_path / "R1.fq.gz"), str(tmp_path / "R2.fq.gz"))

        with gzip.open(tmp_path / "R2.fq.gz", "rt") as f:
            assert f.read() == "@r1 b\nTGCA\n+\nIIII"


class TestCreateReadsChunking:
    """reseq 1.1 hangs on >10,000 fragments; create_reads must chunk the input."""

    def test_chunks_large_input_with_distinct_seeds(self, mocker, tmp_path, tools_dict):
        from muc_one_up.read_simulator.wrappers import reseq_wrapper

        mocker.patch.object(reseq_wrapper, "RESEQ_MAX_FRAGMENTS", 2)
        mocker.patch.object(reseq_wrapper, "get_tool_version", return_value="ReSeq version 1.1")
        fragments = tmp_path / "fragments.fa"
        fragments.write_text("".join(f">{i} x\nACGT\n" for i in range(5)))
        output_fq = tmp_path / "reads.fq"
        calls = []

        def fake_run(cmd, **kwargs):
            args = [str(c) for c in cmd]
            chunk_in = args[args.index("-i") + 1]
            names = [
                ln[1:].split()[0]
                for ln in Path(chunk_in).read_text().splitlines()
                if ln.startswith(">")
            ]
            calls.append((names, args[args.index("--seed") + 1]))
            Path(args[args.index("-o") + 1]).write_text(
                "".join(f"@{n}\nACGT\n+\nIIII\n" for n in names)
            )

        mocker.patch.object(reseq_wrapper, "run_command", side_effect=fake_run)

        create_reads(str(fragments), "model", str(output_fq), 2, tools_dict, seed=10)

        assert calls == [(["0", "1"], "10"), (["2", "3"], "11"), (["4"], "12")]
        assert [ln for ln in output_fq.read_text().splitlines() if ln.startswith("@")] == [
            "@0",
            "@1",
            "@2",
            "@3",
            "@4",
        ]
        assert sorted(p.name for p in tmp_path.iterdir()) == ["fragments.fa", "reads.fq"]

    def test_reseq2_gets_whole_input_in_one_call(self, mocker, tmp_path, tools_dict):
        from muc_one_up.read_simulator.wrappers import reseq_wrapper

        mocker.patch.object(reseq_wrapper, "RESEQ_MAX_FRAGMENTS", 2)
        mocker.patch.object(reseq_wrapper, "get_tool_version", return_value="ReSeq2 version 2.0.3")
        fragments = tmp_path / "fragments.fa"
        fragments.write_text("".join(f">{i} x\nACGT\n" for i in range(5)))
        output_fq = tmp_path / "reads.fq"
        calls = []

        def fake_run(cmd, **kwargs):
            args = [str(c) for c in cmd]
            calls.append((args[args.index("-i") + 1], args[args.index("--seed") + 1]))
            Path(args[args.index("-o") + 1]).write_text("@0\nACGT\n+\nIIII\n")

        mocker.patch.object(reseq_wrapper, "run_command", side_effect=fake_run)

        create_reads(str(fragments), "model", str(output_fq), 2, tools_dict, seed=10)

        assert calls == [(str(fragments), "10")]

    @pytest.mark.parametrize(
        ("version", "legacy"),
        [
            ("ReSeq version 1.1", True),
            ("ReSeq2 version 2.0.3", False),
            ("ReSeq2 version 10.0", False),
            ("N/A", True),
        ],
    )
    def test_is_legacy_reseq(self, mocker, version, legacy):
        from muc_one_up.read_simulator.wrappers import reseq_wrapper

        mocker.patch.object(reseq_wrapper, "get_tool_version", return_value=version)
        assert reseq_wrapper._is_legacy_reseq("reseq") is legacy
