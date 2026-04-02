"""Tests for BWA wrapper - focuses on OUR command construction logic.

Following Phase 2 testing principles:
- Mock ONLY at system boundary (run_command, run_pipeline)
- Test OUR code's logic (command construction, error handling, file validation)
- NOT testing external tools (BWA, samtools) - that's integration testing

These tests verify that our wrapper correctly:
1. Constructs BWA mem commands
2. Chains BWA -> samtools view using run_pipeline
3. Sorts and indexes BAM files
4. Validates input files
5. Handles tool failures appropriately
"""

from pathlib import Path

import pytest

from muc_one_up.exceptions import ExternalToolError, FileOperationError
from muc_one_up.read_simulator.utils.common_utils import RunResult
from muc_one_up.read_simulator.wrappers.bwa_wrapper import align_reads

# Common mock target paths
_PIPELINE = "muc_one_up.read_simulator.wrappers.bwa_wrapper.run_pipeline"
_RUN_CMD = "muc_one_up.read_simulator.wrappers.bwa_wrapper.run_command"


def _ok_result(cmd: str = "ok") -> RunResult:
    return RunResult(returncode=0, stdout="", stderr="", command=cmd)


class TestBWAWrapperCommandConstruction:
    """Test BWA wrapper command construction (tests OUR code)."""

    def test_constructs_correct_bwa_mem_command(self, mocker, tmp_path, tools_dict):
        """Test that BWA mem command is constructed with correct arguments."""
        ref = tmp_path / "ref.fa"
        r1 = tmp_path / "R1.fq"
        r2 = tmp_path / "R2.fq"
        out = tmp_path / "out.bam"

        ref.write_text(">chr1\nACGT")
        r1.write_text("@read1\nACGT\n+\nIIII")
        r2.write_text("@read1\nACGT\n+\nIIII")

        mock_pipeline = mocker.patch(_PIPELINE, return_value=_ok_result())

        # Mock run_command for sort and index - create output files
        def mock_run_cmd(cmd, **kwargs):
            if "sort" in cmd:
                out.write_bytes(b"SORTED_BAM_DATA")
            elif "index" in cmd:
                Path(str(out) + ".bai").write_bytes(b"INDEX_DATA")
            return _ok_result()

        mocker.patch(_RUN_CMD, side_effect=mock_run_cmd)

        align_reads(
            read1=str(r1),
            read2=str(r2),
            human_reference=str(ref),
            output_bam=str(out),
            tools=tools_dict,
            threads=4,
        )

        # Verify run_pipeline was called with correct BWA mem command
        assert mock_pipeline.call_count == 1
        pipeline_cmds = mock_pipeline.call_args[0][0]
        bwa_cmd = pipeline_cmds[0]
        assert bwa_cmd == ["bwa", "mem", "-t", "4", str(ref), str(r1), str(r2)]

    def test_constructs_correct_samtools_view_command(self, mocker, tmp_path, tools_dict):
        """Test that samtools view command is constructed correctly."""
        ref = tmp_path / "ref.fa"
        r1 = tmp_path / "R1.fq"
        r2 = tmp_path / "R2.fq"
        out = tmp_path / "out.bam"

        ref.write_text(">chr1\nACGT")
        r1.write_text("@read1\nACGT\n+\nIIII")
        r2.write_text("@read1\nACGT\n+\nIIII")

        mock_pipeline = mocker.patch(_PIPELINE, return_value=_ok_result())

        def mock_run_cmd(cmd, **kwargs):
            if "sort" in cmd:
                out.write_bytes(b"SORTED_BAM_DATA")
            elif "index" in cmd:
                Path(str(out) + ".bai").write_bytes(b"INDEX_DATA")
            return _ok_result()

        mocker.patch(_RUN_CMD, side_effect=mock_run_cmd)

        align_reads(str(r1), str(r2), str(ref), str(out), tools_dict, threads=4)

        # Verify samtools view command includes -o for direct file output
        pipeline_cmds = mock_pipeline.call_args[0][0]
        samtools_cmd = pipeline_cmds[1]

        unsorted_bam = str(out.parent / f"{out.stem}_unsorted.bam")
        assert samtools_cmd == ["samtools", "view", "-bS", "-o", unsorted_bam, "-"]

    def test_uses_correct_thread_count(self, mocker, tmp_path, tools_dict):
        """Test that thread count is correctly passed to BWA and samtools sort."""
        ref = tmp_path / "ref.fa"
        r1 = tmp_path / "R1.fq"
        r2 = tmp_path / "R2.fq"
        out = tmp_path / "out.bam"

        ref.write_text(">chr1\nACGT")
        r1.write_text("@read1\nACGT\n+\nIIII")
        r2.write_text("@read1\nACGT\n+\nIIII")

        mock_pipeline = mocker.patch(_PIPELINE, return_value=_ok_result())

        run_cmd_calls = []

        def mock_run_cmd(cmd, **kwargs):
            run_cmd_calls.append(cmd)
            if "sort" in cmd:
                out.write_bytes(b"SORTED_BAM_DATA")
            elif "index" in cmd:
                Path(str(out) + ".bai").write_bytes(b"INDEX_DATA")
            return _ok_result()

        mocker.patch(_RUN_CMD, side_effect=mock_run_cmd)

        align_reads(str(r1), str(r2), str(ref), str(out), tools_dict, threads=8)

        # BWA command should have -t 8
        bwa_cmd = mock_pipeline.call_args[0][0][0]
        assert "-t" in bwa_cmd
        t_index = bwa_cmd.index("-t")
        assert bwa_cmd[t_index + 1] == "8"

        # samtools sort should have -@ 8
        sort_cmd = run_cmd_calls[0]
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

    def test_raises_error_when_pipeline_fails(self, mocker, tmp_path, tools_dict):
        """Test that pipeline failure (bwa or samtools view) is properly reported."""
        ref = tmp_path / "ref.fa"
        r1 = tmp_path / "R1.fq"
        r2 = tmp_path / "R2.fq"

        ref.write_text(">chr1\nACGT")
        r1.write_text("@read1\nACGT\n+\nIIII")
        r2.write_text("@read1\nACGT\n+\nIIII")

        mocker.patch(
            _PIPELINE,
            side_effect=ExternalToolError(
                tool="pipeline",
                exit_code=1,
                stderr="Pipeline stage 0 (bwa mem) failed with exit code 1",
                cmd="bwa mem | samtools view",
            ),
        )

        with pytest.raises(ExternalToolError, match=r"pipeline.*failed"):
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
        ref = tmp_path / "ref.fa"
        r1 = tmp_path / "R1.fq"
        r2 = tmp_path / "R2.fq"

        ref.write_text(">chr1\nACGT")
        r1.write_text("@read1\nACGT\n+\nIIII")
        r2.write_text("@read1\nACGT\n+\nIIII")

        mocker.patch(_PIPELINE, return_value=_ok_result())

        # run_command for sort raises ExternalToolError
        mocker.patch(
            _RUN_CMD,
            side_effect=ExternalToolError(
                tool="command",
                exit_code=1,
                stderr="sort failed",
                cmd="samtools sort",
            ),
        )

        with pytest.raises(ExternalToolError, match=r"command.*failed"):
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
        ref = tmp_path / "ref.fa"
        r1 = tmp_path / "R1.fq"
        r2 = tmp_path / "R2.fq"
        out = tmp_path / "out.bam"

        ref.write_text(">chr1\nACGT")
        r1.write_text("@read1\nACGT\n+\nIIII")
        r2.write_text("@read1\nACGT\n+\nIIII")

        mocker.patch(_PIPELINE, return_value=_ok_result())

        def mock_run_cmd(cmd, **kwargs):
            if "sort" in cmd:
                out.write_bytes(b"SORTED_BAM_DATA")
            elif "index" in cmd:
                Path(str(out) + ".bai").write_bytes(b"INDEX_DATA")
            return _ok_result()

        mocker.patch(_RUN_CMD, side_effect=mock_run_cmd)

        align_reads(str(r1), str(r2), str(ref), str(out), tools_dict, threads=4)

        assert out.exists()
        assert Path(str(out) + ".bai").exists()

    def test_raises_error_when_output_bam_missing(self, mocker, tmp_path, tools_dict):
        """Test that missing output BAM raises FileOperationError."""
        ref = tmp_path / "ref.fa"
        r1 = tmp_path / "R1.fq"
        r2 = tmp_path / "R2.fq"
        out = tmp_path / "out.bam"

        ref.write_text(">chr1\nACGT")
        r1.write_text("@read1\nACGT\n+\nIIII")
        r2.write_text("@read1\nACGT\n+\nIIII")

        mocker.patch(_PIPELINE, return_value=_ok_result())

        # run_command succeeds but doesn't create output files
        mocker.patch(_RUN_CMD, return_value=_ok_result())

        with pytest.raises(FileOperationError, match="output file missing or empty"):
            align_reads(str(r1), str(r2), str(ref), str(out), tools_dict, threads=4)
