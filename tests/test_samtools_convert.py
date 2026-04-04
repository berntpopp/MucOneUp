"""Tests for samtools_convert paired FASTQ conversion."""

from pathlib import Path
from unittest.mock import patch

import pytest

from muc_one_up.read_simulator.utils.common_utils import RunResult
from muc_one_up.read_simulator.wrappers.samtools_convert import (
    FastqConversionOptions,
    convert_bam_to_paired_fastq,
)

FASTQ_RECORD = "@read1\nACGT\n+\nIIII\n"


def _write_fake_outputs(fq1: Path, fq2: Path):
    """Return a side_effect that writes fake FASTQ to the output paths."""

    def side_effect(cmds, *, capture=True, timeout=None, cwd=None):
        fq1.write_text(FASTQ_RECORD)
        fq2.write_text(FASTQ_RECORD)
        return RunResult(returncode=0, stdout="", stderr="", command="mock")

    return side_effect


@pytest.fixture
def bam_and_outputs(tmp_path):
    """Create a fake BAM (non-empty) and output FASTQ paths."""
    bam = tmp_path / "input.bam"
    bam.write_bytes(b"\x1f\x8b" + b"\x00" * 100)  # fake gzip header
    fq1 = tmp_path / "R1.fastq"
    fq2 = tmp_path / "R2.fastq"
    return bam, fq1, fq2


class TestConvertBamToPairedFastq:
    """Exercise the collation pipeline path in convert_bam_to_paired_fastq."""

    def test_collation_pipeline_builds_correct_commands(self, bam_and_outputs):
        """Collation path must call run_pipeline with collate|fastq commands."""
        bam, fq1, fq2 = bam_and_outputs
        captured_cmds = []

        def capture_pipeline(cmds, *, capture=True, timeout=None, cwd=None):
            captured_cmds.extend(cmds)
            fq1.write_text(FASTQ_RECORD)
            fq2.write_text(FASTQ_RECORD)
            return RunResult(returncode=0, stdout="", stderr="", command="mock")

        opts = FastqConversionOptions(
            validate_pairs=False, collate_before_conversion=True, threads=4
        )

        with patch(
            "muc_one_up.read_simulator.wrappers.samtools_convert.run_pipeline",
            side_effect=capture_pipeline,
        ):
            convert_bam_to_paired_fastq("samtools", bam, fq1, fq2, opts)

        assert len(captured_cmds) == 2
        collate_cmd, fastq_cmd = captured_cmds
        assert collate_cmd[0] == "samtools"
        assert collate_cmd[1] == "collate"
        assert str(bam) in collate_cmd
        assert fastq_cmd[0] == "samtools"
        assert fastq_cmd[1] == "fastq"
        assert str(fq1) in fastq_cmd
        assert str(fq2) in fastq_cmd
        # stdin marker must be last arg
        assert fastq_cmd[-1] == "-"

    def test_collation_pipeline_with_conda_wrapper(self, bam_and_outputs):
        """Conda-wrapped samtools must include wrapper prefix in both commands."""
        bam, fq1, fq2 = bam_and_outputs
        captured_cmds = []

        def capture_pipeline(cmds, *, capture=True, timeout=None, cwd=None):
            captured_cmds.extend(cmds)
            fq1.write_text(FASTQ_RECORD)
            fq2.write_text(FASTQ_RECORD)
            return RunResult(returncode=0, stdout="", stderr="", command="mock")

        opts = FastqConversionOptions(validate_pairs=False, collate_before_conversion=True)

        with patch(
            "muc_one_up.read_simulator.wrappers.samtools_convert.run_pipeline",
            side_effect=capture_pipeline,
        ):
            convert_bam_to_paired_fastq("mamba run -n myenv samtools", bam, fq1, fq2, opts)

        collate_cmd, fastq_cmd = captured_cmds
        # Both commands must start with the full wrapper prefix
        assert collate_cmd[:4] == ["mamba", "run", "-n", "myenv"]
        assert collate_cmd[4] == "samtools"
        assert collate_cmd[5] == "collate"
        assert fastq_cmd[:4] == ["mamba", "run", "-n", "myenv"]
        assert fastq_cmd[4] == "samtools"
        assert fastq_cmd[5] == "fastq"

    def test_collation_pipeline_stderr_logging(self, bam_and_outputs):
        """Pipeline stderr should be captured without error."""
        bam, fq1, fq2 = bam_and_outputs

        def pipeline_with_stderr(cmds, *, capture=True, timeout=None, cwd=None):
            fq1.write_text(FASTQ_RECORD)
            fq2.write_text(FASTQ_RECORD)
            return RunResult(
                returncode=0, stdout="", stderr="[W::samtools] warning", command="mock"
            )

        opts = FastqConversionOptions(validate_pairs=False, collate_before_conversion=True)

        with patch(
            "muc_one_up.read_simulator.wrappers.samtools_convert.run_pipeline",
            side_effect=pipeline_with_stderr,
        ):
            result = convert_bam_to_paired_fastq("samtools", bam, fq1, fq2, opts)

        assert result == (str(fq1), str(fq2))

    def test_no_collation_uses_run_command(self, bam_and_outputs):
        """Without collation, should use run_command instead of run_pipeline."""
        bam, fq1, fq2 = bam_and_outputs

        def fake_run_command(cmd, **kwargs):
            fq1.write_text(FASTQ_RECORD)
            fq2.write_text(FASTQ_RECORD)
            return RunResult(returncode=0, stdout="", stderr="", command="mock")

        opts = FastqConversionOptions(validate_pairs=False, collate_before_conversion=False)

        with patch(
            "muc_one_up.read_simulator.wrappers.samtools_convert.run_command",
            side_effect=fake_run_command,
        ):
            result = convert_bam_to_paired_fastq("samtools", bam, fq1, fq2, opts)

        assert result == (str(fq1), str(fq2))
