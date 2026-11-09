"""Tests for shared pipeline_utils module.

Tests cover:
- resolve_pipeline_outputs: output dir/base resolution with fallback chain
- resolve_human_reference: single precedence for the alignment reference
- create_pipeline_metadata: thin wrapper around write_metadata_file
- cleanup_intermediates: safe removal of intermediate files
"""

from __future__ import annotations

import logging
from datetime import datetime
from unittest.mock import patch

import pytest

from muc_one_up.exceptions import ConfigurationError
from muc_one_up.read_simulator.assembly_context import AssemblyContext
from muc_one_up.read_simulator.output_config import OutputConfig
from muc_one_up.read_simulator.pipeline_utils import (
    cleanup_intermediates,
    create_pipeline_metadata,
    intermediate_work_dir,
    resolve_human_reference,
    resolve_pipeline_outputs,
)


class TestResolvePipelineOutputs:
    """Tests for resolve_pipeline_outputs function."""

    def test_with_output_config(self, tmp_path):
        """OutputConfig takes precedence over rs_config and input file parent."""
        out_dir = tmp_path / "explicit_out"
        out_dir.mkdir()
        output_config = OutputConfig(out_dir=out_dir, out_base="my_base")

        input_fa = str(tmp_path / "input" / "sample.fa")
        rs_config = {"output_dir": str(tmp_path / "rs_out")}

        resolved_dir, resolved_base = resolve_pipeline_outputs(
            input_fa=input_fa,
            rs_config=rs_config,
            output_config=output_config,
        )

        assert resolved_dir == out_dir
        assert resolved_base == "my_base"

    def test_without_output_config_uses_rs_config(self, tmp_path):
        """Falls back to rs_config['output_dir'] when output_config is None."""
        rs_out = tmp_path / "rs_out"
        input_fa = str(tmp_path / "input" / "sample.fa")
        rs_config = {"output_dir": str(rs_out)}

        resolved_dir, resolved_base = resolve_pipeline_outputs(
            input_fa=input_fa,
            rs_config=rs_config,
            output_config=None,
        )

        assert resolved_dir == rs_out
        # base name should be derived from input file stem
        assert resolved_base == "sample"

    def test_without_output_config_or_rs_config(self, tmp_path):
        """Falls back to input file parent when neither output_config nor rs_config dir."""
        parent_dir = tmp_path / "input"
        parent_dir.mkdir()
        input_fa = str(parent_dir / "sample.fa")

        resolved_dir, resolved_base = resolve_pipeline_outputs(
            input_fa=input_fa,
            rs_config={},
            output_config=None,
        )

        assert resolved_dir == parent_dir
        assert resolved_base == "sample"

    def test_creates_output_directory(self, tmp_path):
        """mkdir is called for the resolved output directory if it doesn't exist."""
        new_dir = tmp_path / "new_subdir" / "deeply_nested"
        assert not new_dir.exists()

        output_config = OutputConfig(out_dir=new_dir, out_base="base")
        input_fa = str(tmp_path / "sample.fa")

        resolve_pipeline_outputs(
            input_fa=input_fa,
            rs_config={},
            output_config=output_config,
        )

        assert new_dir.exists()


def _ctx(human_reference=None):
    return AssemblyContext(
        assembly_name="hg38",
        left_constant="ACGT",
        right_constant="TGCA",
        human_reference=human_reference,
    )


class TestResolveHumanReference:
    """One precedence: human_reference_<asm> > human_reference > reference_genomes."""

    def _genomes(self, tmp_path):
        ref = tmp_path / "hg38.fa"
        ref.write_text(">chr1\nACGT\n")
        return ref, {
            "reference_assembly": "hg38",
            "reference_genomes": {"hg38": {"fasta_path": str(ref)}},
        }

    def test_explicit_human_reference_wins_over_reference_genomes(self, tmp_path):
        _, config = self._genomes(tmp_path)
        explicit = tmp_path / "explicit.fa"
        explicit.write_text(">c\nA\n")
        config["read_simulation"] = {"human_reference": str(explicit)}

        assert resolve_human_reference(config, _ctx(), "bwa") == str(explicit)

    def test_assembly_keyed_reference_wins(self, tmp_path):
        _, config = self._genomes(tmp_path)
        asm = tmp_path / "asm.fa"
        asm.write_text(">c\nA\n")
        config["read_simulation"] = {"human_reference": "/explicit.fa"}

        assert resolve_human_reference(config, _ctx(str(asm)), "bwa") == str(asm)

    def test_missing_explicit_reference_raises(self, tmp_path):
        _, config = self._genomes(tmp_path)
        config["read_simulation"] = {"human_reference": str(tmp_path / "gone.fa")}

        with pytest.raises(ConfigurationError, match=r"not found.*gone\.fa"):
            resolve_human_reference(config, _ctx(), "bwa")

    def test_falls_back_to_reference_genomes(self, tmp_path):
        ref, config = self._genomes(tmp_path)

        assert resolve_human_reference(config, _ctx(), "minimap2") == str(ref)

    def test_broken_reference_genomes_entry_raises_with_reason(self, tmp_path):
        config = {"reference_genomes": {"hg38": {"fasta_path": str(tmp_path / "typo.fa")}}}

        with pytest.raises(ConfigurationError, match=r"typo\.fa"):
            resolve_human_reference(config, _ctx(), "bwa")

    def test_raises_when_no_reference(self):
        with pytest.raises(ConfigurationError, match="No human reference"):
            resolve_human_reference({"reference_assembly": "hg38"}, _ctx(), "bwa")


class TestCreatePipelineMetadata:
    """Tests for create_pipeline_metadata function."""

    def test_writes_metadata_file(self, tmp_path):
        """Delegates to write_metadata_file with correct arguments including Path conversion."""
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        start_time = datetime(2026, 1, 1, 10, 0, 0)
        end_time = datetime(2026, 1, 1, 11, 0, 0)
        config = {"tools": {}}
        tools_used = ["bwa", "samtools"]

        with patch(
            "muc_one_up.read_simulator.pipeline_utils.write_metadata_file",
            return_value=str(output_dir / "base_metadata.tsv"),
        ) as mock_write:
            result = create_pipeline_metadata(
                output_dir=output_dir,
                output_base="base",
                config=config,
                start_time=start_time,
                end_time=end_time,
                platform="Illumina",
                tools_used=tools_used,
            )

        mock_write.assert_called_once_with(
            str(output_dir),
            "base",
            config,
            start_time,
            end_time,
            "Illumina",
            tools_used,
            extra_rows=None,
            seed=None,
        )
        assert result == str(output_dir / "base_metadata.tsv")


class TestCleanupIntermediates:
    """Tests for cleanup_intermediates function."""

    def test_removes_existing_files(self, tmp_path):
        """Files that exist are removed."""
        f1 = tmp_path / "file1.tmp"
        f2 = tmp_path / "file2.tmp"
        f1.write_text("data1")
        f2.write_text("data2")

        cleanup_intermediates([str(f1), str(f2)])

        assert not f1.exists()
        assert not f2.exists()

    def test_ignores_missing_files(self, tmp_path):
        """No error is raised for nonexistent files."""
        missing = str(tmp_path / "does_not_exist.tmp")

        # Should not raise
        cleanup_intermediates([missing])

    def test_ignores_none_and_empty(self, tmp_path):
        """None and empty string entries are silently skipped."""
        existing = tmp_path / "real.tmp"
        existing.write_text("data")

        # Should not raise
        cleanup_intermediates([None, "", str(existing), None])

        assert not existing.exists()

    def test_logs_warning_on_failure(self, tmp_path, caplog):
        """A warning is logged (not an exception) when unlink raises an error."""
        target = tmp_path / "locked.tmp"
        target.write_text("data")

        with (
            patch("pathlib.Path.unlink", side_effect=OSError("permission denied")),
            caplog.at_level(logging.WARNING),
        ):
            cleanup_intermediates([str(target)])

        assert any("permission denied" in record.message for record in caplog.records)


class TestCleanupUnlessKept:
    """keep_intermediate_files is honoured uniformly by amplicon pipelines (#110)."""

    def test_removes_by_default(self, tmp_path):
        from muc_one_up.read_simulator.pipeline_utils import cleanup_unless_kept

        f = tmp_path / "x.tmp"
        f.write_text("x")
        cleanup_unless_kept({}, [str(f)])
        assert not f.exists()

    def test_keeps_when_configured(self, tmp_path):
        from muc_one_up.read_simulator.pipeline_utils import cleanup_unless_kept

        f = tmp_path / "x.tmp"
        f.write_text("x")
        cleanup_unless_kept({"read_simulation": {"keep_intermediate_files": True}}, [str(f)])
        assert f.exists()


@pytest.mark.parametrize("keep", [False, True])
def test_intermediate_work_dir_kept_under_output_dir(tmp_path, keep):
    config = {"read_simulation": {"keep_intermediate_files": keep}}

    with intermediate_work_dir(config, tmp_path, "sim_") as work_dir:
        (work_dir / "x").write_text("x")

    assert (work_dir.parent == tmp_path) is keep
    assert work_dir.exists() is keep
