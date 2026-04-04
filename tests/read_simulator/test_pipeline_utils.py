"""Tests for shared pipeline_utils module.

Tests cover:
- resolve_pipeline_outputs: output dir/base resolution with fallback chain
- resolve_human_reference: 3-path fallback for human reference resolution
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


class TestResolveHumanReference:
    """Tests for resolve_human_reference function."""

    def test_from_reference_genomes(self, tmp_path):
        """Uses reference_genomes section as first priority."""
        ref_file = tmp_path / "hg38.fa"
        ref_file.write_text(">chr1\nACGT\n")
        # Create required index files so validation passes
        for ext in [".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"]:
            (tmp_path / f"hg38.fa{ext}").touch()

        config = {
            "reference_assembly": "hg38",
            "reference_genomes": {
                "hg38": {
                    "fasta_path": str(ref_file),
                    "vntr_region": "chr1:155188487-155192239",
                }
            },
        }
        assembly_ctx = AssemblyContext(
            assembly_name="hg38",
            left_constant="ACGT",
            right_constant="TGCA",
            human_reference=str(tmp_path / "fallback.fa"),
        )

        result = resolve_human_reference(config=config, assembly_ctx=assembly_ctx, aligner="bwa")

        assert result == str(ref_file)

    def test_fallback_to_assembly_context(self, tmp_path):
        """Falls back to assembly_ctx.human_reference when reference_genomes lookup fails."""
        fallback_ref = tmp_path / "assembly_ctx_ref.fa"
        fallback_ref.write_text(">chr1\nACGT\n")

        config = {
            "reference_assembly": "hg38",
            # No reference_genomes section — will cause get_reference_path_for_assembly to raise
        }
        assembly_ctx = AssemblyContext(
            assembly_name="hg38",
            left_constant="ACGT",
            right_constant="TGCA",
            human_reference=str(fallback_ref),
        )

        # Mock get_reference_path_for_assembly to simulate failure
        with patch(
            "muc_one_up.read_simulator.pipeline_utils.get_reference_path_for_assembly",
            side_effect=Exception("no reference_genomes"),
        ):
            result = resolve_human_reference(
                config=config, assembly_ctx=assembly_ctx, aligner="bwa"
            )

        assert result == str(fallback_ref)

    def test_fallback_to_rs_config(self, tmp_path):
        """Falls back to config['read_simulation']['human_reference'] (legacy)."""
        legacy_ref = tmp_path / "legacy_ref.fa"
        legacy_ref.write_text(">chr1\nACGT\n")

        config = {
            "reference_assembly": "hg38",
            "read_simulation": {
                "human_reference": str(legacy_ref),
            },
            # No reference_genomes
        }
        assembly_ctx = AssemblyContext(
            assembly_name="hg38",
            left_constant="ACGT",
            right_constant="TGCA",
            human_reference=None,  # No assembly_ctx fallback
        )

        with patch(
            "muc_one_up.read_simulator.pipeline_utils.get_reference_path_for_assembly",
            side_effect=Exception("no reference_genomes"),
        ):
            result = resolve_human_reference(
                config=config, assembly_ctx=assembly_ctx, aligner="bwa"
            )

        assert result == str(legacy_ref)

    def test_raises_when_no_reference(self, tmp_path):
        """Raises ConfigurationError when no reference can be found via any path."""
        config = {
            "reference_assembly": "hg38",
            # No reference_genomes, no read_simulation.human_reference
        }
        assembly_ctx = AssemblyContext(
            assembly_name="hg38",
            left_constant="ACGT",
            right_constant="TGCA",
            human_reference=None,
        )

        with (
            patch(
                "muc_one_up.read_simulator.pipeline_utils.get_reference_path_for_assembly",
                side_effect=Exception("no reference_genomes"),
            ),
            pytest.raises(ConfigurationError),
        ):
            resolve_human_reference(config=config, assembly_ctx=assembly_ctx, aligner="bwa")


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
