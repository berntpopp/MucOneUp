"""Tests for fragment_preparation stage (Illumina pipeline stages 1-8).

Following project testing principles:
- Mock at the module level where functions are imported into the stage module
- Test OUR code's logic (argument passing, error handling, result construction)
- NOT testing the underlying wrappers themselves
"""

from __future__ import annotations

from pathlib import Path

import pytest

from muc_one_up.exceptions import ConfigurationError
from muc_one_up.read_simulator.assembly_context import AssemblyContext
from muc_one_up.read_simulator.stages import FragmentResult
from muc_one_up.read_simulator.stages.fragment_preparation import prepare_fragments

MODULE = "muc_one_up.read_simulator.stages.fragment_preparation"


@pytest.fixture
def tools():
    return {
        "reseq": "reseq",
        "faToTwoBit": "faToTwoBit",
        "samtools": "samtools",
        "pblat": "pblat",
    }


@pytest.fixture
def rs_config():
    return {
        "reseq_model": "/path/to/model.reseq",
        "sample_bam": "/path/to/sample.bam",
        "read_number": 1000,
        "fragment_size": 350,
        "fragment_sd": 50,
        "min_fragment": 200,
        "binding_min": 0.5,
        "threads": 4,
        "pblat_threads": 4,
        "pblat_min_score": 95,
        "pblat_min_identity": 95,
        "seqtoillumina_timeout": 120,
    }


@pytest.fixture
def assembly_ctx():
    return AssemblyContext(
        assembly_name="hg38",
        left_constant="ACGT",
        right_constant="TGCA",
        sample_bam="/path/to/assembly_sample.bam",
    )


@pytest.fixture
def assembly_ctx_no_bam():
    return AssemblyContext(
        assembly_name="hg38",
        left_constant="ACGT",
        right_constant="TGCA",
        sample_bam=None,
    )


class TestPrepareFragments:
    """Tests for the prepare_fragments function."""

    def test_returns_fragment_result(self, mocker, tmp_path, tools, rs_config, assembly_ctx):
        """Mock all wrapper calls; verify returns FragmentResult with correct paths."""
        mock_replace_ns = mocker.patch(f"{MODULE}.replace_Ns")
        mock_gen_syser = mocker.patch(f"{MODULE}.generate_systematic_errors")
        mock_fa_to_2bit = mocker.patch(f"{MODULE}.fa_to_twobit")
        mock_extract_subset = mocker.patch(
            f"{MODULE}.extract_subset_reference",
            return_value=str(tmp_path / "_test_collated.bam"),
        )
        mock_run_pblat = mocker.patch(f"{MODULE}.run_pblat")
        mock_simulate = mocker.patch(f"{MODULE}.simulate_fragments")
        mock_create_reads = mocker.patch(f"{MODULE}.create_reads")
        mock_split_reads = mocker.patch(f"{MODULE}.split_reads")

        input_fa = str(tmp_path / "input.fa")
        output_base = "test_sample"

        result = prepare_fragments(
            tools=tools,
            rs_config=rs_config,
            input_fa=input_fa,
            output_dir=tmp_path,
            output_base=output_base,
            assembly_ctx=assembly_ctx,
        )

        assert isinstance(result, FragmentResult)
        assert result.r1_fastq == str(tmp_path / f"{output_base}_R1.fastq.gz")
        assert result.r2_fastq == str(tmp_path / f"{output_base}_R2.fastq.gz")
        assert len(result.intermediate_files) > 0

        # Verify all stage functions were called
        mock_replace_ns.assert_called_once()
        mock_gen_syser.assert_called_once()
        mock_fa_to_2bit.assert_called_once()
        mock_extract_subset.assert_called_once()
        mock_run_pblat.assert_called_once()
        mock_simulate.assert_called_once()
        mock_create_reads.assert_called_once()
        mock_split_reads.assert_called_once()

    def test_intermediate_files_have_underscore_prefix(
        self, mocker, tmp_path, tools, rs_config, assembly_ctx
    ):
        """Intermediate files should use underscore prefix naming convention."""
        mocker.patch(f"{MODULE}.replace_Ns")
        mocker.patch(f"{MODULE}.generate_systematic_errors")
        mocker.patch(f"{MODULE}.fa_to_twobit")
        collated_bam = str(tmp_path / "_test_collated.bam")
        mocker.patch(f"{MODULE}.extract_subset_reference", return_value=collated_bam)
        mocker.patch(f"{MODULE}.run_pblat")
        mocker.patch(f"{MODULE}.simulate_fragments")
        mocker.patch(f"{MODULE}.create_reads")
        mocker.patch(f"{MODULE}.split_reads")

        input_fa = str(tmp_path / "input.fa")
        output_base = "mysample"

        result = prepare_fragments(
            tools=tools,
            rs_config=rs_config,
            input_fa=input_fa,
            output_dir=tmp_path,
            output_base=output_base,
            assembly_ctx=assembly_ctx,
        )

        # All intermediate files with names derived from output_base should use underscore prefix
        named_intermediates = [
            f for f in result.intermediate_files if output_base in Path(f).name
        ]
        for f in named_intermediates:
            assert Path(f).name.startswith("_"), (
                f"Expected underscore prefix for intermediate file: {f}"
            )

    def test_raises_on_missing_reseq_model(self, mocker, tmp_path, tools, assembly_ctx):
        """ConfigurationError is raised when reseq_model is absent from rs_config."""
        mocker.patch(f"{MODULE}.replace_Ns")

        rs_config_no_model: dict = {}

        with pytest.raises(ConfigurationError, match="reseq_model"):
            prepare_fragments(
                tools=tools,
                rs_config=rs_config_no_model,
                input_fa=str(tmp_path / "input.fa"),
                output_dir=tmp_path,
                output_base="test",
                assembly_ctx=assembly_ctx,
            )

    def test_raises_on_missing_sample_bam(
        self, mocker, tmp_path, tools, rs_config, assembly_ctx_no_bam
    ):
        """ConfigurationError raised when no sample_bam in ctx or rs_config."""
        mocker.patch(f"{MODULE}.replace_Ns")
        mocker.patch(f"{MODULE}.generate_systematic_errors")
        mocker.patch(f"{MODULE}.fa_to_twobit")

        # rs_config without sample_bam fallback
        rs_config_no_bam = {k: v for k, v in rs_config.items() if k != "sample_bam"}

        with pytest.raises(ConfigurationError, match=r"[Ss]ample [Bb][Aa][Mm]|sample_bam"):
            prepare_fragments(
                tools=tools,
                rs_config=rs_config_no_bam,
                input_fa=str(tmp_path / "input.fa"),
                output_dir=tmp_path,
                output_base="test",
                assembly_ctx=assembly_ctx_no_bam,
            )

    def test_assembly_ctx_bam_takes_priority_over_rs_config(
        self, mocker, tmp_path, tools, rs_config, assembly_ctx
    ):
        """assembly_ctx.sample_bam is used in preference to rs_config sample_bam."""
        mocker.patch(f"{MODULE}.replace_Ns")
        mocker.patch(f"{MODULE}.generate_systematic_errors")
        mocker.patch(f"{MODULE}.fa_to_twobit")
        mock_extract = mocker.patch(
            f"{MODULE}.extract_subset_reference",
            return_value=str(tmp_path / "_collated.bam"),
        )
        mocker.patch(f"{MODULE}.run_pblat")
        mocker.patch(f"{MODULE}.simulate_fragments")
        mocker.patch(f"{MODULE}.create_reads")
        mocker.patch(f"{MODULE}.split_reads")

        prepare_fragments(
            tools=tools,
            rs_config=rs_config,
            input_fa=str(tmp_path / "input.fa"),
            output_dir=tmp_path,
            output_base="test",
            assembly_ctx=assembly_ctx,
        )

        # The first argument to extract_subset_reference should be the ctx bam
        call_args = mock_extract.call_args
        assert call_args[0][0] == assembly_ctx.sample_bam

    def test_fragment_origins_path_set_when_source_tracker_present(
        self, mocker, tmp_path, tools, rs_config, assembly_ctx
    ):
        """fragment_origins_path is passed to simulate_fragments when source_tracker given."""
        mocker.patch(f"{MODULE}.replace_Ns")
        mocker.patch(f"{MODULE}.generate_systematic_errors")
        mocker.patch(f"{MODULE}.fa_to_twobit")
        mocker.patch(
            f"{MODULE}.extract_subset_reference",
            return_value=str(tmp_path / "_collated.bam"),
        )
        mocker.patch(f"{MODULE}.run_pblat")
        mock_simulate = mocker.patch(f"{MODULE}.simulate_fragments")
        mocker.patch(f"{MODULE}.create_reads")
        mocker.patch(f"{MODULE}.split_reads")

        output_base = "mysample"
        prepare_fragments(
            tools=tools,
            rs_config=rs_config,
            input_fa=str(tmp_path / "input.fa"),
            output_dir=tmp_path,
            output_base=output_base,
            assembly_ctx=assembly_ctx,
            source_tracker=object(),  # any truthy object
        )

        call_kwargs = mock_simulate.call_args[1]
        assert call_kwargs.get("fragment_origins_path") is not None
        expected = str(tmp_path / f"{output_base}_fragment_origins.tsv")
        assert call_kwargs["fragment_origins_path"] == expected

    def test_fragment_origins_path_none_without_source_tracker(
        self, mocker, tmp_path, tools, rs_config, assembly_ctx
    ):
        """fragment_origins_path is None when source_tracker is not provided."""
        mocker.patch(f"{MODULE}.replace_Ns")
        mocker.patch(f"{MODULE}.generate_systematic_errors")
        mocker.patch(f"{MODULE}.fa_to_twobit")
        mocker.patch(
            f"{MODULE}.extract_subset_reference",
            return_value=str(tmp_path / "_collated.bam"),
        )
        mocker.patch(f"{MODULE}.run_pblat")
        mock_simulate = mocker.patch(f"{MODULE}.simulate_fragments")
        mocker.patch(f"{MODULE}.create_reads")
        mocker.patch(f"{MODULE}.split_reads")

        prepare_fragments(
            tools=tools,
            rs_config=rs_config,
            input_fa=str(tmp_path / "input.fa"),
            output_dir=tmp_path,
            output_base="mysample",
            assembly_ctx=assembly_ctx,
            source_tracker=None,
        )

        call_kwargs = mock_simulate.call_args[1]
        assert call_kwargs.get("fragment_origins_path") is None
