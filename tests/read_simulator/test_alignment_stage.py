"""Tests for alignment_stage (Illumina pipeline stages 9-11).

Following project testing principles:
- Mock at the module level where functions are imported into the stage module
- Test OUR code's logic (argument passing, error handling, result construction)
- NOT testing the underlying wrappers themselves
"""

from __future__ import annotations

import pytest

from muc_one_up.exceptions import ConfigurationError
from muc_one_up.read_simulator.assembly_context import AssemblyContext
from muc_one_up.read_simulator.stages import AlignmentResult
from muc_one_up.read_simulator.stages.alignment import align_and_refine

MODULE = "muc_one_up.read_simulator.stages.alignment"


@pytest.fixture
def tools():
    return {
        "bwa": "bwa",
        "samtools": "samtools",
    }


@pytest.fixture
def rs_config_basic():
    """Minimal rs_config: VNTR disabled, no downsampling."""
    return {
        "threads": 4,
        "vntr_capture_efficiency": {"enabled": False},
    }


@pytest.fixture
def assembly_ctx():
    return AssemblyContext(
        assembly_name="hg38",
        left_constant="",
        right_constant="",
        vntr_region="chr1:100000-200000",
    )


class TestAlignAndRefine:
    """Tests for the align_and_refine function."""

    def test_returns_alignment_result_basic(
        self, mocker, tmp_path, tools, rs_config_basic, assembly_ctx
    ):
        """VNTR disabled, no downsampling: returns AlignmentResult with correct final_bam."""
        mock_align = mocker.patch(f"{MODULE}.align_reads")

        output_base = "test_sample"
        result = align_and_refine(
            tools=tools,
            rs_config=rs_config_basic,
            r1=str(tmp_path / "r1.fastq.gz"),
            r2=str(tmp_path / "r2.fastq.gz"),
            human_ref="/ref/hg38.fa",
            output_dir=tmp_path,
            output_base=output_base,
            assembly_ctx=assembly_ctx,
        )

        assert isinstance(result, AlignmentResult)
        assert result.final_bam.endswith(f"{output_base}.bam")
        assert result.intermediate_bams == []

        # align_reads must have been called once
        mock_align.assert_called_once()

    def test_align_reads_called_with_correct_args(
        self, mocker, tmp_path, tools, rs_config_basic, assembly_ctx
    ):
        """align_reads receives r1, r2, human_ref, output_bam, tools, threads."""
        mock_align = mocker.patch(f"{MODULE}.align_reads")

        output_base = "mysample"
        r1 = str(tmp_path / "r1.fastq.gz")
        r2 = str(tmp_path / "r2.fastq.gz")
        human_ref = "/ref/hg38.fa"

        align_and_refine(
            tools=tools,
            rs_config=rs_config_basic,
            r1=r1,
            r2=r2,
            human_ref=human_ref,
            output_dir=tmp_path,
            output_base=output_base,
            assembly_ctx=assembly_ctx,
        )

        call_args, _call_kwargs = mock_align.call_args
        # Positional args: r1, r2, human_ref, output_bam, tools, threads
        assert call_args[0] == r1
        assert call_args[1] == r2
        assert call_args[2] == human_ref
        assert call_args[3].endswith(f"{output_base}.bam")
        assert call_args[4] == tools

    def test_vntr_bias_produces_intermediate_bam(self, mocker, tmp_path, tools, assembly_ctx):
        """VNTR enabled: final_bam contains 'vntr_biased', original BAM is intermediate."""
        mocker.patch(f"{MODULE}.align_reads")
        mocker.patch("shutil.rmtree")

        mock_model = mocker.MagicMock()
        mock_model.apply_efficiency_bias.return_value = {"reads_kept": 100}
        mocker.patch(
            f"{MODULE}.VNTREfficiencyModel",
            return_value=mock_model,
        )

        rs_config_vntr = {
            "threads": 4,
            "vntr_capture_efficiency": {
                "enabled": True,
                "penalty_factor": 0.375,
                "seed": 42,
                "output_fastq": {"enabled": False},
                "validation": {"report_statistics": False},
            },
        }

        output_base = "test_sample"
        result = align_and_refine(
            tools=tools,
            rs_config=rs_config_vntr,
            r1=str(tmp_path / "r1.fastq.gz"),
            r2=str(tmp_path / "r2.fastq.gz"),
            human_ref="/ref/hg38.fa",
            output_dir=tmp_path,
            output_base=output_base,
            assembly_ctx=assembly_ctx,
        )

        assert isinstance(result, AlignmentResult)
        assert "vntr_biased" in result.final_bam
        assert len(result.intermediate_bams) == 1
        # The original aligned BAM should be the intermediate
        assert result.intermediate_bams[0].endswith(f"{output_base}.bam")

    def test_vntr_failure_falls_back_to_original_bam(self, mocker, tmp_path, tools, assembly_ctx):
        """If VNTR model raises, continue with the original aligned BAM."""
        mocker.patch(f"{MODULE}.align_reads")

        mocker.patch(
            f"{MODULE}.VNTREfficiencyModel",
            side_effect=RuntimeError("vntr failed"),
        )

        rs_config_vntr = {
            "threads": 4,
            "vntr_capture_efficiency": {"enabled": True},
        }

        output_base = "test_sample"
        result = align_and_refine(
            tools=tools,
            rs_config=rs_config_vntr,
            r1=str(tmp_path / "r1.fastq.gz"),
            r2=str(tmp_path / "r2.fastq.gz"),
            human_ref="/ref/hg38.fa",
            output_dir=tmp_path,
            output_base=output_base,
            assembly_ctx=assembly_ctx,
        )

        # Falls back: final_bam is the original aligned BAM
        assert result.final_bam.endswith(f"{output_base}.bam")
        assert result.intermediate_bams == []

    def test_downsampling_vntr_mode_adds_depth_file(self, mocker, tmp_path, tools, assembly_ctx):
        """Downsampling in 'vntr' mode adds depth file to intermediate_files."""
        mocker.patch(f"{MODULE}.align_reads")

        depth_file = str(tmp_path / "depth.txt")
        mocker.patch(
            f"{MODULE}.calculate_vntr_coverage",
            side_effect=[
                (150.0, depth_file),          # Pre-downsampling measurement
                (100.0, "/tmp/depth_post.txt"),  # Post-downsampling validation
            ],
        )
        mocker.patch(f"{MODULE}.downsample_bam")

        rs_config_ds = {
            "threads": 4,
            "vntr_capture_efficiency": {"enabled": False},
            "coverage": 100,
            "downsample_mode": "vntr",
        }

        result = align_and_refine(
            tools=tools,
            rs_config=rs_config_ds,
            r1=str(tmp_path / "r1.fastq.gz"),
            r2=str(tmp_path / "r2.fastq.gz"),
            human_ref="/ref/hg38.fa",
            output_dir=tmp_path,
            output_base="test_sample",
            assembly_ctx=assembly_ctx,
        )

        assert depth_file in result.intermediate_files
        assert "downsampled" in result.final_bam

    def test_downsampling_no_downsample_when_coverage_below_target(
        self, mocker, tmp_path, tools, assembly_ctx
    ):
        """No downsampling when current coverage is already below target."""
        mocker.patch(f"{MODULE}.align_reads")

        depth_file = str(tmp_path / "depth.txt")
        mocker.patch(
            f"{MODULE}.calculate_vntr_coverage",
            return_value=(50.0, depth_file),  # below target of 100
        )
        mock_downsample = mocker.patch(f"{MODULE}.downsample_bam")

        rs_config_ds = {
            "threads": 4,
            "vntr_capture_efficiency": {"enabled": False},
            "coverage": 100,
            "downsample_mode": "vntr",
        }

        output_base = "test_sample"
        result = align_and_refine(
            tools=tools,
            rs_config=rs_config_ds,
            r1=str(tmp_path / "r1.fastq.gz"),
            r2=str(tmp_path / "r2.fastq.gz"),
            human_ref="/ref/hg38.fa",
            output_dir=tmp_path,
            output_base=output_base,
            assembly_ctx=assembly_ctx,
        )

        mock_downsample.assert_not_called()
        assert result.final_bam.endswith(f"{output_base}.bam")

    def test_downsampling_non_vntr_mode(self, mocker, tmp_path, tools, assembly_ctx):
        """Downsampling in 'non_vntr' mode calls downsample_entire_bam."""
        mocker.patch(f"{MODULE}.align_reads")

        depth_file = str(tmp_path / "depth.txt")
        mocker.patch(
            f"{MODULE}.calculate_target_coverage",
            side_effect=[
                (200.0, depth_file),             # Pre-downsampling measurement
                (100.0, "/tmp/depth_post.txt"),   # Post-downsampling validation
            ],
        )
        mock_downsample_entire = mocker.patch(f"{MODULE}.downsample_entire_bam")

        rs_config_ds = {
            "threads": 4,
            "vntr_capture_efficiency": {"enabled": False},
            "coverage": 100,
            "downsample_mode": "non_vntr",
            "sample_target_bed": "/path/to/targets.bed",
        }

        result = align_and_refine(
            tools=tools,
            rs_config=rs_config_ds,
            r1=str(tmp_path / "r1.fastq.gz"),
            r2=str(tmp_path / "r2.fastq.gz"),
            human_ref="/ref/hg38.fa",
            output_dir=tmp_path,
            output_base="test_sample",
            assembly_ctx=assembly_ctx,
        )

        mock_downsample_entire.assert_called_once()
        assert "downsampled" in result.final_bam

    def test_raises_configuration_error_invalid_downsample_mode(
        self, mocker, tmp_path, tools, assembly_ctx
    ):
        """ConfigurationError raised when downsample_mode is invalid."""
        mocker.patch(f"{MODULE}.align_reads")

        rs_config_bad = {
            "threads": 4,
            "vntr_capture_efficiency": {"enabled": False},
            "coverage": 100,
            "downsample_mode": "bad_mode",
        }

        with pytest.raises(ConfigurationError, match="bad_mode"):
            align_and_refine(
                tools=tools,
                rs_config=rs_config_bad,
                r1=str(tmp_path / "r1.fastq.gz"),
                r2=str(tmp_path / "r2.fastq.gz"),
                human_ref="/ref/hg38.fa",
                output_dir=tmp_path,
                output_base="test_sample",
                assembly_ctx=assembly_ctx,
            )

    def test_raises_configuration_error_vntr_mode_no_vntr_region(self, mocker, tmp_path, tools):
        """ConfigurationError raised when vntr mode selected but no vntr_region in context."""
        mocker.patch(f"{MODULE}.align_reads")

        assembly_ctx_no_vntr = AssemblyContext(
            assembly_name="hg38",
            left_constant="",
            right_constant="",
            vntr_region=None,
        )

        rs_config_ds = {
            "threads": 4,
            "vntr_capture_efficiency": {"enabled": False},
            "coverage": 100,
            "downsample_mode": "vntr",
        }

        with pytest.raises(ConfigurationError, match=r"[Vv][Nn][Tt][Rr]"):
            align_and_refine(
                tools=tools,
                rs_config=rs_config_ds,
                r1=str(tmp_path / "r1.fastq.gz"),
                r2=str(tmp_path / "r2.fastq.gz"),
                human_ref="/ref/hg38.fa",
                output_dir=tmp_path,
                output_base="test_sample",
                assembly_ctx=assembly_ctx_no_vntr,
            )

    def test_raises_configuration_error_non_vntr_mode_no_bed(
        self, mocker, tmp_path, tools, assembly_ctx
    ):
        """ConfigurationError raised when non_vntr mode but no sample_target_bed."""
        mocker.patch(f"{MODULE}.align_reads")

        rs_config_ds = {
            "threads": 4,
            "vntr_capture_efficiency": {"enabled": False},
            "coverage": 100,
            "downsample_mode": "non_vntr",
            # no sample_target_bed
        }

        with pytest.raises(ConfigurationError, match=r"sample_target_bed|[Bb][Ee][Dd]"):
            align_and_refine(
                tools=tools,
                rs_config=rs_config_ds,
                r1=str(tmp_path / "r1.fastq.gz"),
                r2=str(tmp_path / "r2.fastq.gz"),
                human_ref="/ref/hg38.fa",
                output_dir=tmp_path,
                output_base="test_sample",
                assembly_ctx=assembly_ctx,
            )

    def test_post_downsampling_validation_logged(self, mocker, tmp_path, tools, assembly_ctx):
        """After downsampling, VNTR coverage is re-measured and logged."""
        mocker.patch(f"{MODULE}.align_reads")

        # Return coverages: first call = pre-downsample (200x), second call = post-downsample (155x)
        mock_calc_vntr = mocker.patch(
            f"{MODULE}.calculate_vntr_coverage",
            side_effect=[(200.0, "/tmp/depth1.txt"), (155.0, "/tmp/depth2.txt")],
        )
        mocker.patch(f"{MODULE}.downsample_bam")

        rs_config = {
            "threads": 4,
            "coverage": 150,
            "downsample_mode": "vntr",
            "downsample_seed": 42,
            "vntr_capture_efficiency": {"enabled": False},
        }
        assembly_ctx_with_vntr = mocker.MagicMock()
        assembly_ctx_with_vntr.vntr_region = "chr1:155188487-155192239"
        assembly_ctx_with_vntr.assembly_name = "hg38"

        result = align_and_refine(
            tools=tools,
            rs_config=rs_config,
            r1=str(tmp_path / "r1.fastq.gz"),
            r2=str(tmp_path / "r2.fastq.gz"),
            human_ref="/ref/hg38.fa",
            output_dir=tmp_path,
            output_base="test_sample",
            assembly_ctx=assembly_ctx_with_vntr,
        )

        # calculate_vntr_coverage called twice: pre-downsample + post-validation
        assert mock_calc_vntr.call_count == 2

    def test_post_downsampling_validation_warns_on_deviation(
        self, mocker, tmp_path, tools, assembly_ctx, caplog
    ):
        """Warning logged when post-downsampling coverage deviates >20% from target."""
        import logging

        mocker.patch(f"{MODULE}.align_reads")
        # Pre-downsample: 500x, Post-downsample: 250x (67% deviation from target 150)
        mocker.patch(
            f"{MODULE}.calculate_vntr_coverage",
            side_effect=[(500.0, "/tmp/d1.txt"), (250.0, "/tmp/d2.txt")],
        )
        mocker.patch(f"{MODULE}.downsample_bam")

        rs_config = {
            "threads": 4,
            "coverage": 150,
            "downsample_mode": "vntr",
            "downsample_seed": 42,
            "vntr_capture_efficiency": {"enabled": False},
        }
        assembly_ctx_with_vntr = mocker.MagicMock()
        assembly_ctx_with_vntr.vntr_region = "chr1:155188487-155192239"
        assembly_ctx_with_vntr.assembly_name = "hg38"

        with caplog.at_level(logging.WARNING):
            align_and_refine(
                tools=tools,
                rs_config=rs_config,
                r1=str(tmp_path / "r1.fastq.gz"),
                r2=str(tmp_path / "r2.fastq.gz"),
                human_ref="/ref/hg38.fa",
                output_dir=tmp_path,
                output_base="test_sample",
                assembly_ctx=assembly_ctx_with_vntr,
            )

        assert any("deviates" in record.message for record in caplog.records)
