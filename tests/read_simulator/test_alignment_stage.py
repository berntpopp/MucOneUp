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
                "penalty_factor": 0.39,
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

    def test_downsampling_vntr_mode_two_step(self, mocker, tmp_path, tools, assembly_ctx):
        """VNTR mode uses two-step: downsample entire BAM then VNTR region."""
        mocker.patch(f"{MODULE}.align_reads")

        # Step 1: non-VNTR BED coverage (200x → fraction to hit 100x target)
        # Step 2: VNTR coverage after step 1, then post-validation
        target_depth = str(tmp_path / "target_depth.txt")
        mocker.patch(
            f"{MODULE}.calculate_target_coverage",
            side_effect=[
                (200.0, target_depth),  # Step 1: measure non-VNTR BED
                (100.0, "/tmp/non_vntr_post1.txt"),  # Step 1: re-measure after downsample
                (99.0, "/tmp/post_flank.txt"),  # Post-validation: final non-VNTR
            ],
        )
        mocker.patch(
            f"{MODULE}.calculate_vntr_coverage",
            side_effect=[
                (500.0, "/tmp/vntr_post1.txt"),  # Step 2: measure VNTR after step 1
                (140.0, "/tmp/post_vntr.txt"),  # Post-validation: final VNTR
            ],
        )
        mocker.patch(f"{MODULE}.downsample_entire_bam")  # Step 1
        mocker.patch(f"{MODULE}.downsample_bam")  # Step 2

        rs_config_ds = {
            "threads": 4,
            "vntr_capture_efficiency": {"enabled": False},
            "coverage": 100,
            "downsample_mode": "vntr",
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

        assert target_depth in result.intermediate_files
        assert "downsampled" in result.final_bam

    def test_downsampling_vntr_mode_skips_when_below_target(
        self, mocker, tmp_path, tools, assembly_ctx
    ):
        """No downsampling when non-VNTR BED coverage is already below target."""
        mocker.patch(f"{MODULE}.align_reads")

        depth_file = str(tmp_path / "depth.txt")
        mocker.patch(
            f"{MODULE}.calculate_target_coverage",
            side_effect=[
                (50.0, depth_file),  # Step 1: non-VNTR BED below target (skip)
                (50.0, "/tmp/re_flank.txt"),  # Re-measure non-VNTR for ratio calc
                (50.0, "/tmp/post_flank.txt"),  # Post-validation non-VNTR
            ],
        )
        mocker.patch(
            f"{MODULE}.calculate_vntr_coverage",
            side_effect=[
                (60.0, "/tmp/vntr.txt"),  # Step 2: VNTR measurement (60 < 50*1.4=70, skip)
                (60.0, "/tmp/post_vntr.txt"),  # Post-validation VNTR
            ],
        )
        mock_ds_entire = mocker.patch(f"{MODULE}.downsample_entire_bam")
        mock_ds_bam = mocker.patch(f"{MODULE}.downsample_bam")

        rs_config_ds = {
            "threads": 4,
            "vntr_capture_efficiency": {"enabled": False},
            "coverage": 100,
            "downsample_mode": "vntr",
            "sample_target_bed": "/path/to/targets.bed",
        }

        output_base = "test_sample"
        align_and_refine(
            tools=tools,
            rs_config=rs_config_ds,
            r1=str(tmp_path / "r1.fastq.gz"),
            r2=str(tmp_path / "r2.fastq.gz"),
            human_ref="/ref/hg38.fa",
            output_dir=tmp_path,
            output_base=output_base,
            assembly_ctx=assembly_ctx,
        )

        mock_ds_entire.assert_not_called()
        mock_ds_bam.assert_not_called()

    def test_downsampling_non_vntr_mode(self, mocker, tmp_path, tools, assembly_ctx):
        """Downsampling in 'non_vntr' mode calls downsample_entire_bam."""
        mocker.patch(f"{MODULE}.align_reads")

        depth_file = str(tmp_path / "depth.txt")
        mocker.patch(
            f"{MODULE}.calculate_target_coverage",
            side_effect=[
                (200.0, depth_file),  # Pre-downsampling measurement
                (100.0, "/tmp/depth_post.txt"),  # Post-downsampling validation
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

    def test_post_downsampling_vntr_mode_logs_ratio(self, mocker, tmp_path, tools, assembly_ctx):
        """VNTR mode logs final VNTR and non-VNTR coverage after both steps."""
        mocker.patch(f"{MODULE}.align_reads")

        # Step 1: non-VNTR BED 300x → downsample entire BAM
        # Step 2: VNTR 400x after step 1 → downsample VNTR region
        # Post-validation: measure both
        mocker.patch(
            f"{MODULE}.calculate_target_coverage",
            side_effect=[
                (300.0, "/tmp/t1.txt"),  # Step 1: measure non-VNTR
                (150.0, "/tmp/t2.txt"),  # Step 2: re-measure non-VNTR after step 1
                (149.0, "/tmp/t3.txt"),  # Post-validation: final non-VNTR
            ],
        )
        mock_calc_vntr = mocker.patch(
            f"{MODULE}.calculate_vntr_coverage",
            side_effect=[
                (400.0, "/tmp/v1.txt"),  # Step 2: VNTR after step 1
                (200.0, "/tmp/v2.txt"),  # Post-validation: final VNTR
            ],
        )
        mocker.patch(f"{MODULE}.downsample_entire_bam")
        mocker.patch(f"{MODULE}.downsample_bam")

        rs_config = {
            "threads": 4,
            "coverage": 150,
            "downsample_mode": "vntr",
            "downsample_seed": 42,
            "sample_target_bed": "/path/to/targets.bed",
            "vntr_capture_efficiency": {"enabled": False},
        }
        assembly_ctx_with_vntr = mocker.MagicMock()
        assembly_ctx_with_vntr.vntr_region = "chr1:155188487-155192239"
        assembly_ctx_with_vntr.assembly_name = "hg38"

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

        # calculate_vntr_coverage called: step 2 measure + post-validation
        assert mock_calc_vntr.call_count == 2

    def test_post_downsampling_validation_warns_on_deviation(
        self, mocker, tmp_path, tools, assembly_ctx, caplog
    ):
        """Warning logged when post-downsampling coverage deviates >20% from target (non_vntr mode)."""
        import logging

        mocker.patch(f"{MODULE}.align_reads")
        mocker.patch(
            f"{MODULE}.calculate_target_coverage",
            side_effect=[(500.0, "/tmp/d1.txt"), (250.0, "/tmp/d2.txt")],
        )
        mocker.patch(f"{MODULE}.downsample_entire_bam")

        rs_config = {
            "threads": 4,
            "coverage": 150,
            "downsample_mode": "non_vntr",
            "downsample_seed": 42,
            "sample_target_bed": "/path/to/targets.bed",
            "vntr_capture_efficiency": {"enabled": False},
        }

        with caplog.at_level(logging.WARNING):
            align_and_refine(
                tools=tools,
                rs_config=rs_config,
                r1=str(tmp_path / "r1.fastq.gz"),
                r2=str(tmp_path / "r2.fastq.gz"),
                human_ref="/ref/hg38.fa",
                output_dir=tmp_path,
                output_base="test_sample",
                assembly_ctx=assembly_ctx,
            )

        assert any("deviates" in record.message for record in caplog.records)

    def test_post_downsampling_validation_non_vntr_mode(
        self, mocker, tmp_path, tools, assembly_ctx
    ):
        """Post-downsampling validation also works in non_vntr mode."""
        mocker.patch(f"{MODULE}.align_reads")

        mock_calc_target = mocker.patch(
            f"{MODULE}.calculate_target_coverage",
            side_effect=[(200.0, "/tmp/depth1.txt"), (105.0, "/tmp/depth2.txt")],
        )
        mocker.patch(f"{MODULE}.downsample_entire_bam")

        rs_config = {
            "threads": 4,
            "coverage": 100,
            "downsample_mode": "non_vntr",
            "downsample_seed": 42,
            "sample_target_bed": "/path/to/targets.bed",
            "vntr_capture_efficiency": {"enabled": False},
        }

        align_and_refine(
            tools=tools,
            rs_config=rs_config,
            r1=str(tmp_path / "r1.fastq.gz"),
            r2=str(tmp_path / "r2.fastq.gz"),
            human_ref="/ref/hg38.fa",
            output_dir=tmp_path,
            output_base="test_sample",
            assembly_ctx=assembly_ctx,
        )

        # calculate_target_coverage called twice: pre-downsample + post-validation
        assert mock_calc_target.call_count == 2
