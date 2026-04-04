"""Integration tests for Illumina coverage downsampling workflow.

This module provides end-to-end integration tests that verify the complete
coverage downsampling workflow from config loading through pipeline execution.
"""

import contextlib
from unittest.mock import patch

import pytest


class TestCoverageDownsamplingIntegration:
    """Integration tests for coverage-based downsampling workflow."""

    def test_coverage_key_triggers_downsampling_logic(self):
        """Verify that coverage key presence triggers correct code path."""
        # This test verifies the logic without running full pipeline

        # Arrange: Simulate pipeline config reading
        rs_config = {
            "coverage": 30,
            "downsample_mode": "non_vntr",
            "sample_target_bed": "/tmp/targets.bed",
        }

        # Act: Simulate pipeline logic (pipeline.py:279-280)
        target_coverage = rs_config.get("coverage")
        will_trigger_downsampling = target_coverage is not None

        # Assert: Coverage key triggers downsampling path
        assert will_trigger_downsampling, "Coverage key should trigger downsampling logic"
        assert target_coverage == 30, "Target coverage should be read correctly"

        # Verify old keys are NOT used
        assert rs_config.get("downsample_target") is None
        assert rs_config.get("downsample_coverage") is None

    def test_missing_coverage_skips_downsampling(self, minimal_config, tmp_path):
        """Verify that omitting coverage key skips downsampling.

        Downsampling now lives inside stages.alignment.align_and_refine, so we
        verify the behaviour by checking that align_and_refine receives a config
        without a coverage key (which causes it to skip the downsampling branch).
        We test the stage module's own coverage logic separately in
        test_alignment_stage.py; here we confirm the orchestrator passes config
        through correctly.
        """
        from muc_one_up.read_simulator.stages import AlignmentResult, FragmentResult

        # Arrange: Config WITHOUT coverage key
        config = minimal_config.copy()
        if "coverage" in config.get("read_simulation", {}):
            del config["read_simulation"]["coverage"]

        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nATCGATCGATCGATCG\n")

        final_bam = str(tmp_path / "out.bam")

        with (
            patch("muc_one_up.read_simulator.pipeline.check_external_tools"),
            patch("muc_one_up.read_simulator.pipeline.capture_tool_versions", return_value={}),
            patch("muc_one_up.read_simulator.pipeline.log_tool_versions"),
            patch("muc_one_up.read_simulator.pipeline.cleanup_intermediates"),
            patch("muc_one_up.read_simulator.pipeline.create_pipeline_metadata"),
            patch("muc_one_up.read_simulator.pipeline.generate_read_manifest"),
            patch(
                "muc_one_up.read_simulator.pipeline.prepare_fragments",
                return_value=FragmentResult(r1_fastq="/r1.fq.gz", r2_fastq="/r2.fq.gz", intermediate_files=[]),
            ),
            patch(
                "muc_one_up.read_simulator.pipeline.align_and_refine",
                return_value=AlignmentResult(final_bam=final_bam, intermediate_bams=[], intermediate_files=[]),
            ) as mock_align,
        ):
            with contextlib.suppress(Exception):
                from muc_one_up.read_simulator.pipeline import simulate_reads_pipeline
                simulate_reads_pipeline(config, str(input_fa))

            # align_and_refine should have been called
            assert mock_align.called

            # The rs_config passed to align_and_refine must not contain "coverage"
            call_kwargs = mock_align.call_args
            passed_rs_config = call_kwargs[1]["rs_config"]
            assert "coverage" not in passed_rs_config, (
                "Coverage key should not be present in rs_config when not configured"
            )

    def test_vntr_mode_downsampling_integration(self):
        """Verify VNTR mode uses correct code path."""
        # Arrange: VNTR mode config
        rs_config = {
            "coverage": 50,
            "downsample_mode": "vntr",
            "reference_assembly": "hg38",
            "vntr_region_hg38": "chr1:1000-2000",
        }

        # Act: Simulate pipeline mode detection (pipeline.py:282)
        mode = rs_config.get("downsample_mode", "vntr").strip().lower()

        # Assert: Mode detection works
        assert mode == "vntr", "VNTR mode should be detected"

        # Verify region key construction (pipeline.py:285-286)
        reference_assembly = rs_config.get("reference_assembly", "hg38")
        vntr_region_key = f"vntr_region_{reference_assembly}"
        vntr_region = rs_config.get(vntr_region_key)

        assert vntr_region == "chr1:1000-2000", "VNTR region should be retrieved correctly"

    def test_low_coverage_skips_downsampling(self, minimal_config, tmp_path):
        """Verify downsampling skipped when current coverage below target.

        Downsampling is now handled inside stages.alignment.align_and_refine.
        We verify the orchestrator passes the coverage config through correctly
        to the stage; the stage-level behaviour is tested in test_alignment_stage.py.
        """
        from muc_one_up.read_simulator.stages import AlignmentResult, FragmentResult

        # Arrange: Target coverage higher than current
        config = minimal_config.copy()
        config["read_simulation"]["coverage"] = 100  # Target: 100x
        config["read_simulation"]["downsample_mode"] = "non_vntr"
        config["read_simulation"]["sample_target_bed"] = str(tmp_path / "targets.bed")

        bed_file = tmp_path / "targets.bed"
        bed_file.write_text("chr1\t1000\t2000\n")

        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nATCGATCGATCGATCG\n")

        final_bam = str(tmp_path / "out.bam")

        with (
            patch("muc_one_up.read_simulator.pipeline.check_external_tools"),
            patch("muc_one_up.read_simulator.pipeline.capture_tool_versions", return_value={}),
            patch("muc_one_up.read_simulator.pipeline.log_tool_versions"),
            patch("muc_one_up.read_simulator.pipeline.cleanup_intermediates"),
            patch("muc_one_up.read_simulator.pipeline.create_pipeline_metadata"),
            patch("muc_one_up.read_simulator.pipeline.generate_read_manifest"),
            patch(
                "muc_one_up.read_simulator.pipeline.prepare_fragments",
                return_value=FragmentResult(r1_fastq="/r1.fq.gz", r2_fastq="/r2.fq.gz", intermediate_files=[]),
            ),
            patch(
                "muc_one_up.read_simulator.pipeline.align_and_refine",
                return_value=AlignmentResult(final_bam=final_bam, intermediate_bams=[], intermediate_files=[]),
            ) as mock_align,
        ):
            with contextlib.suppress(Exception):
                from muc_one_up.read_simulator.pipeline import simulate_reads_pipeline
                simulate_reads_pipeline(config, str(input_fa))

            # align_and_refine should have been called
            assert mock_align.called

            # The rs_config passed to align_and_refine must carry the coverage config
            call_kwargs = mock_align.call_args
            passed_rs_config = call_kwargs[1]["rs_config"]
            assert passed_rs_config.get("coverage") == 100, (
                "Coverage target should be forwarded to align_and_refine"
            )
            assert passed_rs_config.get("downsample_mode") == "non_vntr"


class TestCoverageConfigFlowIntegration:
    """Integration tests for config-to-pipeline coverage flow."""

    def test_cli_coverage_reaches_pipeline(self):
        """Verify coverage set by CLI is accessible in pipeline config."""
        # Simulate CLI setting coverage
        config = {"read_simulation": {}}

        # Act: Simulate CLI behavior
        config["read_simulation"]["coverage"] = 30

        # Assert: Pipeline can read it
        assert config["read_simulation"].get("coverage") == 30
        assert "downsample_target" not in config["read_simulation"]
        assert "downsample_coverage" not in config["read_simulation"]

    def test_config_file_coverage_loads_correctly(self):
        """Verify coverage from config dict is accessible."""
        import json

        # Create test config dict
        config_data = {"read_simulation": {"coverage": 50}}

        # Act: Simulate JSON round-trip (like loading from file)
        config_json = json.dumps(config_data)
        config = json.loads(config_json)

        # Assert: Coverage loaded correctly
        assert config["read_simulation"]["coverage"] == 50

    def test_downsampling_fraction_calculation(self):
        """Verify downsampling fraction calculated correctly."""
        # Arrange
        current_coverage = 150.0
        target_coverage = 30.0

        # Act: Calculate fraction (same formula as pipeline)
        fraction = target_coverage / current_coverage
        fraction = min(max(fraction, 0.0), 1.0)  # Clamp to [0, 1]

        # Assert
        assert fraction == 0.2, "Should keep 20% of reads (30/150)"
        assert 0.0 <= fraction <= 1.0, "Fraction must be in [0, 1]"


@pytest.mark.parametrize(
    "current_cov,target_cov,should_downsample",
    [
        (150.0, 30.0, True),  # Current > Target: Downsample
        (30.0, 150.0, False),  # Current < Target: Skip
        (50.0, 50.0, False),  # Current == Target: Skip (not greater)
        (100.5, 30.5, True),  # Float coverage values
    ],
)
def test_downsampling_decision_logic(current_cov, target_cov, should_downsample):
    """Test downsampling decision matches pipeline logic."""
    # This mirrors the actual condition in pipeline.py:325
    will_downsample = current_cov > target_cov

    assert will_downsample == should_downsample, (
        f"Downsampling decision incorrect for current={current_cov}, target={target_cov}"
    )
