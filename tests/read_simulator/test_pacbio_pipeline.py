"""Integration tests for PacBio HiFi read simulation pipeline.

Following Phase 2 testing principles:
- Mock ONLY at wrapper boundaries (pbsim3, ccs, minimap2, samtools wrappers)
- Test OUR pipeline orchestration logic
- NOT testing individual tools themselves

These tests verify that the pipeline correctly:
1. Orchestrates the complete workflow (pbsim3 → CCS → FASTQ → minimap2)
2. Passes parameters correctly between stages
3. Handles intermediate file cleanup
4. Supports both alignment and FASTQ-only modes
5. Validates configuration parameters
"""

import pytest

from muc_one_up.read_simulator.pacbio_pipeline import simulate_pacbio_hifi_reads


class TestPacBioHiFiPipeline:
    """Test simulate_pacbio_hifi_reads pipeline orchestration."""

    def test_full_pipeline_with_alignment(self, mocker, tmp_path):
        """Test complete pipeline: pbsim3 → CCS → FASTQ → minimap2 alignment."""
        # Arrange
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGTACGTACGT\n")

        human_ref = tmp_path / "hg38.fa"
        human_ref.write_text(">chr1\nACGTACGTACGTACGTACGT\n")

        config = {
            "tools": {"samtools": "samtools", "minimap2": "minimap2"},
            "pacbio_params": {
                "pbsim3_cmd": "pbsim",
                "ccs_cmd": "ccs",
                "model_type": "qshmm",
                "model_file": str(tmp_path / "model.model"),
                "coverage": 30,
                "pass_num": 3,
                "min_passes": 3,
                "min_rq": 0.99,
                "threads": 4,
                "seed": 42,
            },
        }

        # Create model file
        (tmp_path / "model.model").write_text("model data")

        # Expected intermediate files
        clr_bam = tmp_path / "input_clr.bam"
        hifi_bam = tmp_path / "input_hifi.bam"
        hifi_fastq = tmp_path / "input_hifi.fastq"
        aligned_bam = tmp_path / "input_aligned.bam"

        # Mock all wrapper functions
        mock_pbsim3 = mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.run_pbsim3_simulation"
        )
        mock_pbsim3.return_value = str(clr_bam)

        def create_clr_bam(*args, **kwargs):
            clr_bam.write_bytes(b"CLR BAM data" * 100)
            return str(clr_bam)

        mock_pbsim3.side_effect = create_clr_bam

        mock_ccs = mocker.patch("muc_one_up.read_simulator.pacbio_pipeline.run_ccs_consensus")
        mock_ccs.return_value = str(hifi_bam)

        def create_hifi_bam(*args, **kwargs):
            hifi_bam.write_bytes(b"HiFi BAM data" * 100)
            return str(hifi_bam)

        mock_ccs.side_effect = create_hifi_bam

        mock_convert_fastq = mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.convert_bam_to_fastq"
        )
        mock_convert_fastq.return_value = str(hifi_fastq)

        def create_hifi_fastq(*args, **kwargs):
            hifi_fastq.write_text("@read1\nACGT\n+\nIIII\n")
            return str(hifi_fastq)

        mock_convert_fastq.side_effect = create_hifi_fastq

        mock_align = mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.align_reads_with_minimap2"
        )
        mock_align.return_value = str(aligned_bam)

        def create_aligned_bam(*args, **kwargs):
            aligned_bam.write_bytes(b"Aligned BAM data" * 100)
            return str(aligned_bam)

        mock_align.side_effect = create_aligned_bam

        # Mock cleanup
        mock_cleanup = mocker.patch("muc_one_up.read_simulator.pacbio_pipeline.cleanup_files")

        # Act
        result = simulate_pacbio_hifi_reads(
            config=config, input_fa=str(input_fa), human_reference=str(human_ref)
        )

        # Assert: Verify all stages were called
        assert mock_pbsim3.called
        assert mock_ccs.called
        assert mock_convert_fastq.called
        assert mock_align.called
        assert mock_cleanup.called

        # Verify pbsim3 was called with correct parameters
        pbsim3_call = mock_pbsim3.call_args
        assert pbsim3_call[1]["pbsim3_cmd"] == "pbsim"
        assert pbsim3_call[1]["reference"] == str(input_fa)
        assert pbsim3_call[1]["model_type"] == "qshmm"
        assert pbsim3_call[1]["coverage"] == 30
        assert pbsim3_call[1]["pass_num"] == 3
        assert pbsim3_call[1]["seed"] == 42

        # Verify CCS was called with correct parameters
        ccs_call = mock_ccs.call_args
        assert ccs_call[1]["ccs_cmd"] == "ccs"
        assert ccs_call[1]["input_bam"] == str(clr_bam)
        assert ccs_call[1]["min_passes"] == 3
        assert ccs_call[1]["min_rq"] == 0.99
        assert ccs_call[1]["seed"] == 43  # Offset by 1

        # Verify alignment was called with HiFi preset
        align_call = mock_align.call_args
        assert align_call[1]["reference"] == str(human_ref)
        assert align_call[1]["reads_fastq"] == str(hifi_fastq)
        assert align_call[1]["preset"] == "map-hifi"

        # Verify result is aligned BAM
        assert result == str(aligned_bam)

    def test_pipeline_without_alignment(self, mocker, tmp_path):
        """Test pipeline stops at FASTQ when no human_reference provided."""
        # Arrange
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGTACGTACGT\n")

        config = {
            "tools": {"samtools": "samtools"},
            "pacbio_params": {
                "pbsim3_cmd": "pbsim",
                "ccs_cmd": "ccs",
                "model_type": "qshmm",
                "model_file": str(tmp_path / "model.model"),
                "coverage": 30,
                "pass_num": 3,
                "min_passes": 3,
                "min_rq": 0.99,
            },
        }

        # Create model file
        (tmp_path / "model.model").write_text("model data")

        # Expected intermediate files
        clr_bam = tmp_path / "input_clr.bam"
        hifi_bam = tmp_path / "input_hifi.bam"
        hifi_fastq = tmp_path / "input_hifi.fastq"

        # Mock wrapper functions
        mock_pbsim3 = mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.run_pbsim3_simulation"
        )

        def create_clr_bam(*args, **kwargs):
            clr_bam.write_bytes(b"CLR BAM data" * 100)
            return str(clr_bam)

        mock_pbsim3.side_effect = create_clr_bam

        mock_ccs = mocker.patch("muc_one_up.read_simulator.pacbio_pipeline.run_ccs_consensus")

        def create_hifi_bam(*args, **kwargs):
            hifi_bam.write_bytes(b"HiFi BAM data" * 100)
            return str(hifi_bam)

        mock_ccs.side_effect = create_hifi_bam

        mock_convert_fastq = mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.convert_bam_to_fastq"
        )

        def create_hifi_fastq(*args, **kwargs):
            hifi_fastq.write_text("@read1\nACGT\n+\nIIII\n")
            return str(hifi_fastq)

        mock_convert_fastq.side_effect = create_hifi_fastq

        mock_align = mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.align_reads_with_minimap2"
        )

        mock_cleanup = mocker.patch("muc_one_up.read_simulator.pacbio_pipeline.cleanup_files")

        # Act - No human_reference provided
        result = simulate_pacbio_hifi_reads(
            config=config, input_fa=str(input_fa), human_reference=None
        )

        # Assert: Verify pbsim3, CCS, and FASTQ conversion were called
        assert mock_pbsim3.called
        assert mock_ccs.called
        assert mock_convert_fastq.called

        # Verify alignment was NOT called (no human reference)
        assert not mock_align.called

        # Verify cleanup was called
        assert mock_cleanup.called

        # Verify result is FASTQ (not aligned BAM)
        assert result == str(hifi_fastq)

    def test_validates_missing_parameters(self, tmp_path):
        """Test that missing required parameters raise RuntimeError."""
        # Arrange
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGTACGTACGT\n")

        # Missing model_type
        config = {
            "tools": {"samtools": "samtools"},
            "pacbio_params": {
                "pbsim3_cmd": "pbsim",
                "ccs_cmd": "ccs",
                # Missing model_type!
                "model_file": str(tmp_path / "model.model"),
                "coverage": 30,
                "pass_num": 3,
                "min_passes": 3,
                "min_rq": 0.99,
            },
        }

        # Act & Assert
        with pytest.raises(RuntimeError, match="Missing required PacBio parameters"):
            simulate_pacbio_hifi_reads(config=config, input_fa=str(input_fa))

    def test_validates_parameter_ranges(self, mocker, tmp_path):
        """Test that invalid parameter ranges raise RuntimeError."""
        # Arrange
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGTACGTACGT\n")

        config = {
            "tools": {"samtools": "samtools"},
            "pacbio_params": {
                "pbsim3_cmd": "pbsim",
                "ccs_cmd": "ccs",
                "model_type": "qshmm",
                "model_file": str(tmp_path / "model.model"),
                "coverage": 30,
                "pass_num": 1,  # Invalid! Must be ≥ 2
                "min_passes": 3,
                "min_rq": 0.99,
            },
        }

        # Create model file
        (tmp_path / "model.model").write_text("model data")

        # Act & Assert
        with pytest.raises(RuntimeError, match="Invalid PacBio parameters"):
            simulate_pacbio_hifi_reads(config=config, input_fa=str(input_fa))

    def test_uses_default_parameters(self, mocker, tmp_path):
        """Test that default parameters are used when not specified."""
        # Arrange
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGTACGTACGT\n")

        config = {
            "tools": {"samtools": "samtools"},
            "pacbio_params": {
                "pbsim3_cmd": "pbsim",
                "ccs_cmd": "ccs",
                "model_type": "qshmm",
                "model_file": str(tmp_path / "model.model"),
                "coverage": 30,
                "pass_num": 3,
                "min_passes": 3,
                "min_rq": 0.99,
                # No accuracy/length parameters specified - should use defaults
            },
        }

        # Create model file
        (tmp_path / "model.model").write_text("model data")

        clr_bam = tmp_path / "input_clr.bam"
        hifi_bam = tmp_path / "input_hifi.bam"
        hifi_fastq = tmp_path / "input_hifi.fastq"

        # Mock wrappers
        mock_pbsim3 = mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.run_pbsim3_simulation"
        )

        def create_clr_bam(*args, **kwargs):
            clr_bam.write_bytes(b"CLR BAM data" * 100)
            return str(clr_bam)

        mock_pbsim3.side_effect = create_clr_bam

        mock_ccs = mocker.patch("muc_one_up.read_simulator.pacbio_pipeline.run_ccs_consensus")

        def create_hifi_bam(*args, **kwargs):
            hifi_bam.write_bytes(b"HiFi BAM data" * 100)
            return str(hifi_bam)

        mock_ccs.side_effect = create_hifi_bam

        mock_convert_fastq = mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.convert_bam_to_fastq"
        )

        def create_hifi_fastq(*args, **kwargs):
            hifi_fastq.write_text("@read1\nACGT\n+\nIIII\n")
            return str(hifi_fastq)

        mock_convert_fastq.side_effect = create_hifi_fastq

        mocker.patch("muc_one_up.read_simulator.pacbio_pipeline.cleanup_files")

        # Act
        simulate_pacbio_hifi_reads(config=config, input_fa=str(input_fa), human_reference=None)

        # Assert: Verify pbsim3 was called with default parameters
        from muc_one_up.read_simulator.constants import (
            DEFAULT_PBSIM3_ACCURACY_MEAN,
            DEFAULT_PBSIM3_ACCURACY_MIN,
            DEFAULT_PBSIM3_ACCURACY_SD,
            DEFAULT_PBSIM3_LENGTH_MAX,
            DEFAULT_PBSIM3_LENGTH_MEAN,
            DEFAULT_PBSIM3_LENGTH_MIN,
            DEFAULT_PBSIM3_LENGTH_SD,
        )

        pbsim3_call = mock_pbsim3.call_args
        assert pbsim3_call[1]["accuracy_mean"] == DEFAULT_PBSIM3_ACCURACY_MEAN
        assert pbsim3_call[1]["accuracy_sd"] == DEFAULT_PBSIM3_ACCURACY_SD
        assert pbsim3_call[1]["accuracy_min"] == DEFAULT_PBSIM3_ACCURACY_MIN
        assert pbsim3_call[1]["length_mean"] == DEFAULT_PBSIM3_LENGTH_MEAN
        assert pbsim3_call[1]["length_sd"] == DEFAULT_PBSIM3_LENGTH_SD
        assert pbsim3_call[1]["length_min"] == DEFAULT_PBSIM3_LENGTH_MIN
        assert pbsim3_call[1]["length_max"] == DEFAULT_PBSIM3_LENGTH_MAX

    def test_derives_output_paths_from_input(self, mocker, tmp_path):
        """Test that output paths are correctly derived from input filename."""
        # Arrange
        input_fa = tmp_path / "sample.001.simulated.fa"
        input_fa.write_text(">chr1\nACGTACGTACGT\n")

        config = {
            "tools": {"samtools": "samtools"},
            "pacbio_params": {
                "pbsim3_cmd": "pbsim",
                "ccs_cmd": "ccs",
                "model_type": "qshmm",
                "model_file": str(tmp_path / "model.model"),
                "coverage": 30,
                "pass_num": 3,
                "min_passes": 3,
                "min_rq": 0.99,
            },
        }

        # Create model file
        (tmp_path / "model.model").write_text("model data")

        # Expected paths based on input stem
        expected_clr_prefix = tmp_path / "sample.001.simulated_clr"
        expected_hifi_bam = tmp_path / "sample.001.simulated_hifi.bam"
        expected_hifi_fastq = tmp_path / "sample.001.simulated_hifi.fastq"

        # Mock wrappers
        mock_pbsim3 = mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.run_pbsim3_simulation"
        )

        def create_clr_bam(*args, **kwargs):
            clr_bam = tmp_path / "sample.001.simulated_clr.bam"
            clr_bam.write_bytes(b"CLR BAM data" * 100)
            return str(clr_bam)

        mock_pbsim3.side_effect = create_clr_bam

        mock_ccs = mocker.patch("muc_one_up.read_simulator.pacbio_pipeline.run_ccs_consensus")

        def create_hifi_bam(*args, **kwargs):
            expected_hifi_bam.write_bytes(b"HiFi BAM data" * 100)
            return str(expected_hifi_bam)

        mock_ccs.side_effect = create_hifi_bam

        mock_convert_fastq = mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.convert_bam_to_fastq"
        )

        def create_hifi_fastq(*args, **kwargs):
            expected_hifi_fastq.write_text("@read1\nACGT\n+\nIIII\n")
            return str(expected_hifi_fastq)

        mock_convert_fastq.side_effect = create_hifi_fastq

        mocker.patch("muc_one_up.read_simulator.pacbio_pipeline.cleanup_files")

        # Act
        simulate_pacbio_hifi_reads(config=config, input_fa=str(input_fa), human_reference=None)

        # Assert: Verify output_prefix uses correct base name
        pbsim3_call = mock_pbsim3.call_args
        assert pbsim3_call[1]["output_prefix"] == str(expected_clr_prefix)

        # Verify CCS output_bam uses correct name
        ccs_call = mock_ccs.call_args
        assert ccs_call[1]["output_bam"] == str(expected_hifi_bam)

        # Verify FASTQ output uses correct name
        fastq_call = mock_convert_fastq.call_args
        assert fastq_call[1]["output_fastq"] == str(expected_hifi_fastq)
