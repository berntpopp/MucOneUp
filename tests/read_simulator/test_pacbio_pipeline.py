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
            "tools": {
                "pbsim3": "pbsim",
                "ccs": "ccs",
                "samtools": "samtools",
                "minimap2": "minimap2",
            },
            "pacbio_params": {
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
        mock_pbsim3.return_value = [str(clr_bam)]

        def create_clr_bam(*args, **kwargs):
            clr_bam.write_bytes(b"CLR BAM data" * 100)
            return [str(clr_bam)]

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
        mock_cleanup = mocker.patch("muc_one_up.read_simulator.pacbio_pipeline.cleanup_unless_kept")

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
            "tools": {"pbsim3": "pbsim", "ccs": "ccs", "samtools": "samtools"},
            "pacbio_params": {
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
            return [str(clr_bam)]

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

        mock_cleanup = mocker.patch("muc_one_up.read_simulator.pacbio_pipeline.cleanup_unless_kept")

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

    def test_diploid_reference_divides_pbsim3_coverage(self, mocker, tmp_path):
        """Test that coverage passed to pbsim3 is divided by number of sequences."""
        input_fa = tmp_path / "input_diploid.fa"
        input_fa.write_text(">hap1\nACGTACGTACGT\n>hap2\nACGTACGTACGT\n")

        config = {
            "tools": {"pbsim3": "pbsim", "ccs": "ccs", "samtools": "samtools"},
            "pacbio_params": {
                "model_type": "qshmm",
                "model_file": str(tmp_path / "model.model"),
                "coverage": 30,
                "pass_num": 3,
                "min_passes": 3,
                "min_rq": 0.99,
            },
        }
        (tmp_path / "model.model").write_text("model data")
        clr_bam = tmp_path / "input_diploid_clr.bam"
        hifi_bam = tmp_path / "input_diploid_hifi.bam"
        hifi_fastq = tmp_path / "input_diploid_hifi.fastq"

        mock_pbsim3 = mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.run_pbsim3_simulation"
        )
        mock_pbsim3.return_value = [str(clr_bam)]
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.run_ccs_consensus",
            return_value=str(hifi_bam),
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.convert_bam_to_fastq",
            return_value=str(hifi_fastq),
        )
        mocker.patch("muc_one_up.read_simulator.pacbio_pipeline.cleanup_unless_kept")

        simulate_pacbio_hifi_reads(config=config, input_fa=str(input_fa), human_reference=None)

        assert mock_pbsim3.called
        assert mock_pbsim3.call_args[1]["coverage"] == 15.0

    def test_pipeline_without_alignment_executes_source_tracker_and_metadata(
        self, mocker, tmp_path
    ):
        """Test pipeline generates source tracking manifest and metadata even when human_reference is None."""
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGTACGTACGT\n")

        config = {
            "tools": {"pbsim3": "pbsim", "ccs": "ccs", "samtools": "samtools"},
            "pacbio_params": {
                "model_type": "qshmm",
                "model_file": str(tmp_path / "model.model"),
                "coverage": 30,
                "pass_num": 3,
                "min_passes": 3,
                "min_rq": 0.99,
            },
        }
        (tmp_path / "model.model").write_text("model data")
        clr_bam = tmp_path / "input_clr_0001.bam"
        hifi_bam = tmp_path / "input_hifi.bam"
        hifi_fastq = tmp_path / "input_hifi.fastq"
        hifi_fastq.write_text("")

        maf_file = tmp_path / "input_clr_0001.maf"
        maf_file.write_text("maf data")

        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.run_pbsim3_simulation",
            return_value=[str(clr_bam)],
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.run_ccs_consensus",
            return_value=str(hifi_bam),
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.convert_bam_to_fastq",
            return_value=str(hifi_fastq),
        )
        mocker.patch(
            "muc_one_up.read_simulator.parsers.pacbio_parser.parse_pacbio_reads",
            return_value=[],
        )
        mock_meta = mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.create_pipeline_metadata"
        )

        mock_tracker = mocker.Mock()
        mock_tracker.annotate_reads.return_value = []

        result = simulate_pacbio_hifi_reads(
            config=config,
            input_fa=str(input_fa),
            human_reference=None,
            source_tracker=mock_tracker,
        )

        assert result == str(hifi_fastq)
        assert mock_tracker.write_manifest.called
        assert not maf_file.exists()
        assert mock_meta.called
        assert mock_meta.call_args[1]["tools_used"] == ["pbsim3", "ccs", "samtools"]

    def test_validates_missing_parameters(self, tmp_path):
        """Test that missing required parameters raise RuntimeError."""
        # Arrange
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGTACGTACGT\n")

        # Missing model_type
        config = {
            "tools": {"pbsim3": "pbsim", "ccs": "ccs", "samtools": "samtools"},
            "pacbio_params": {
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
            "tools": {"pbsim3": "pbsim", "ccs": "ccs", "samtools": "samtools"},
            "pacbio_params": {
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

    def test_validates_threads_parameter_via_ccs(self, tmp_path):
        """Test that invalid threads (e.g. 0) is rejected by early CCS validation."""
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGTACGTACGT\n")

        config = {
            "tools": {"pbsim3": "pbsim", "ccs": "ccs", "samtools": "samtools"},
            "pacbio_params": {
                "model_type": "qshmm",
                "model_file": str(tmp_path / "model.model"),
                "coverage": 30,
                "pass_num": 3,
                "min_passes": 3,
                "min_rq": 0.99,
                "threads": 0,  # Invalid threads! Must be >= 1
            },
        }

        (tmp_path / "model.model").write_text("model data")

        with pytest.raises(RuntimeError, match="Invalid PacBio parameters: Invalid threads"):
            simulate_pacbio_hifi_reads(config=config, input_fa=str(input_fa))

    def test_uses_default_parameters(self, mocker, tmp_path):
        """Test that default parameters are used when not specified."""
        # Arrange
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGTACGTACGT\n")

        config = {
            "tools": {"pbsim3": "pbsim", "ccs": "ccs", "samtools": "samtools"},
            "pacbio_params": {
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
            return [str(clr_bam)]

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

        mocker.patch("muc_one_up.read_simulator.pacbio_pipeline.cleanup_unless_kept")

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
            "tools": {"pbsim3": "pbsim", "ccs": "ccs", "samtools": "samtools"},
            "pacbio_params": {
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
            return [str(clr_bam)]

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

        mocker.patch("muc_one_up.read_simulator.pacbio_pipeline.cleanup_unless_kept")

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

    def test_keep_intermediate_files_honored_and_sidecars_registered(self, mocker, tmp_path):
        """Test Issue #150: honor keep_intermediate_files and register CCS sidecars."""
        from muc_one_up.read_simulator.pipeline_utils import cleanup_unless_kept

        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">hap1\nACGTACGTACGT\n>hap2\nACGTACGTACGT\n")

        config = {
            "tools": {"pbsim3": "pbsim", "ccs": "ccs", "samtools": "samtools"},
            "read_simulation": {"keep_intermediate_files": True},
            "pacbio_params": {
                "model_type": "qshmm",
                "model_file": str(tmp_path / "model.model"),
                "coverage": 30,
                "pass_num": 3,
                "min_passes": 3,
                "min_rq": 0.99,
            },
        }
        (tmp_path / "model.model").write_text("model data")

        clr_bam1 = tmp_path / "input_clr_0001.bam"
        clr_bam2 = tmp_path / "input_clr_0002.bam"
        hifi_bam1 = tmp_path / "input_hifi_0001.bam"
        hifi_bam2 = tmp_path / "input_hifi_0002.bam"
        hifi_bam = tmp_path / "input_hifi.bam"
        hifi_fastq = tmp_path / "input_hifi.fastq"
        pbi_sidecar = tmp_path / "input_hifi_0001.bam.pbi"
        maf_sidecar = tmp_path / "input_clr_0001.maf"
        for p in (clr_bam1, clr_bam2, hifi_bam1, hifi_bam2, hifi_bam, pbi_sidecar, maf_sidecar):
            p.write_bytes(b"DATA")

        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.run_pbsim3_simulation",
            return_value=[str(clr_bam1), str(clr_bam2)],
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.run_ccs_consensus",
            side_effect=[str(hifi_bam1), str(hifi_bam2)],
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.merge_bam_files",
            return_value=str(hifi_bam),
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.convert_bam_to_fastq",
            return_value=str(hifi_fastq),
        )
        mock_cleanup = mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.cleanup_unless_kept",
            side_effect=cleanup_unless_kept,
        )

        simulate_pacbio_hifi_reads(config=config, input_fa=str(input_fa), human_reference=None)

        assert mock_cleanup.called
        assert clr_bam1.exists()
        assert pbi_sidecar.exists()
        assert maf_sidecar.exists()
        # Check that registered intermediate files contain sidecars for both individual HiFi BAMs and merged BAM
        registered_files = mock_cleanup.call_args[0][1]
        assert str(hifi_bam1) in registered_files
        assert f"{hifi_bam1}.pbi" in registered_files
        assert f"{hifi_bam1}.gzi" in registered_files
        assert f"{str(hifi_bam1).removesuffix('.bam')}.ccs_report.txt" in registered_files
        assert f"{str(hifi_bam1).removesuffix('.bam')}.zmw_metrics.json.gz" in registered_files

        assert str(hifi_bam2) in registered_files
        assert f"{hifi_bam2}.pbi" in registered_files
        assert f"{hifi_bam2}.gzi" in registered_files

        assert str(hifi_bam) in registered_files
        assert f"{hifi_bam}.pbi" in registered_files

        # Verify MAF files are registered in intermediate files
        assert f"{str(clr_bam1).removesuffix('.bam')}.maf" in registered_files
        assert f"{str(clr_bam1).removesuffix('.bam')}.maf.gz" in registered_files
        assert f"{str(clr_bam1).removesuffix('.bam')}.ref" in registered_files

    def test_intermediate_files_cleaned_up_on_failure(self, mocker, tmp_path):
        """Test Issue #150: partial files are registered before tools run and cleaned on failure."""
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGTACGTACGT\n")

        config = {
            "tools": {"pbsim3": "pbsim", "ccs": "ccs", "samtools": "samtools"},
            "pacbio_params": {
                "model_type": "qshmm",
                "model_file": str(tmp_path / "model.model"),
                "coverage": 30,
                "pass_num": 3,
                "min_passes": 3,
                "min_rq": 0.99,
            },
        }
        (tmp_path / "model.model").write_text("model data")
        partial_sam = tmp_path / "input_clr_0001.sam"
        partial_maf = tmp_path / "input_clr_0001.maf"

        def fail_pbsim3(**kwargs):
            partial_sam.write_text("@HD\tVN:1.6\n")
            partial_maf.write_text("a\n")
            raise RuntimeError("pbsim3 crashed partway")

        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.run_pbsim3_simulation",
            side_effect=fail_pbsim3,
        )

        with pytest.raises(RuntimeError, match="pbsim3 crashed partway"):
            simulate_pacbio_hifi_reads(config=config, input_fa=str(input_fa), human_reference=None)

        assert not partial_sam.exists()
        assert not partial_maf.exists()

    def test_propagates_external_tool_error_directly(self, mocker, tmp_path):
        """Test Issue #150: ExternalToolError propagates directly without RuntimeError wrapping."""
        from muc_one_up.exceptions import ExternalToolError

        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGTACGTACGT\n")

        config = {
            "tools": {"pbsim3": "pbsim", "ccs": "ccs", "samtools": "samtools"},
            "pacbio_params": {
                "model_type": "qshmm",
                "model_file": str(tmp_path / "model.model"),
                "coverage": 30,
                "pass_num": 3,
                "min_passes": 3,
                "min_rq": 0.99,
            },
        }
        (tmp_path / "model.model").write_text("model data")

        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.run_pbsim3_simulation",
            side_effect=ExternalToolError(tool="pbsim3", exit_code=1, stderr="pbsim3 failed"),
        )
        mocker.patch("muc_one_up.read_simulator.pacbio_pipeline.cleanup_unless_kept")

        with pytest.raises(ExternalToolError, match="pbsim3 failed"):
            simulate_pacbio_hifi_reads(config=config, input_fa=str(input_fa), human_reference=None)

    def test_propagates_file_operation_error_directly(self, mocker, tmp_path):
        """Test Issue #150: FileOperationError propagates directly without RuntimeError wrapping."""
        from muc_one_up.exceptions import FileOperationError

        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGTACGTACGT\n")

        config = {
            "tools": {"pbsim3": "pbsim", "ccs": "ccs", "samtools": "samtools"},
            "pacbio_params": {
                "model_type": "qshmm",
                "model_file": str(tmp_path / "model.model"),
                "coverage": 30,
                "pass_num": 3,
                "min_passes": 3,
                "min_rq": 0.99,
            },
        }
        (tmp_path / "model.model").write_text("model data")
        clr_bam = tmp_path / "input_clr.bam"
        hifi_bam = tmp_path / "input_hifi.bam"

        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.run_pbsim3_simulation",
            return_value=[str(clr_bam)],
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.run_ccs_consensus",
            return_value=str(hifi_bam),
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.convert_bam_to_fastq",
            side_effect=FileOperationError("Could not write FASTQ"),
        )
        mocker.patch("muc_one_up.read_simulator.pacbio_pipeline.cleanup_unless_kept")

        with pytest.raises(FileOperationError, match="Could not write FASTQ"):
            simulate_pacbio_hifi_reads(config=config, input_fa=str(input_fa), human_reference=None)

    def test_resolve_pipeline_outputs_receives_read_simulation_config(self, mocker, tmp_path):
        """Test Issue #150: resolve_pipeline_outputs receives config.get('read_simulation', {})."""
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGTACGTACGT\n")

        rs_cfg = {"output_dir": str(tmp_path / "custom_out")}
        config = {
            "tools": {"pbsim3": "pbsim", "ccs": "ccs", "samtools": "samtools"},
            "read_simulation": rs_cfg,
            "pacbio_params": {
                "model_type": "qshmm",
                "model_file": str(tmp_path / "model.model"),
                "coverage": 30,
                "pass_num": 3,
                "min_passes": 3,
                "min_rq": 0.99,
            },
        }
        (tmp_path / "model.model").write_text("model data")

        mock_resolve = mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.resolve_pipeline_outputs",
            return_value=(tmp_path, "sample"),
        )
        clr_bam = tmp_path / "sample_clr.bam"
        hifi_bam = tmp_path / "sample_hifi.bam"
        hifi_fastq = tmp_path / "sample_hifi.fastq"

        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.run_pbsim3_simulation",
            return_value=[str(clr_bam)],
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.run_ccs_consensus",
            return_value=str(hifi_bam),
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.convert_bam_to_fastq",
            return_value=str(hifi_fastq),
        )
        mocker.patch("muc_one_up.read_simulator.pacbio_pipeline.cleanup_unless_kept")

        simulate_pacbio_hifi_reads(config=config, input_fa=str(input_fa), human_reference=None)

        assert mock_resolve.called
        assert mock_resolve.call_args[0][1] == rs_cfg

    def test_pacbio_manifest_reconciles_ccs_reads_and_streams(self, mocker, tmp_path):
        """Test Issue #127 and #150: manifest reconciles CCS names and streams without list materialization."""
        from types import GeneratorType

        from muc_one_up.read_simulator.source_tracking import ReadOrigin

        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGTACGTACGT\n")

        config = {
            "tools": {"pbsim3": "pbsim", "ccs": "ccs", "samtools": "samtools"},
            "pacbio_params": {
                "model_type": "qshmm",
                "model_file": str(tmp_path / "model.model"),
                "coverage": 30,
                "pass_num": 3,
                "min_passes": 3,
                "min_rq": 0.99,
            },
        }
        (tmp_path / "model.model").write_text("model data")

        clr_bam = tmp_path / "input_clr_0001.bam"
        maf_file = tmp_path / "input_clr_0001.maf"
        maf_file.write_text("dummy maf")
        hifi_bam = tmp_path / "input_hifi.bam"
        hifi_fastq = tmp_path / "input_hifi.fastq"
        # FASTQ contains surviving read movie/1/ccs, but NOT movie/2/ccs (which was rejected by CCS)
        hifi_fastq.write_text("@movie/1/ccs\nACGT\n+\nIIII\n")

        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.run_pbsim3_simulation",
            return_value=[str(clr_bam)],
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.run_ccs_consensus",
            return_value=str(hifi_bam),
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.convert_bam_to_fastq",
            return_value=str(hifi_fastq),
        )
        # MAF parser returns raw subreads for ZMW 1 and rejected ZMW 2
        raw_origins = [
            ReadOrigin(read_id="movie/1/0_500", haplotype=1, ref_start=10, ref_end=510, strand="+"),
            ReadOrigin(
                read_id="movie/1/550_1050", haplotype=1, ref_start=10, ref_end=510, strand="-"
            ),
            ReadOrigin(
                read_id="movie/2/0_500", haplotype=1, ref_start=100, ref_end=600, strand="+"
            ),
        ]
        mocker.patch(
            "muc_one_up.read_simulator.parsers.pacbio_parser.parse_pacbio_reads",
            return_value=raw_origins,
        )
        mocker.patch("muc_one_up.read_simulator.pacbio_pipeline.cleanup_unless_kept")

        mock_tracker = mocker.Mock()

        def mock_annotate(origins):
            yield from origins

        mock_tracker.annotate_reads.side_effect = mock_annotate

        simulate_pacbio_hifi_reads(
            config=config,
            input_fa=str(input_fa),
            human_reference=None,
            source_tracker=mock_tracker,
        )

        assert mock_tracker.write_manifest.called
        manifest_call = mock_tracker.write_manifest.call_args
        # First arg passed to write_manifest is generator (streaming), not materialized list
        assert isinstance(manifest_call[0][0], GeneratorType)
        # Consume generator and verify contents
        annotated_items = list(manifest_call[0][0])
        assert len(annotated_items) == 1
        assert annotated_items[0].read_id == "movie/1/ccs"
        assert annotated_items[0].ref_start == 10
        assert annotated_items[0].ref_end == 510
