"""Tests for reference assembly management integration in read simulation pipelines."""

import pytest

from muc_one_up.exceptions import ConfigurationError


class TestIlluminaPipelineReferenceIntegration:
    """Test reference assembly management in Illumina pipeline."""

    def test_uses_reference_from_reference_genomes(self, tmp_path, mocker):
        """Test pipeline uses reference_genomes section when available."""
        ref_file = tmp_path / "hg38.fa"
        ref_file.write_text(">chr1\nATCG\n")

        # Create mock reference indices
        for ext in [".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"]:
            (tmp_path / f"hg38.fa{ext}").touch()

        config = {
            "tools": {
                "reseq": "reseq",
                "faToTwoBit": "faToTwoBit",
                "samtools": "samtools",
                "pblat": "pblat",
                "bwa": "bwa",
            },
            "read_simulation": {
                "reseq_model": str(tmp_path / "model"),
                "sample_bam": str(tmp_path / "sample.bam"),
                "output_dir": str(tmp_path),
            },
            "reference_assembly": "hg38",
            "reference_genomes": {
                "hg38": {
                    "fasta_path": str(ref_file),
                    "vntr_region": "chr1:155188487-155192239",
                }
            },
        }

        # Mock all the pipeline steps
        mocker.patch("muc_one_up.read_simulator.pipeline.check_external_tools")
        mocker.patch("muc_one_up.read_simulator.pipeline.replace_Ns")
        mocker.patch("muc_one_up.read_simulator.pipeline.generate_systematic_errors")
        mocker.patch("muc_one_up.read_simulator.pipeline.fa_to_twobit")
        mocker.patch("muc_one_up.read_simulator.pipeline.extract_subset_reference")
        mocker.patch("muc_one_up.read_simulator.pipeline.run_pblat")
        mocker.patch("muc_one_up.read_simulator.pipeline.simulate_fragments")
        mocker.patch("muc_one_up.read_simulator.pipeline.create_reads")
        mocker.patch("muc_one_up.read_simulator.pipeline.split_reads")

        # Mock align_reads to capture the reference used
        mock_align = mocker.patch("muc_one_up.read_simulator.pipeline.align_reads")

        # Mock cleanup
        mocker.patch("muc_one_up.read_simulator.pipeline.cleanup_files")

        # Import after mocking to avoid import-time side effects
        from muc_one_up.read_simulator.pipeline import simulate_reads_pipeline

        input_fa = tmp_path / "test.fa"
        input_fa.write_text(">test\nATCG\n")

        simulate_reads_pipeline(config, str(input_fa))

        # Verify align_reads was called with reference from reference_genomes
        assert mock_align.called
        call_args = mock_align.call_args
        assert call_args[0][2] == str(ref_file)  # human_reference argument

    def test_fallback_to_old_human_reference(self, tmp_path, mocker):
        """Test pipeline falls back to old human_reference config."""
        ref_file = tmp_path / "old_ref.fa"
        ref_file.write_text(">chr1\nATCG\n")

        config = {
            "tools": {
                "reseq": "reseq",
                "faToTwoBit": "faToTwoBit",
                "samtools": "samtools",
                "pblat": "pblat",
                "bwa": "bwa",
            },
            "read_simulation": {
                "reseq_model": str(tmp_path / "model"),
                "sample_bam": str(tmp_path / "sample.bam"),
                "human_reference": str(ref_file),  # Old config format
                "output_dir": str(tmp_path),
            },
        }

        # Mock all pipeline steps
        mocker.patch("muc_one_up.read_simulator.pipeline.check_external_tools")
        mocker.patch("muc_one_up.read_simulator.pipeline.replace_Ns")
        mocker.patch("muc_one_up.read_simulator.pipeline.generate_systematic_errors")
        mocker.patch("muc_one_up.read_simulator.pipeline.fa_to_twobit")
        mocker.patch("muc_one_up.read_simulator.pipeline.extract_subset_reference")
        mocker.patch("muc_one_up.read_simulator.pipeline.run_pblat")
        mocker.patch("muc_one_up.read_simulator.pipeline.simulate_fragments")
        mocker.patch("muc_one_up.read_simulator.pipeline.create_reads")
        mocker.patch("muc_one_up.read_simulator.pipeline.split_reads")
        mock_align = mocker.patch("muc_one_up.read_simulator.pipeline.align_reads")
        mocker.patch("muc_one_up.read_simulator.pipeline.cleanup_files")

        from muc_one_up.read_simulator.pipeline import simulate_reads_pipeline

        input_fa = tmp_path / "test.fa"
        input_fa.write_text(">test\nATCG\n")

        simulate_reads_pipeline(config, str(input_fa))

        # Verify fallback to old human_reference
        assert mock_align.called
        call_args = mock_align.call_args
        assert call_args[0][2] == str(ref_file)

    def test_raises_error_when_no_reference_configured(self, tmp_path, mocker):
        """Test pipeline raises error when no reference configured."""
        config = {
            "tools": {
                "reseq": "reseq",
                "faToTwoBit": "faToTwoBit",
                "samtools": "samtools",
                "pblat": "pblat",
                "bwa": "bwa",
            },
            "read_simulation": {
                "reseq_model": str(tmp_path / "model"),
                "sample_bam": str(tmp_path / "sample.bam"),
                "output_dir": str(tmp_path),
                # No reference configured
            },
        }

        mocker.patch("muc_one_up.read_simulator.pipeline.check_external_tools")
        mocker.patch("muc_one_up.read_simulator.pipeline.replace_Ns")
        mocker.patch("muc_one_up.read_simulator.pipeline.generate_systematic_errors")
        mocker.patch("muc_one_up.read_simulator.pipeline.fa_to_twobit")
        mocker.patch("muc_one_up.read_simulator.pipeline.extract_subset_reference")
        mocker.patch("muc_one_up.read_simulator.pipeline.run_pblat")
        mocker.patch("muc_one_up.read_simulator.pipeline.simulate_fragments")
        mocker.patch("muc_one_up.read_simulator.pipeline.create_reads")
        mocker.patch("muc_one_up.read_simulator.pipeline.split_reads")

        from muc_one_up.read_simulator.pipeline import simulate_reads_pipeline

        input_fa = tmp_path / "test.fa"
        input_fa.write_text(">test\nATCG\n")

        with pytest.raises(ConfigurationError, match="human_reference not specified"):
            simulate_reads_pipeline(config, str(input_fa))


class TestONTPipelineReferenceIntegration:
    """Test reference assembly management in ONT pipeline."""

    def test_uses_reference_from_reference_genomes(self, tmp_path, mocker):
        """Test ONT pipeline uses reference_genomes section when available."""
        ref_file = tmp_path / "hg38.fa"
        ref_file.write_text(">chr1\nATCG\n")

        # Create minimap2 index
        (tmp_path / "hg38.fa.mmi").touch()

        config = {
            "tools": {
                "nanosim": "simulator.py",
                "minimap2": "minimap2",
                "samtools": "samtools",
            },
            "nanosim_params": {
                "training_data_path": str(tmp_path / "model"),
                "coverage": 30,
            },
            "read_simulation": {
                "output_dir": str(tmp_path),
            },
            "reference_assembly": "hg38",
            "reference_genomes": {
                "hg38": {
                    "fasta_path": str(ref_file),
                    "vntr_region": "chr1:155188487-155192239",
                }
            },
        }

        # Mock pipeline steps
        mock_nanosim = mocker.patch("muc_one_up.read_simulator.ont_pipeline.run_nanosim_simulation")
        mock_nanosim.return_value = str(tmp_path / "reads.fastq")

        mock_align = mocker.patch(
            "muc_one_up.read_simulator.ont_pipeline.align_ont_reads_with_minimap2"
        )

        # Mock diploid check to avoid split-simulation complexity
        mocker.patch(
            "muc_one_up.read_simulator.ont_pipeline.is_diploid_reference", return_value=False
        )

        from muc_one_up.read_simulator.ont_pipeline import simulate_ont_reads_pipeline

        input_fa = tmp_path / "test.fa"
        input_fa.write_text(">test\nATCG\n")

        simulate_ont_reads_pipeline(config, str(input_fa))

        # Verify align was called with reference from reference_genomes
        assert mock_align.called
        call_args = mock_align.call_args
        assert call_args[1]["human_reference"] == str(ref_file)

    def test_fallback_to_input_fa_on_error(self, tmp_path, mocker):
        """Test ONT pipeline falls back to input_fa when reference loading fails."""
        config = {
            "tools": {
                "nanosim": "simulator.py",
                "minimap2": "minimap2",
                "samtools": "samtools",
            },
            "nanosim_params": {
                "training_data_path": str(tmp_path / "model"),
                "coverage": 30,
            },
            "read_simulation": {
                "output_dir": str(tmp_path),
            },
            "reference_assembly": "nonexistent",  # Will fail
        }

        mock_nanosim = mocker.patch("muc_one_up.read_simulator.ont_pipeline.run_nanosim_simulation")
        mock_nanosim.return_value = str(tmp_path / "reads.fastq")

        mock_align = mocker.patch(
            "muc_one_up.read_simulator.ont_pipeline.align_ont_reads_with_minimap2"
        )
        mocker.patch(
            "muc_one_up.read_simulator.ont_pipeline.is_diploid_reference", return_value=False
        )

        from muc_one_up.read_simulator.ont_pipeline import simulate_ont_reads_pipeline

        input_fa = tmp_path / "test.fa"
        input_fa.write_text(">test\nATCG\n")

        simulate_ont_reads_pipeline(config, str(input_fa))

        # Verify fallback to input_fa
        assert mock_align.called
        call_args = mock_align.call_args
        assert call_args[1]["human_reference"] == str(input_fa)

    def test_uses_user_provided_reference(self, tmp_path, mocker):
        """Test ONT pipeline uses user-provided reference when specified."""
        user_ref = tmp_path / "user_ref.fa"
        user_ref.write_text(">chr1\nATCG\n")

        config = {
            "tools": {
                "nanosim": "simulator.py",
                "minimap2": "minimap2",
                "samtools": "samtools",
            },
            "nanosim_params": {
                "training_data_path": str(tmp_path / "model"),
                "coverage": 30,
            },
            "read_simulation": {
                "output_dir": str(tmp_path),
            },
        }

        mock_nanosim = mocker.patch("muc_one_up.read_simulator.ont_pipeline.run_nanosim_simulation")
        mock_nanosim.return_value = str(tmp_path / "reads.fastq")

        mock_align = mocker.patch(
            "muc_one_up.read_simulator.ont_pipeline.align_ont_reads_with_minimap2"
        )
        mocker.patch(
            "muc_one_up.read_simulator.ont_pipeline.is_diploid_reference", return_value=False
        )

        from muc_one_up.read_simulator.ont_pipeline import simulate_ont_reads_pipeline

        input_fa = tmp_path / "test.fa"
        input_fa.write_text(">test\nATCG\n")

        simulate_ont_reads_pipeline(config, str(input_fa), human_reference=str(user_ref))

        # Verify user-provided reference is used
        assert mock_align.called
        call_args = mock_align.call_args
        assert call_args[1]["human_reference"] == str(user_ref)
