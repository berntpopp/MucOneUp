"""Contract tests for pipeline output naming, metadata paths, and return values.

These tests validate the orchestration logic of all four simulation pipelines
(Illumina, ONT, PacBio, Amplicon) by mocking external tools and verifying:

- Output file naming follows OutputConfig and input-stem conventions
- Metadata helpers are called with the correct output_base
- Return values are strings (not Path objects)
- BAM vs FASTQ return depends on whether human_reference is supplied
"""

from __future__ import annotations

from muc_one_up.read_simulator.output_config import OutputConfig
from muc_one_up.read_simulator.stages import AlignmentResult, FragmentResult

# ---------------------------------------------------------------------------
# Helper: build a minimal Illumina config that passes all validation
# ---------------------------------------------------------------------------

def _illumina_config(tmp_path):
    """Minimal Illumina config with real reference file + index stubs."""
    ref_file = tmp_path / "hg38.fa"
    ref_file.write_text(">chr1\nATCG\n")
    for ext in [".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"]:
        (tmp_path / f"hg38.fa{ext}").touch()

    return {
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
            "vntr_capture_efficiency": {"enabled": False},
        },
        "reference_assembly": "hg38",
        "reference_genomes": {
            "hg38": {
                "fasta_path": str(ref_file),
                "vntr_region": "chr1:155188487-155192239",
            }
        },
    }


def _patch_illumina(mocker, tmp_path, final_bam: str):
    """Patch all Illumina pipeline dependencies."""
    mocker.patch("muc_one_up.read_simulator.pipeline.check_external_tools")
    mocker.patch(
        "muc_one_up.read_simulator.pipeline.capture_tool_versions", return_value={}
    )
    mocker.patch("muc_one_up.read_simulator.pipeline.log_tool_versions")
    mocker.patch("muc_one_up.read_simulator.pipeline.cleanup_intermediates")
    mocker.patch("muc_one_up.read_simulator.pipeline.create_pipeline_metadata")
    mocker.patch("muc_one_up.read_simulator.pipeline.generate_read_manifest")

    mocker.patch(
        "muc_one_up.read_simulator.pipeline.prepare_fragments",
        return_value=FragmentResult(r1_fastq="/r1.fq.gz", r2_fastq="/r2.fq.gz"),
    )
    mocker.patch(
        "muc_one_up.read_simulator.pipeline.align_and_refine",
        return_value=AlignmentResult(final_bam=final_bam),
    )


# ---------------------------------------------------------------------------
# 1. TestIlluminaOutputNaming
# ---------------------------------------------------------------------------


class TestIlluminaOutputNaming:
    """Verify Illumina pipeline output naming and return-value type."""

    def test_output_bam_with_output_config(self, mocker, tmp_path):
        """BAM path reflects OutputConfig(out_dir, out_base)."""
        from muc_one_up.read_simulator.pipeline import simulate_reads_pipeline

        custom_dir = tmp_path / "custom"
        custom_dir.mkdir()
        expected_bam = str(custom_dir / "myreads.bam")

        config = _illumina_config(tmp_path)
        _patch_illumina(mocker, tmp_path, final_bam=expected_bam)

        input_fa = tmp_path / "sample_001.fa"
        input_fa.write_text(">test\nATCG\n")

        oc = OutputConfig(out_dir=custom_dir, out_base="myreads")
        result = simulate_reads_pipeline(config, str(input_fa), output_config=oc)

        assert "myreads" in result
        assert result.endswith(".bam")

    def test_output_bam_without_output_config(self, mocker, tmp_path):
        """BAM path uses input file stem when no OutputConfig provided."""
        from muc_one_up.read_simulator.pipeline import simulate_reads_pipeline

        expected_bam = str(tmp_path / "sample_001.bam")

        config = _illumina_config(tmp_path)
        _patch_illumina(mocker, tmp_path, final_bam=expected_bam)

        input_fa = tmp_path / "sample_001.fa"
        input_fa.write_text(">test\nATCG\n")

        result = simulate_reads_pipeline(config, str(input_fa))

        assert "sample_001" in result

    def test_returns_string(self, mocker, tmp_path):
        """Return value is str, not Path."""
        from muc_one_up.read_simulator.pipeline import simulate_reads_pipeline

        expected_bam = str(tmp_path / "out.bam")

        config = _illumina_config(tmp_path)
        _patch_illumina(mocker, tmp_path, final_bam=expected_bam)

        input_fa = tmp_path / "sample.fa"
        input_fa.write_text(">test\nATCG\n")

        result = simulate_reads_pipeline(config, str(input_fa))

        assert isinstance(result, str)


# ---------------------------------------------------------------------------
# 2. TestONTOutputNaming
# ---------------------------------------------------------------------------


def _ont_config():
    """Minimal ONT config."""
    return {
        "tools": {
            "nanosim": "simulator.py",
            "minimap2": "minimap2",
            "samtools": "samtools",
        },
        "nanosim_params": {
            "training_data_path": "/models/ont_model",
            "coverage": 30,
        },
        "read_simulation": {},
    }


class TestONTOutputNaming:
    """Regression tests for the input_basename bug in the ONT pipeline."""

    def test_output_with_output_config(self, mocker, tmp_path):
        """create_pipeline_metadata is called with output_base from OutputConfig."""
        from muc_one_up.read_simulator.ont_pipeline import simulate_ont_reads_pipeline

        reads_fq = str(tmp_path / "reads.fq")

        mocker.patch(
            "muc_one_up.read_simulator.ont_pipeline.run_nanosim_simulation",
            return_value=reads_fq,
        )
        mocker.patch(
            "muc_one_up.read_simulator.ont_pipeline.align_ont_reads_with_minimap2"
        )
        mocker.patch(
            "muc_one_up.read_simulator.ont_pipeline.is_diploid_reference",
            return_value=False,
        )
        mock_meta = mocker.patch(
            "muc_one_up.read_simulator.ont_pipeline.create_pipeline_metadata"
        )

        input_fa = tmp_path / "original_name.fa"
        input_fa.write_text(">seq1\nACGT\n")

        oc = OutputConfig(out_dir=tmp_path / "out", out_base="custom_base")
        (tmp_path / "out").mkdir()

        simulate_ont_reads_pipeline(_ont_config(), str(input_fa), output_config=oc)

        # The fix: output_base should come from output_config, not input_basename
        assert mock_meta.called
        call_kwargs = mock_meta.call_args
        output_base_arg = call_kwargs[1].get("output_base") or call_kwargs[0][1]
        assert output_base_arg == "custom_base", (
            f"Expected output_base='custom_base', got '{output_base_arg}'. "
            "This is the regression for the input_basename bug."
        )

    def test_output_bam_uses_output_config_base(self, mocker, tmp_path):
        """Returned BAM path uses OutputConfig base, not input stem."""
        from muc_one_up.read_simulator.ont_pipeline import simulate_ont_reads_pipeline

        out_dir = tmp_path / "results"
        out_dir.mkdir()
        reads_fq = str(out_dir / "reads.fq")

        mocker.patch(
            "muc_one_up.read_simulator.ont_pipeline.run_nanosim_simulation",
            return_value=reads_fq,
        )
        mocker.patch(
            "muc_one_up.read_simulator.ont_pipeline.align_ont_reads_with_minimap2"
        )
        mocker.patch(
            "muc_one_up.read_simulator.ont_pipeline.is_diploid_reference",
            return_value=False,
        )
        mocker.patch("muc_one_up.read_simulator.ont_pipeline.create_pipeline_metadata")

        input_fa = tmp_path / "original_name.fa"
        input_fa.write_text(">seq1\nACGT\n")

        oc = OutputConfig(out_dir=out_dir, out_base="myont")
        result = simulate_ont_reads_pipeline(_ont_config(), str(input_fa), output_config=oc)

        assert "myont" in result
        assert "original_name" not in result


# ---------------------------------------------------------------------------
# 3. TestPacBioOutputNaming
# ---------------------------------------------------------------------------


def _pacbio_config(tmp_path):
    """Minimal PacBio config."""
    model_file = tmp_path / "model.model"
    model_file.write_text("model data")
    return {
        "tools": {
            "pbsim3": "pbsim",
            "ccs": "ccs",
            "samtools": "samtools",
            "minimap2": "minimap2",
        },
        "pacbio_params": {
            "model_type": "qshmm",
            "model_file": str(model_file),
            "coverage": 30,
            "pass_num": 3,
            "min_passes": 3,
            "min_rq": 0.99,
            "threads": 4,
            "seed": 42,
        },
    }


def _patch_pacbio(mocker, tmp_path):
    """Patch all PacBio pipeline dependencies."""
    mocker.patch("muc_one_up.read_simulator.pacbio_pipeline.validate_pbsim3_parameters")
    mocker.patch("muc_one_up.read_simulator.pacbio_pipeline.validate_ccs_parameters")
    mocker.patch(
        "muc_one_up.read_simulator.pacbio_pipeline.run_pbsim3_simulation",
        return_value=[str(tmp_path / "clr.bam")],
    )
    mocker.patch(
        "muc_one_up.read_simulator.pacbio_pipeline.run_ccs_consensus",
        return_value=str(tmp_path / "hifi.bam"),
    )
    mocker.patch(
        "muc_one_up.read_simulator.pacbio_pipeline.convert_bam_to_fastq",
        return_value=str(tmp_path / "hifi.fq"),
    )
    mocker.patch(
        "muc_one_up.read_simulator.pacbio_pipeline.align_reads_with_minimap2",
        return_value=str(tmp_path / "aligned.bam"),
    )
    mocker.patch("muc_one_up.read_simulator.pacbio_pipeline.cleanup_intermediates")
    mocker.patch("muc_one_up.read_simulator.pacbio_pipeline.create_pipeline_metadata")


class TestPacBioOutputNaming:
    """Verify PacBio pipeline return values based on human_reference."""

    def test_returns_bam_with_reference(self, mocker, tmp_path):
        """Returns .bam when human_reference is provided."""
        from muc_one_up.read_simulator.pacbio_pipeline import simulate_pacbio_hifi_reads

        _patch_pacbio(mocker, tmp_path)

        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">seq1\nACGT\n")
        human_ref = tmp_path / "hg38.fa"
        human_ref.write_text(">chr1\nACGT\n")

        result = simulate_pacbio_hifi_reads(
            _pacbio_config(tmp_path),
            str(input_fa),
            human_reference=str(human_ref),
        )

        assert isinstance(result, str)
        assert result.endswith(".bam")

    def test_returns_fastq_without_reference(self, mocker, tmp_path):
        """Returns .fastq when human_reference is None."""
        from muc_one_up.read_simulator.pacbio_pipeline import simulate_pacbio_hifi_reads

        _patch_pacbio(mocker, tmp_path)

        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">seq1\nACGT\n")

        result = simulate_pacbio_hifi_reads(
            _pacbio_config(tmp_path),
            str(input_fa),
            human_reference=None,
        )

        assert isinstance(result, str)
        assert result.endswith(".fastq") or result.endswith(".fq")

    def test_returns_string_not_path(self, mocker, tmp_path):
        """Return value is str, not Path."""
        from muc_one_up.read_simulator.pacbio_pipeline import simulate_pacbio_hifi_reads

        _patch_pacbio(mocker, tmp_path)

        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">seq1\nACGT\n")

        result = simulate_pacbio_hifi_reads(
            _pacbio_config(tmp_path),
            str(input_fa),
            human_reference=None,
        )

        assert isinstance(result, str)


# ---------------------------------------------------------------------------
# 4. TestAmpliconOutputNaming
# ---------------------------------------------------------------------------


def _amplicon_config(tmp_path):
    """Minimal Amplicon config."""
    model_file = tmp_path / "model.model"
    model_file.write_text("model data")
    return {
        "tools": {
            "pbsim3": "pbsim",
            "ccs": "ccs",
            "samtools": "samtools",
            "minimap2": "minimap2",
        },
        "amplicon_params": {
            "forward_primer": "ACGTACGT",
            "reverse_primer": "TGCATGCA",
        },
        "pacbio_params": {
            "model_type": "qshmm",
            "model_file": str(model_file),
            "pass_num": 3,
            "min_passes": 3,
            "min_rq": 0.99,
            "threads": 4,
        },
        "read_simulation": {
            "coverage": 30,
        },
    }


def _patch_amplicon(mocker, tmp_path):
    """Patch all Amplicon pipeline dependencies."""
    mocker.patch(
        "muc_one_up.read_simulator.amplicon_pipeline.is_diploid_reference",
        return_value=False,
    )

    # Mock AmpliconExtractor as a class returning an object with extract()
    mock_extractor = mocker.MagicMock()
    amp_result = mocker.MagicMock()
    amp_result.length = 500
    amp_result.fasta_path = str(tmp_path / "amp.fa")
    mock_extractor.extract.return_value = amp_result
    mocker.patch(
        "muc_one_up.read_simulator.amplicon_pipeline.AmpliconExtractor",
        return_value=mock_extractor,
    )

    mocker.patch("muc_one_up.read_simulator.amplicon_pipeline.generate_template_fasta")
    mocker.patch(
        "muc_one_up.read_simulator.amplicon_pipeline.run_pbsim3_template_simulation",
        return_value=[str(tmp_path / "clr.bam")],
    )
    mocker.patch(
        "muc_one_up.read_simulator.amplicon_pipeline.run_ccs_consensus",
        return_value=str(tmp_path / "hifi.bam"),
    )
    # shutil.copy is called when there is only one haplotype BAM (haploid path).
    # Patch it so no real file copy is attempted.
    mocker.patch("muc_one_up.read_simulator.amplicon_pipeline.shutil.copy")
    mocker.patch(
        "muc_one_up.read_simulator.amplicon_pipeline.merge_bam_files",
        return_value=str(tmp_path / "merged.bam"),
    )
    mocker.patch(
        "muc_one_up.read_simulator.amplicon_pipeline.convert_bam_to_fastq",
        return_value=str(tmp_path / "hifi.fq"),
    )
    mocker.patch(
        "muc_one_up.read_simulator.amplicon_pipeline.align_reads_with_minimap2",
        return_value=str(tmp_path / "aligned.bam"),
    )
    mocker.patch("muc_one_up.read_simulator.amplicon_pipeline.cleanup_intermediates")
    mocker.patch("muc_one_up.read_simulator.amplicon_pipeline.create_pipeline_metadata")


class TestAmpliconOutputNaming:
    """Verify Amplicon pipeline return values based on human_reference."""

    def test_returns_bam_with_reference(self, mocker, tmp_path):
        """Returns .bam when human_reference is provided."""
        from muc_one_up.read_simulator.amplicon_pipeline import (
            simulate_amplicon_reads_pipeline,
        )

        _patch_amplicon(mocker, tmp_path)

        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">seq1\nACGT\n")
        human_ref = tmp_path / "hg38.fa"
        human_ref.write_text(">chr1\nACGT\n")

        result = simulate_amplicon_reads_pipeline(
            _amplicon_config(tmp_path),
            str(input_fa),
            human_reference=str(human_ref),
        )

        assert isinstance(result, str)
        assert result.endswith(".bam")

    def test_returns_fastq_without_reference(self, mocker, tmp_path):
        """Returns FASTQ when human_reference is None."""
        from muc_one_up.read_simulator.amplicon_pipeline import (
            simulate_amplicon_reads_pipeline,
        )

        _patch_amplicon(mocker, tmp_path)

        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">seq1\nACGT\n")

        result = simulate_amplicon_reads_pipeline(
            _amplicon_config(tmp_path),
            str(input_fa),
            human_reference=None,
        )

        assert isinstance(result, str)
        assert result.endswith(".fastq") or result.endswith(".fq")

    def test_returns_string_not_path(self, mocker, tmp_path):
        """Return value is str, not Path."""
        from muc_one_up.read_simulator.amplicon_pipeline import (
            simulate_amplicon_reads_pipeline,
        )

        _patch_amplicon(mocker, tmp_path)

        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">seq1\nACGT\n")

        result = simulate_amplicon_reads_pipeline(
            _amplicon_config(tmp_path),
            str(input_fa),
            human_reference=None,
        )

        assert isinstance(result, str)
