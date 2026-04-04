"""Integration tests for amplicon simulation pipeline."""

from unittest.mock import MagicMock, patch

import pytest

from muc_one_up.read_simulator.amplicon_pipeline import simulate_amplicon_reads_pipeline


@pytest.fixture
def muc1_primers():
    return {
        "forward": "GGAGAAAAGGAGACTTCGGCTACCCAG",
        "reverse": "GCCGTTGTGCACCAGAGTAGAAGCTGA",
    }


@pytest.fixture
def diploid_fasta_with_primers(tmp_path, muc1_primers):
    """Create a diploid FASTA with primer sites in both haplotypes."""
    from Bio.Seq import Seq

    fwd = muc1_primers["forward"]
    rev_rc = str(Seq(muc1_primers["reverse"]).reverse_complement())

    # Haplotype 1: shorter VNTR
    vntr1 = "ACGT" * 400  # 1600bp
    seq1 = "N" * 50 + fwd + vntr1 + rev_rc + "N" * 50

    # Haplotype 2: longer VNTR
    vntr2 = "ACGT" * 800  # 3200bp
    seq2 = "N" * 50 + fwd + vntr2 + rev_rc + "N" * 50

    fasta = tmp_path / "diploid.fa"
    fasta.write_text(f">hap1\n{seq1}\n>hap2\n{seq2}\n")
    return fasta


@pytest.fixture
def haploid_fasta_with_primers(tmp_path, muc1_primers):
    """Create a haploid FASTA with primer sites."""
    from Bio.Seq import Seq

    fwd = muc1_primers["forward"]
    rev_rc = str(Seq(muc1_primers["reverse"]).reverse_complement())
    vntr = "ACGT" * 500
    seq = "N" * 50 + fwd + vntr + rev_rc + "N" * 50

    fasta = tmp_path / "haploid.fa"
    fasta.write_text(f">hap1\n{seq}\n")
    return fasta


@pytest.fixture
def amplicon_config(tmp_path, muc1_primers):
    """Minimal config for amplicon pipeline tests."""
    model_file = tmp_path / "test.model"
    model_file.write_text("model data")
    return {
        "tools": {
            "pbsim3": "pbsim",
            "ccs": "ccs",
            "samtools": "samtools",
            "minimap2": "minimap2",
        },
        "amplicon_params": {
            "forward_primer": muc1_primers["forward"],
            "reverse_primer": muc1_primers["reverse"],
            "pcr_bias": {"preset": "default"},
        },
        "pacbio_params": {
            "model_type": "qshmm",
            "model_file": str(model_file),
            "coverage": 100,
            "pass_num": 3,
            "min_passes": 3,
            "min_rq": 0.99,
            "threads": 4,
        },
        "read_simulation": {
            "coverage": 100,
        },
    }


class TestAmpliconPipelineDiploid:
    """Tests for diploid amplicon pipeline."""

    @patch("muc_one_up.read_simulator.amplicon_pipeline.run_pbsim3_template_simulation")
    @patch("muc_one_up.read_simulator.amplicon_pipeline.run_ccs_consensus")
    @patch("muc_one_up.read_simulator.amplicon_pipeline.merge_bam_files")
    @patch("muc_one_up.read_simulator.amplicon_pipeline.convert_bam_to_fastq")
    @patch("muc_one_up.read_simulator.amplicon_pipeline.align_reads_with_minimap2")
    def test_diploid_calls_pbsim3_twice(
        self,
        mock_align,
        mock_convert,
        mock_merge,
        mock_ccs,
        mock_pbsim3,
        diploid_fasta_with_primers,
        amplicon_config,
        tmp_path,
    ):
        """Diploid input should invoke PBSIM3 template mode once per allele."""
        mock_pbsim3.side_effect = lambda **kw: [kw["output_prefix"] + ".bam"]
        mock_ccs.return_value = str(tmp_path / "hifi.bam")
        mock_merge.return_value = str(tmp_path / "merged.bam")
        mock_convert.return_value = str(tmp_path / "hifi.fastq")
        mock_align.return_value = str(tmp_path / "aligned.bam")

        for name in ["hifi.bam", "merged.bam", "hifi.fastq", "aligned.bam"]:
            (tmp_path / name).write_bytes(b"FAKE")

        simulate_amplicon_reads_pipeline(
            config=amplicon_config,
            input_fa=str(diploid_fasta_with_primers),
            human_reference=str(tmp_path / "ref.fa"),
        )

        assert mock_pbsim3.call_count == 2

    @patch("muc_one_up.read_simulator.amplicon_pipeline.run_pbsim3_template_simulation")
    @patch("muc_one_up.read_simulator.amplicon_pipeline.run_ccs_consensus")
    @patch("muc_one_up.read_simulator.amplicon_pipeline.merge_bam_files")
    @patch("muc_one_up.read_simulator.amplicon_pipeline.convert_bam_to_fastq")
    @patch("muc_one_up.read_simulator.amplicon_pipeline.align_reads_with_minimap2")
    def test_pcr_bias_affects_template_counts(
        self,
        mock_align,
        mock_convert,
        mock_merge,
        mock_ccs,
        mock_pbsim3,
        diploid_fasta_with_primers,
        amplicon_config,
        tmp_path,
    ):
        """PCR bias should give shorter allele more template copies."""
        template_fastas = []

        def capture_template(**kwargs):
            template_fastas.append(kwargs["template_fasta"])
            return [kwargs["output_prefix"] + ".bam"]

        mock_pbsim3.side_effect = capture_template
        mock_ccs.return_value = str(tmp_path / "hifi.bam")
        mock_merge.return_value = str(tmp_path / "merged.bam")
        mock_convert.return_value = str(tmp_path / "hifi.fastq")
        mock_align.return_value = str(tmp_path / "aligned.bam")

        for name in ["hifi.bam", "merged.bam", "hifi.fastq", "aligned.bam"]:
            (tmp_path / name).write_bytes(b"FAKE")

        simulate_amplicon_reads_pipeline(
            config=amplicon_config,
            input_fa=str(diploid_fasta_with_primers),
            human_reference=str(tmp_path / "ref.fa"),
        )

        from Bio import SeqIO

        counts = [len(list(SeqIO.parse(f, "fasta"))) for f in template_fastas]
        assert counts[0] > counts[1]
        assert sum(counts) == 100


class TestAmpliconPipelineHaploid:
    """Tests for haploid amplicon pipeline."""

    @patch("muc_one_up.read_simulator.amplicon_pipeline.run_pbsim3_template_simulation")
    @patch("muc_one_up.read_simulator.amplicon_pipeline.run_ccs_consensus")
    @patch("muc_one_up.read_simulator.amplicon_pipeline.convert_bam_to_fastq")
    @patch("muc_one_up.read_simulator.amplicon_pipeline.align_reads_with_minimap2")
    def test_haploid_calls_pbsim3_once(
        self,
        mock_align,
        mock_convert,
        mock_ccs,
        mock_pbsim3,
        haploid_fasta_with_primers,
        amplicon_config,
        tmp_path,
    ):
        """Haploid input should invoke PBSIM3 once with full coverage."""
        mock_pbsim3.side_effect = lambda **kw: [kw["output_prefix"] + ".bam"]
        mock_ccs.return_value = str(tmp_path / "hifi.bam")
        mock_convert.return_value = str(tmp_path / "hifi.fastq")
        mock_align.return_value = str(tmp_path / "aligned.bam")

        for name in ["hifi.bam", "hifi.fastq", "aligned.bam"]:
            (tmp_path / name).write_bytes(b"FAKE")

        simulate_amplicon_reads_pipeline(
            config=amplicon_config,
            input_fa=str(haploid_fasta_with_primers),
            human_reference=str(tmp_path / "ref.fa"),
        )

        assert mock_pbsim3.call_count == 1


class TestAmpliconPipelineTrackingRejection:
    """Tests for read-source tracking rejection."""

    def test_track_read_source_raises_error(
        self,
        haploid_fasta_with_primers,
        amplicon_config,
    ):
        """--track-read-source should raise for amplicon mode."""
        tracker = MagicMock()

        with pytest.raises(RuntimeError, match="not yet supported"):
            simulate_amplicon_reads_pipeline(
                config=amplicon_config,
                input_fa=str(haploid_fasta_with_primers),
                source_tracker=tracker,
            )
