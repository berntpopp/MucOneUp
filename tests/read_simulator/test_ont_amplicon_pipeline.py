"""Tests for ONT amplicon simulation pipeline."""

from pathlib import Path
from unittest.mock import patch

import pytest

from muc_one_up.exceptions import ReadSimulationError
from muc_one_up.read_simulator.amplicon_common import AmpliconPrep


@pytest.fixture
def muc1_primers():
    return {
        "forward": "GGAGAAAAGGAGACTTCGGCTACCCAG",
        "reverse": "GCCGTTGTGCACCAGAGTAGAAGCTGA",
    }


@pytest.fixture
def ont_amplicon_config(tmp_path, muc1_primers):
    model_file = tmp_path / "ERRHMM-ONT-HQ.model"
    model_file.write_text("model")
    return {
        "tools": {
            "pbsim3": "pbsim",
            "samtools": "samtools",
            "minimap2": "minimap2",
        },
        "amplicon_params": {
            "forward_primer": muc1_primers["forward"],
            "reverse_primer": muc1_primers["reverse"],
            "pcr_bias": {"preset": "default"},
        },
        "pacbio_params": {
            "model_type": "errhmm",
            "model_file": str(model_file),
            "threads": 4,
        },
        "read_simulation": {
            "simulator": "ont-amplicon",
            "coverage": 50,
        },
    }


def _make_prep(tmp_path, diploid=True):
    """Create an AmpliconPrep with template files on disk."""
    templates = []
    coverages = []
    count = 2 if diploid else 1
    for i in range(1, count + 1):
        t = tmp_path / f"template_hap{i}.fa"
        t.write_text(">copy\nACGT\n")
        templates.append(t)
        coverages.append(25)
    return AmpliconPrep(
        allele_templates=templates,
        allele_coverages=coverages,
        output_dir=tmp_path,
        output_base="amplicon",
        intermediate_files=[],
        is_diploid=diploid,
    )


class TestOntAmpliconPipeline:
    @patch("muc_one_up.read_simulator.ont_amplicon_pipeline.extract_and_prepare_amplicons")
    @patch("muc_one_up.read_simulator.ont_amplicon_pipeline.run_pbsim3_template_simulation")
    @patch("muc_one_up.read_simulator.ont_amplicon_pipeline.convert_bam_to_fastq")
    @patch("muc_one_up.read_simulator.ont_amplicon_pipeline.align_reads_with_minimap2")
    @patch("muc_one_up.read_simulator.ont_amplicon_pipeline.create_pipeline_metadata")
    @patch("muc_one_up.read_simulator.ont_amplicon_pipeline.cleanup_intermediates")
    def test_calls_pbsim3_with_pass_num_1(
        self,
        mock_cleanup,
        mock_metadata,
        mock_align,
        mock_convert,
        mock_pbsim3,
        mock_prep,
        ont_amplicon_config,
        tmp_path,
    ):
        """ONT amplicon must call pbsim3 with pass_num=1."""
        mock_prep.return_value = _make_prep(tmp_path)
        mock_pbsim3.side_effect = lambda **kw: [kw["output_prefix"] + ".bam"]

        def convert_side_effect(**kw):
            Path(kw["output_fastq"]).write_text("@read\nACGT\n+\nIIII\n")
            return kw["output_fastq"]

        mock_convert.side_effect = convert_side_effect
        mock_align.return_value = str(tmp_path / "aligned.bam")

        fasta = tmp_path / "input.fa"
        fasta.write_text(">seq\nACGT\n")

        from muc_one_up.read_simulator.ont_amplicon_pipeline import (
            simulate_ont_amplicon_pipeline,
        )

        simulate_ont_amplicon_pipeline(
            ont_amplicon_config,
            str(fasta),
            human_reference=str(fasta),
        )

        # Verify pass_num=1 in all pbsim3 calls
        for call in mock_pbsim3.call_args_list:
            assert call.kwargs["pass_num"] == 1

    @patch("muc_one_up.read_simulator.ont_amplicon_pipeline.extract_and_prepare_amplicons")
    @patch("muc_one_up.read_simulator.ont_amplicon_pipeline.run_pbsim3_template_simulation")
    @patch("muc_one_up.read_simulator.ont_amplicon_pipeline.convert_bam_to_fastq")
    @patch("muc_one_up.read_simulator.ont_amplicon_pipeline.align_reads_with_minimap2")
    @patch("muc_one_up.read_simulator.ont_amplicon_pipeline.create_pipeline_metadata")
    @patch("muc_one_up.read_simulator.ont_amplicon_pipeline.cleanup_intermediates")
    def test_uses_map_ont_preset(
        self,
        mock_cleanup,
        mock_metadata,
        mock_align,
        mock_convert,
        mock_pbsim3,
        mock_prep,
        ont_amplicon_config,
        tmp_path,
    ):
        """ONT amplicon must align with map-ont preset."""
        mock_prep.return_value = _make_prep(tmp_path)
        mock_pbsim3.side_effect = lambda **kw: [kw["output_prefix"] + ".bam"]

        def convert_side_effect(**kw):
            Path(kw["output_fastq"]).write_text("@read\nACGT\n+\nIIII\n")
            return kw["output_fastq"]

        mock_convert.side_effect = convert_side_effect
        mock_align.return_value = str(tmp_path / "aligned.bam")

        fasta = tmp_path / "input.fa"
        fasta.write_text(">seq\nACGT\n")

        from muc_one_up.read_simulator.ont_amplicon_pipeline import (
            simulate_ont_amplicon_pipeline,
        )

        simulate_ont_amplicon_pipeline(
            ont_amplicon_config,
            str(fasta),
            human_reference=str(fasta),
        )

        mock_align.assert_called_once()
        assert mock_align.call_args.kwargs["preset"] == "map-ont"

    def test_rejects_source_tracker(self, ont_amplicon_config, tmp_path):
        """Source tracking not supported — must raise."""
        from muc_one_up.read_simulator.ont_amplicon_pipeline import (
            simulate_ont_amplicon_pipeline,
        )

        fasta = tmp_path / "test.fa"
        fasta.write_text(">seq\nACGT\n")

        with pytest.raises(ReadSimulationError, match="source tracking"):
            simulate_ont_amplicon_pipeline(
                ont_amplicon_config,
                str(fasta),
                source_tracker="not_none",
            )
