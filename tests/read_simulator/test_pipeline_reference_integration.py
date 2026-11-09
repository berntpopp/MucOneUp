"""The alignment reference is resolved once in simulate_reads and shared by every pipeline."""

import logging

import pytest

from muc_one_up.exceptions import ConfigurationError
from muc_one_up.read_simulation import simulate_reads

PIPELINES = {
    "illumina": "muc_one_up.read_simulator.pipeline.simulate_reads_pipeline",
    "ont": "muc_one_up.read_simulator.ont_pipeline.simulate_ont_reads_pipeline",
    "pacbio": "muc_one_up.read_simulator.pacbio_pipeline.simulate_pacbio_hifi_reads",
    "amplicon": "muc_one_up.read_simulator.amplicon_pipeline.simulate_amplicon_reads_pipeline",
    "ont-amplicon": "muc_one_up.read_simulator.ont_amplicon_pipeline.simulate_ont_amplicon_pipeline",
    "ont-fragments": "muc_one_up.read_simulator.ont_fragment_pipeline.simulate_ont_fragment_pipeline",
}


def _config(tmp_path, simulator, **rs):
    genomes_ref = tmp_path / "hg38.fa"
    genomes_ref.write_text(">chr1\nACGT\n")
    return {
        "reference_assembly": "hg38",
        "reference_genomes": {"hg38": {"fasta_path": str(genomes_ref)}},
        "read_simulation": {"simulator": simulator, **rs},
    }


@pytest.mark.parametrize("simulator", sorted(PIPELINES))
def test_every_pipeline_gets_the_same_reference(simulator, tmp_path, mocker):
    """Same config -> same alignment reference, whatever the platform (#2)."""
    pipeline = mocker.patch(PIPELINES[simulator], return_value="out")
    explicit = tmp_path / "GRCh38.fa"
    explicit.write_text(">chr1\nACGT\n")
    config = _config(tmp_path, simulator, human_reference=str(explicit))

    simulate_reads(config, "in.fa")

    assert pipeline.call_args.kwargs["human_reference"] == str(explicit)


@pytest.mark.parametrize("simulator", sorted(PIPELINES))
def test_missing_explicit_reference_raises_before_simulating(simulator, tmp_path, mocker):
    """A typo in human_reference must fail up front, not after simulating (#176)."""
    pipeline = mocker.patch(PIPELINES[simulator], return_value="out")
    config = _config(tmp_path, simulator, human_reference=str(tmp_path / "missing.fa"))

    with pytest.raises(ConfigurationError, match="not found"):
        simulate_reads(config, "in.fa")
    pipeline.assert_not_called()


@pytest.mark.parametrize("simulator", sorted(PIPELINES))
def test_reference_genomes_used_when_no_explicit_reference(simulator, tmp_path, mocker):
    pipeline = mocker.patch(PIPELINES[simulator], return_value="out")

    simulate_reads(_config(tmp_path, simulator), "in.fa")

    assert pipeline.call_args.kwargs["human_reference"] == str(tmp_path / "hg38.fa")


@pytest.mark.parametrize("simulator", sorted(PIPELINES))
def test_skip_alignment_passes_no_reference(simulator, tmp_path, mocker):
    pipeline = mocker.patch(PIPELINES[simulator], return_value="out")
    config = _config(tmp_path, simulator, human_reference="/x.fa", skip_alignment=True)

    simulate_reads(config, "in.fa")

    assert pipeline.call_args.kwargs["human_reference"] is None


def test_broken_reference_is_a_warning_not_a_silent_switch(tmp_path, mocker, caplog):
    pipeline = mocker.patch(PIPELINES["pacbio"], return_value="out")
    config = _config(tmp_path, "pacbio")
    config["reference_genomes"]["hg38"]["fasta_path"] = str(tmp_path / "typo.fa")

    with caplog.at_level(logging.WARNING):
        simulate_reads(config, "in.fa")

    assert pipeline.call_args.kwargs["human_reference"] is None
    assert any(r.levelno == logging.WARNING and "typo.fa" in r.message for r in caplog.records)


def test_illumina_without_reference_raises(tmp_path, mocker):
    mocker.patch(PIPELINES["illumina"], return_value="out")
    config = {"read_simulation": {"simulator": "illumina"}}

    with pytest.raises(ConfigurationError, match="No human reference"):
        simulate_reads(config, "in.fa")


class TestONTAlignmentReference:
    """ONT aligns to the reference it is given and skips alignment without one (#183)."""

    @pytest.fixture
    def run(self, tmp_path, mocker):
        config = {
            "tools": {"nanosim": "simulator.py", "minimap2": "minimap2", "samtools": "samtools"},
            "nanosim_params": {"training_data_path": str(tmp_path / "model"), "coverage": 30},
            "read_simulation": {"output_dir": str(tmp_path)},
        }
        mocker.patch(
            "muc_one_up.read_simulator.ont_pipeline.run_nanosim_simulation",
            return_value=str(tmp_path / "reads.fastq"),
        )
        mocker.patch(
            "muc_one_up.read_simulator.ont_pipeline.is_diploid_reference", return_value=False
        )
        align = mocker.patch("muc_one_up.read_simulator.ont_pipeline.align_reads_with_minimap2")
        input_fa = tmp_path / "test.fa"
        input_fa.write_text(">test\nATCG\n")

        from muc_one_up.read_simulator.ont_pipeline import simulate_ont_reads_pipeline

        def _run(human_reference):
            out = simulate_ont_reads_pipeline(
                config, str(input_fa), human_reference=human_reference
            )
            return align, out

        return _run

    def test_uses_given_reference(self, run):
        align, out = run("/ref/GRCh38.fa")
        assert align.call_args.kwargs["reference"] == "/ref/GRCh38.fa"
        assert out.endswith(".bam")

    def test_without_reference_skips_alignment_and_returns_fastq(self, run, tmp_path):
        align, out = run(None)
        align.assert_not_called()
        assert out == str(tmp_path / "reads.fastq")
