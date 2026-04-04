"""End-to-end pipeline tests with real bioinformatics tools.

Run with: uv run pytest tests/e2e/ -m e2e -v
Requires: real tools + model files (via environment variables)
"""

import os
from pathlib import Path

import pytest


def _get_env_path(var: str, description: str) -> str:
    """Get path from environment variable or skip test."""
    value = os.environ.get(var)
    if not value or not Path(value).exists():
        pytest.skip(f"{description}: set {var} environment variable")
    return value


def _reverse_complement(seq: str) -> str:
    """Reverse complement a DNA sequence."""
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(complement.get(b, b) for b in reversed(seq))


@pytest.mark.e2e
@pytest.mark.requires_tools("reseq", "bwa", "samtools", "pblat", "faToTwoBit")
class TestIlluminaE2E:
    """Full Illumina pipeline: FASTA -> paired FASTQ -> sorted BAM."""

    def test_illumina_produces_valid_bam(
        self, tmp_path, e2e_config_base, e2e_diploid_fasta, minimal_config
    ):
        from muc_one_up.read_simulator.pipeline import simulate_reads_pipeline

        reseq_model = _get_env_path("MUCONEUP_RESEQ_MODEL", "ReSeq model")
        sample_bam = _get_env_path("MUCONEUP_SAMPLE_BAM", "Sample BAM")

        config = {**e2e_config_base, **minimal_config}
        config["read_simulation"] = {
            "simulator": "illumina",
            "reseq_model": reseq_model,
            "sample_bam": sample_bam,
            "human_reference": str(e2e_diploid_fasta),
            "output_dir": str(tmp_path / "illumina_out"),
            "read_number": 100,
            "fragment_size": 250,
            "fragment_sd": 30,
            "min_fragment": 150,
            "threads": 2,
            "vntr_capture_efficiency": {"enabled": False},
        }

        result = simulate_reads_pipeline(config, str(e2e_diploid_fasta))

        assert Path(result).exists(), "BAM file should exist"
        assert Path(f"{result}.bai").exists(), "BAM index should exist"
        assert isinstance(result, str)

        # Try to validate with pysam if available
        try:
            import pysam

            with pysam.AlignmentFile(result, "rb") as bam:
                read_count = sum(1 for _ in bam)
            assert read_count > 0, "BAM should contain reads"
        except ImportError:
            pass  # pysam not installed, skip BAM validation


@pytest.mark.e2e
@pytest.mark.requires_tools("nanosim", "minimap2", "samtools")
class TestONTE2E:
    """Full ONT pipeline: FASTA -> FASTQ -> sorted BAM."""

    def test_ont_produces_valid_bam(
        self, tmp_path, e2e_config_base, e2e_diploid_fasta, minimal_config
    ):
        from muc_one_up.read_simulator.ont_pipeline import simulate_ont_reads_pipeline

        nanosim_model = _get_env_path("MUCONEUP_NANOSIM_MODEL", "NanoSim model")

        config = {**e2e_config_base, **minimal_config}
        config["nanosim_params"] = {
            "training_data_path": nanosim_model,
            "coverage": 5,
            "num_threads": 2,
            "seed": 42,
            "enable_split_simulation": True,
        }
        config["read_simulation"] = {
            "simulator": "ont",
            "output_dir": str(tmp_path / "ont_out"),
        }

        result = simulate_ont_reads_pipeline(
            config,
            str(e2e_diploid_fasta),
            human_reference=str(e2e_diploid_fasta),
        )

        assert Path(result).exists(), "BAM file should exist"
        assert isinstance(result, str)

        try:
            import pysam

            with pysam.AlignmentFile(result, "rb") as bam:
                read_count = sum(1 for _ in bam)
            assert read_count > 0, "BAM should contain reads"
        except ImportError:
            pass


@pytest.mark.e2e
@pytest.mark.requires_tools("pbsim3", "ccs", "minimap2", "samtools")
class TestPacBioE2E:
    """Full PacBio HiFi pipeline: FASTA -> CLR -> CCS -> BAM."""

    def test_pacbio_produces_valid_bam(
        self, tmp_path, e2e_config_base, e2e_diploid_fasta
    ):
        from muc_one_up.read_simulator.pacbio_pipeline import simulate_pacbio_hifi_reads

        pbsim3_model = _get_env_path("MUCONEUP_PBSIM3_MODEL", "pbsim3 model")

        config = {**e2e_config_base}
        config["pacbio_params"] = {
            "model_type": "qshmm",
            "model_file": pbsim3_model,
            "coverage": 5,
            "pass_num": 3,
            "min_passes": 3,
            "min_rq": 0.99,
            "threads": 2,
            "seed": 42,
        }

        result = simulate_pacbio_hifi_reads(
            config,
            str(e2e_diploid_fasta),
            human_reference=str(e2e_diploid_fasta),
        )

        assert Path(result).exists(), "BAM file should exist"
        assert isinstance(result, str)

        try:
            import pysam

            with pysam.AlignmentFile(result, "rb") as bam:
                read_count = sum(1 for _ in bam)
            assert read_count > 0, "BAM should contain reads"
        except ImportError:
            pass


@pytest.mark.e2e
@pytest.mark.requires_tools("pbsim3", "ccs", "minimap2", "samtools")
class TestAmpliconE2E:
    """Full Amplicon pipeline: FASTA -> template -> CLR -> CCS -> BAM."""

    def test_amplicon_produces_valid_bam(
        self, tmp_path, e2e_config_base, e2e_diploid_fasta, minimal_config
    ):
        from muc_one_up.read_simulator.amplicon_pipeline import (
            simulate_amplicon_reads_pipeline,
        )

        pbsim3_model = _get_env_path("MUCONEUP_PBSIM3_MODEL", "pbsim3 model")

        config = {**e2e_config_base, **minimal_config}
        config["pacbio_params"] = {
            "model_type": "qshmm",
            "model_file": pbsim3_model,
            "threads": 2,
            "seed": 42,
        }
        config["amplicon_params"] = {
            "forward_primer": minimal_config["constants"]["hg38"]["left"][:20],
            "reverse_primer": _reverse_complement(
                minimal_config["constants"]["hg38"]["right"][-20:]
            ),
        }
        config["read_simulation"] = {
            "simulator": "amplicon",
            "coverage": 10,
        }

        result = simulate_amplicon_reads_pipeline(
            config,
            str(e2e_diploid_fasta),
            human_reference=str(e2e_diploid_fasta),
        )

        assert Path(result).exists(), "Output file should exist"
        assert isinstance(result, str)

        try:
            import pysam

            with pysam.AlignmentFile(result, "rb") as bam:
                read_count = sum(1 for _ in bam)
            assert read_count > 0, "BAM should contain reads"
        except ImportError:
            pass
