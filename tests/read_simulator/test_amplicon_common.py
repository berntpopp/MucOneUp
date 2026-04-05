"""Tests for shared amplicon extraction and preparation."""

from pathlib import Path

import pytest

from muc_one_up.read_simulator.amplicon_common import AmpliconPrep


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

    vntr1 = "ACGT" * 400
    seq1 = "N" * 50 + fwd + vntr1 + rev_rc + "N" * 50

    vntr2 = "ACGT" * 800
    seq2 = "N" * 50 + fwd + vntr2 + rev_rc + "N" * 50

    fasta = tmp_path / "diploid.fa"
    fasta.write_text(f">hap1\n{seq1}\n>hap2\n{seq2}\n")
    return fasta


class TestAmpliconPrep:
    def test_dataclass_fields(self):
        """AmpliconPrep has required fields."""
        prep = AmpliconPrep(
            allele_templates=[Path("/a.fa")],
            allele_coverages=[100],
            output_dir=Path("/out"),
            output_base="test",
            intermediate_files=[],
            is_diploid=False,
        )
        assert prep.is_diploid is False
        assert len(prep.allele_templates) == 1


class TestExtractAndPrepare:
    def test_diploid_returns_two_templates(
        self, diploid_fasta_with_primers, tmp_path, muc1_primers
    ):
        """Diploid input produces 2 template FASTAs with PCR bias split."""
        from muc_one_up.read_simulator.amplicon_common import extract_and_prepare_amplicons

        prep = extract_and_prepare_amplicons(
            input_fa=str(diploid_fasta_with_primers),
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
            total_coverage=100,
            work_dir=tmp_path,
            expected_product_range=None,
            pcr_bias_config={},
            seed=42,
        )

        assert prep.is_diploid is True
        assert len(prep.allele_templates) == 2
        assert len(prep.allele_coverages) == 2
        assert sum(prep.allele_coverages) >= 2  # each allele gets at least 1
        for t in prep.allele_templates:
            assert t.exists()
