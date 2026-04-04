"""Tests for primer-based amplicon extraction."""

import pytest

from muc_one_up.exceptions import AmpliconExtractionError
from muc_one_up.read_simulator.utils.amplicon_extractor import (
    AmpliconExtractor,
    AmpliconResult,
)


@pytest.fixture
def muc1_primers():
    """MUC1 VNTR amplification primers (Wenzel et al. 2018)."""
    return {
        "forward": "GGAGAAAAGGAGACTTCGGCTACCCAG",
        "reverse": "GCCGTTGTGCACCAGAGTAGAAGCTGA",
    }


@pytest.fixture
def haplotype_fasta_with_primers(tmp_path, muc1_primers):
    """Create a haplotype FASTA containing primer binding sites."""
    from Bio.Seq import Seq

    fwd = muc1_primers["forward"]
    rev_rc = str(Seq(muc1_primers["reverse"]).reverse_complement())
    vntr_region = "ACGT" * 500  # 2000bp fake VNTR
    sequence = "N" * 100 + fwd + vntr_region + rev_rc + "N" * 100
    fasta = tmp_path / "haplotype.fa"
    fasta.write_text(f">hap1\n{sequence}\n")
    return fasta, len(fwd) + len(vntr_region) + len(rev_rc)


@pytest.fixture
def haplotype_fasta_no_primers(tmp_path):
    """Create a haplotype FASTA without primer binding sites."""
    fasta = tmp_path / "no_primers.fa"
    fasta.write_text(">hap1\n" + "ACGT" * 500 + "\n")
    return fasta


class TestAmpliconExtractor:
    """Tests for AmpliconExtractor."""

    def test_successful_extraction(self, haplotype_fasta_with_primers, muc1_primers, tmp_path):
        extractor = AmpliconExtractor(
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
        )
        fasta, expected_len = haplotype_fasta_with_primers
        output = tmp_path / "amplicon.fa"

        result = extractor.extract(str(fasta), str(output))

        assert isinstance(result, AmpliconResult)
        assert result.length == expected_len
        assert result.fasta_path == str(output)
        assert output.exists()

    def test_forward_primer_not_found(self, haplotype_fasta_no_primers, muc1_primers, tmp_path):
        extractor = AmpliconExtractor(
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
        )
        output = tmp_path / "amplicon.fa"

        with pytest.raises(AmpliconExtractionError, match=r"Forward primer.*not found"):
            extractor.extract(str(haplotype_fasta_no_primers), str(output))

    def test_reverse_primer_not_found(self, tmp_path, muc1_primers):
        """Forward primer present but reverse primer missing."""
        fwd = muc1_primers["forward"]
        fasta = tmp_path / "partial.fa"
        fasta.write_text(f">hap1\nNNNN{fwd}ACGTACGT\n")
        output = tmp_path / "amplicon.fa"

        extractor = AmpliconExtractor(
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
        )

        with pytest.raises(AmpliconExtractionError, match=r"Reverse primer.*not found"):
            extractor.extract(str(fasta), str(output))

    def test_multiple_forward_sites_raises(self, tmp_path, muc1_primers):
        """Multiple forward primer binding sites should raise."""
        from Bio.Seq import Seq

        fwd = muc1_primers["forward"]
        rev_rc = str(Seq(muc1_primers["reverse"]).reverse_complement())
        sequence = fwd + "ACGT" * 100 + fwd + "ACGT" * 100 + rev_rc
        fasta = tmp_path / "multi.fa"
        fasta.write_text(f">hap1\n{sequence}\n")
        output = tmp_path / "amplicon.fa"

        extractor = AmpliconExtractor(
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
        )

        with pytest.raises(AmpliconExtractionError, match="Multiple forward primer"):
            extractor.extract(str(fasta), str(output))

    def test_product_range_validation_pass(
        self, haplotype_fasta_with_primers, muc1_primers, tmp_path
    ):
        extractor = AmpliconExtractor(
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
            expected_product_range=(100, 10000),
        )
        fasta, _ = haplotype_fasta_with_primers
        output = tmp_path / "amplicon.fa"

        result = extractor.extract(str(fasta), str(output))
        assert result.length > 0

    def test_product_range_validation_fail(
        self, haplotype_fasta_with_primers, muc1_primers, tmp_path
    ):
        extractor = AmpliconExtractor(
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
            expected_product_range=(10, 50),
        )
        fasta, _ = haplotype_fasta_with_primers
        output = tmp_path / "amplicon.fa"

        with pytest.raises(AmpliconExtractionError, match="outside expected"):
            extractor.extract(str(fasta), str(output))

    def test_case_insensitive_matching(self, tmp_path, muc1_primers):
        """Primers should match regardless of case in the template."""
        from Bio.Seq import Seq

        fwd = muc1_primers["forward"].lower()
        rev_rc = str(Seq(muc1_primers["reverse"]).reverse_complement()).lower()
        sequence = "nnnn" + fwd + "acgt" * 100 + rev_rc + "nnnn"
        fasta = tmp_path / "lowercase.fa"
        fasta.write_text(f">hap1\n{sequence}\n")
        output = tmp_path / "amplicon.fa"

        extractor = AmpliconExtractor(
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
        )

        result = extractor.extract(str(fasta), str(output))
        assert result.length > 0
