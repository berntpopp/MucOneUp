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

    def test_multiple_reverse_sites_raises(self, tmp_path, muc1_primers):
        """Multiple reverse primer binding sites should raise."""
        from Bio.Seq import Seq

        fwd = muc1_primers["forward"]
        rev_rc = str(Seq(muc1_primers["reverse"]).reverse_complement())
        sequence = fwd + "ACGT" * 100 + rev_rc + "ACGT" * 100 + rev_rc
        fasta = tmp_path / "multi_rev.fa"
        fasta.write_text(f">hap1\n{sequence}\n")
        output = tmp_path / "amplicon.fa"

        extractor = AmpliconExtractor(
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
        )

        with pytest.raises(AmpliconExtractionError, match="Multiple reverse primer"):
            extractor.extract(str(fasta), str(output))

    def test_reversed_primer_orientation_raises(self, tmp_path, muc1_primers):
        """Forward primer downstream of reverse primer should raise."""
        from Bio.Seq import Seq

        # Place reverse primer RC BEFORE forward primer
        fwd = muc1_primers["forward"]
        rev_rc = str(Seq(muc1_primers["reverse"]).reverse_complement())
        sequence = "NNNN" + rev_rc + "ACGT" * 100 + fwd + "NNNN"
        fasta = tmp_path / "reversed.fa"
        fasta.write_text(f">hap1\n{sequence}\n")
        output = tmp_path / "amplicon.fa"

        extractor = AmpliconExtractor(
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
        )

        with pytest.raises(AmpliconExtractionError, match="downstream of reverse"):
            extractor.extract(str(fasta), str(output))

    def test_empty_fasta_raises(self, tmp_path, muc1_primers):
        """Empty FASTA should raise."""
        fasta = tmp_path / "empty.fa"
        fasta.write_text("")
        output = tmp_path / "amplicon.fa"

        extractor = AmpliconExtractor(
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
        )

        with pytest.raises(AmpliconExtractionError, match="No sequences found"):
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


def _make_amplicon_fasta(tmp_path, muc1_primers, vntr_bp):
    """Build a FASTA whose amplicon (primers + VNTR fill) is exactly vntr_bp + primer lengths.

    The total amplicon size equals len(fwd) + vntr_bp + len(rev_rc).
    """
    from Bio.Seq import Seq

    fwd = muc1_primers["forward"]
    rev_rc = str(Seq(muc1_primers["reverse"]).reverse_complement())
    # Fill with repeating bases to reach desired VNTR region size
    fill = "ACGT" * ((vntr_bp // 4) + 1)
    vntr_region = fill[:vntr_bp]
    sequence = "N" * 50 + fwd + vntr_region + rev_rc + "N" * 50
    fasta = tmp_path / f"hap_{vntr_bp}.fa"
    fasta.write_text(f">hap1\n{sequence}\n")
    expected_len = len(fwd) + vntr_bp + len(rev_rc)
    return fasta, expected_len


class TestAmpliconProductRangeBoundaries:
    """Boundary tests for expected_product_range with the default [500, 15000] bounds.

    Verifies that short amplicon products are accepted after lowering the
    minimum from 1500 to 500 (#91). Tests use raw bp sizes (not repeat
    counts) since the extractor operates on sequence length, not VNTR
    repeat structure.
    """

    @staticmethod
    def _primer_overhead(muc1_primers):
        """Compute total primer contribution to amplicon length."""
        from Bio.Seq import Seq

        fwd_len = len(muc1_primers["forward"])
        rev_len = len(str(Seq(muc1_primers["reverse"]).reverse_complement()))
        return fwd_len + rev_len

    def test_1200bp_amplicon_passes_with_500_lower_bound(self, tmp_path, muc1_primers):
        """A ~1200bp amplicon (e.g. short VNTR) must pass with [500, 15000]."""
        overhead = self._primer_overhead(muc1_primers)
        target_product = 1200
        fasta, expected_len = _make_amplicon_fasta(
            tmp_path, muc1_primers, vntr_bp=target_product - overhead
        )
        assert abs(expected_len - target_product) <= 2
        output = tmp_path / "amplicon.fa"

        extractor = AmpliconExtractor(
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
            expected_product_range=(500, 15000),
        )

        result = extractor.extract(str(fasta), str(output))
        assert result.length == expected_len

    def test_400bp_amplicon_rejected_below_500(self, tmp_path, muc1_primers):
        """A ~400bp amplicon must be rejected when minimum is 500."""
        overhead = self._primer_overhead(muc1_primers)
        target_product = 400
        fasta, expected_len = _make_amplicon_fasta(
            tmp_path, muc1_primers, vntr_bp=target_product - overhead
        )
        assert abs(expected_len - target_product) <= 2
        output = tmp_path / "amplicon.fa"

        extractor = AmpliconExtractor(
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
            expected_product_range=(500, 15000),
        )

        with pytest.raises(AmpliconExtractionError, match="outside expected"):
            extractor.extract(str(fasta), str(output))

    def test_3000bp_amplicon_passes(self, tmp_path, muc1_primers):
        """A normal-sized ~3000bp amplicon must pass with [500, 15000]."""
        overhead = self._primer_overhead(muc1_primers)
        target_product = 3000
        fasta, expected_len = _make_amplicon_fasta(
            tmp_path, muc1_primers, vntr_bp=target_product - overhead
        )
        assert abs(expected_len - target_product) <= 2
        output = tmp_path / "amplicon.fa"

        extractor = AmpliconExtractor(
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
            expected_product_range=(500, 15000),
        )

        result = extractor.extract(str(fasta), str(output))
        assert result.length == expected_len

    def test_16000bp_amplicon_rejected_above_15000(self, tmp_path, muc1_primers):
        """A ~16000bp amplicon must be rejected when maximum is 15000."""
        overhead = self._primer_overhead(muc1_primers)
        target_product = 16000
        fasta, expected_len = _make_amplicon_fasta(
            tmp_path, muc1_primers, vntr_bp=target_product - overhead
        )
        assert abs(expected_len - target_product) <= 2
        output = tmp_path / "amplicon.fa"

        extractor = AmpliconExtractor(
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
            expected_product_range=(500, 15000),
        )

        with pytest.raises(AmpliconExtractionError, match="outside expected"):
            extractor.extract(str(fasta), str(output))
