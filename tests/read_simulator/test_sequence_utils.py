"""Tests for shared primer binding site detection."""

from muc_one_up.read_simulator.utils.sequence_utils import (
    find_primer_binding_sites,
    reverse_complement,
)


class TestFindPrimerBindingSites:
    """Tests for find_primer_binding_sites()."""

    def test_single_match(self):
        template = "AAACCCGGGTTTAAACCC"
        primer = "CCCGGG"
        sites = find_primer_binding_sites(template, primer)
        assert sites == [3]

    def test_no_match(self):
        template = "AAATTTAAATTT"
        primer = "CCCGGG"
        sites = find_primer_binding_sites(template, primer)
        assert sites == []

    def test_multiple_matches(self):
        template = "ACGTACGTACGT"
        primer = "ACGT"
        sites = find_primer_binding_sites(template, primer)
        assert sites == [0, 4, 8]

    def test_overlapping_matches(self):
        template = "AAAA"
        primer = "AAA"
        sites = find_primer_binding_sites(template, primer)
        assert sites == [0, 1]

    def test_case_normalization(self):
        template = "aaaCCCgggTTT"
        primer = "cccGGG"
        sites = find_primer_binding_sites(template, primer)
        assert sites == [3]

    def test_reverse_complement_search(self):
        # Primer is ACGT, RC is ACGT (palindrome)
        template = "TTTACGTAAA"
        primer = "ACGT"
        sites = find_primer_binding_sites(template, primer, search_reverse_complement=True)
        assert sites == [3]

    def test_reverse_complement_non_palindrome(self):
        # Primer is AAACCC, RC is GGGTTT
        template = "TTTGGGTTTAAA"
        primer = "AAACCC"
        sites = find_primer_binding_sites(template, primer, search_reverse_complement=True)
        assert sites == [3]

    def test_empty_template(self):
        sites = find_primer_binding_sites("", "ACGT")
        assert sites == []

    def test_primer_longer_than_template(self):
        sites = find_primer_binding_sites("ACG", "ACGTACGT")
        assert sites == []


class TestReverseComplement:
    """Tests for reverse_complement()."""

    def test_simple(self):
        assert reverse_complement("ACGT") == "ACGT"  # palindrome

    def test_non_palindrome(self):
        assert reverse_complement("AAACCC") == "GGGTTT"

    def test_case_preserved(self):
        assert reverse_complement("AcGt") == "aCgT"

    def test_empty(self):
        assert reverse_complement("") == ""
