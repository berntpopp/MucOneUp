# tests/test_translate.py

import pytest

from muc_one_up.translate import dna_to_protein, find_orfs_in_memory


@pytest.mark.parametrize(
    "dna,expected",
    [
        ("ATG", "M"),  # single codon
        ("ATGCCC", "MP"),  # two codons
        ("ATGTAA", "M"),  # stop codon => truncated if include_stop=False
        ("atgcccgtt", "MPV"),  # lower-case input
    ],
)
def test_dna_to_protein(dna, expected):
    """Test dna_to_protein with simple codons."""
    prot = dna_to_protein(dna, include_stop=False)
    assert prot == expected


def test_find_orfs_default_start_is_atg_only():
    """Default starts match the CLI (ATG only); TTG/CTG must not open an ORF."""
    assert find_orfs_in_memory("TTGAAATAG", min_len=6) == []
    assert len(find_orfs_in_memory("ATGAAATAG", min_len=6)) == 1
