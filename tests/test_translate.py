# tests/test_translate.py

from pathlib import Path

import pytest

from muc_one_up.translate import (
    dna_to_protein,
    predict_orfs_in_haplotypes,
    reverse_complement,
    run_orf_finder_in_memory,
)


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


def test_reverse_complement():
    """Test the reverse_complement function with a simple case."""
    seq = "ATGC"
    rc = reverse_complement(seq)  # => GCAT
    assert rc == "GCAT"


def test_predict_orfs_in_haplotypes_basic(tmp_path):
    """
    Test predict_orfs_in_haplotypes with a minimal example.
    We'll use a short DNA that has an ORF of ~9 bases => 3 AAs.
    """
    # Prepare a minimal 'results' list => (dna_seq, chain)
    # Suppose chain doesn't matter for this test
    # DNA has an ORF from 0..9 => ATGAAA (Met + Lys) => "MK"
    results = [
        ("ATGAAATAG", ["dummy"]),  # single haplotype
    ]

    # orf_min_aa=2 => minimal 2 AAs => should keep "MK"
    # required_prefix=None => no prefix filter
    haplotype_orfs = predict_orfs_in_haplotypes(
        results,
        min_len=0,  # we can set 0 for orfipy's minimal
        orf_min_aa=2,
        required_prefix=None,
    )

    assert len(haplotype_orfs) == 1
    orfs = next(iter(haplotype_orfs.values()))
    assert len(orfs) == 1  # we expect one kept ORF
    orf_id, peptide, start, stop, strand, desc = orfs[0]
    # "MK" => length=2
    assert peptide == "MK"


def test_predict_orfs_in_haplotypes_prefix_filter():
    """
    Test prefix filtering. We'll create two ORFs:
      - 'MTSSV...' => meets prefix
      - 'MRR...' => does not meet prefix
    We'll require prefix='MTSSV'.
    """
    # We can embed two ORFs in a single DNA:
    # 1) 'ATGACCAGTTCCGTGA' => 'MTSS' then a stop => length=4
    # 2) 'ATGAGA...' => 'MR...' => no prefix
    # We'll ensure the second ORF is also valid length to see it's filtered out by prefix
    dna_seq = "ATGACCAGTTCCGTGAATGAGAAAATAG"  # or something small
    # The first ORF => 'MTSS' => 4 AAs => let's ensure it's ~4 AAs, no prefix 'MTSSV' though?
    # We'll pretend it matches 'MTSSV' if we only check first 4 chars?
    # Let's store it carefully to match exactly "MTSSV" if we want.
    # Actually let's do 'ATGACCAGCAGTG...' => MTSQ => doesn't matter. We'll keep logic simple.

    # For clarity, let's just check the prefix is 'MTSSV' => only the first 5 AAs must match exactly
    # We'll artificially create that 5AA sequence: 'MTSSV' => codons => ATG ACC TCC TCT GTG => "MTSSV"
    dna_seq = "ATGACCTCCTCTGTGTAGATGAGAATTAA"
    # Breaking it down:
    #   ATG ACC TCC TCT GTG => "MTSSV"
    #   then TAA => stop
    # We'll embed a second ORF "MZZ" that doesn't start with 'MTSSV'
    #   ATG GGC => "MG"? etc. Enough to pass length if orf_min_aa=4

    results = [
        (dna_seq, ["whatever"]),
    ]

    # We'll set orf_min_aa=4 => need 4 AAs => "MTSS" wouldn't be enough => must be 5
    # but we have "MTSSV" => 5 AAs => passes length. The second ORF might also pass length but won't start with "MTSSV"
    haplotype_orfs = predict_orfs_in_haplotypes(
        results, min_len=0, orf_min_aa=5, required_prefix="MTSSV"
    )
    # Expect only 1 ORF surviving => the "MTSSV" one
    orfs_list = next(iter(haplotype_orfs.values()))
    assert len(orfs_list) == 1
    assert orfs_list[0][1].startswith("MTSSV"), "peptide must start with 'MTSSV'"


def test_run_orf_finder_in_memory(tmp_path):
    """
    End-to-end test of run_orf_finder_in_memory with a tmp output file.
    We'll create a short DNA with a single valid ORF, then check the output FASTA.
    """
    dna_seq = "ATGCCCCTGTAA"  # => 'MP' + ??? => we'll see
    results = [
        (dna_seq, ["mock"]),
    ]
    # We'll require min 2 AAs => 6 nt => 'MP' => let's do orf_min_aa=2
    out_file = tmp_path / "orf_pep.fa"

    run_orf_finder_in_memory(
        results, output_pep=str(out_file), min_len=0, orf_min_aa=2, required_prefix=None
    )

    assert out_file.exists()
    with Path(out_file).open() as fh:
        lines = fh.read().splitlines()

    # We expect something like:
    # >haplotype_1_ORF1 strand=+ ID=Seq_ORF...
    # MP...
    # Let's confirm we got at least 2 lines
    assert len(lines) >= 2
    assert lines[0].startswith(">haplotype_1_ORF")  # or something
    # second line => the peptide
    pep = lines[1]
    assert len(pep) >= 2  # 'MP' is 2 length
