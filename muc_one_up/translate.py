# muc_one_up/translate.py

import logging
from pathlib import Path

import orfipy_core

# Standard genetic code table.
CODON_TABLE = {
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "ATG": "M",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AAC": "N",
    "AAT": "N",
    "AAA": "K",
    "AAG": "K",
    "AGC": "S",
    "AGT": "S",
    "AGA": "R",
    "AGG": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CAC": "H",
    "CAT": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TTC": "F",
    "TTT": "F",
    "TTA": "L",
    "TTG": "L",
    "TAC": "Y",
    "TAT": "Y",
    "TAA": "_",
    "TAG": "_",
    "TGC": "C",
    "TGT": "C",
    "TGA": "_",
    "TGG": "W",
}


def reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence.

    :param seq: DNA sequence.
    :return: Reverse complement of the sequence.
    """
    rc_map = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "N": "N",
        "a": "t",
        "t": "a",
        "g": "c",
        "c": "g",
        "n": "n",
    }
    return "".join(rc_map.get(base, "N") for base in reversed(seq))


def dna_to_protein(dna, codon_table=CODON_TABLE, include_stop=False):
    """
    Translate a DNA sequence in-frame to protein using a codon table.

    :param dna: DNA sequence.
    :param codon_table: Dict mapping codons to amino acids.
    :param include_stop: If True, include stop codon in translation.
    :return: Translated protein sequence.
    """
    protein = []
    for i in range(0, len(dna), 3):
        codon = dna[i : i + 3].upper()
        if len(codon) < 3:
            break
        aa = codon_table.get(codon, "X")
        if aa == "_" and not include_stop:
            break
        protein.append(aa)
    return "".join(protein)


def find_orfs_in_memory(
    seq,
    min_len=30,
    max_len=1000000,
    strand="b",
    starts=None,
    stops=None,
    include_stop=False,
    partial5=False,
    partial3=False,
):
    """
    Use orfipy_core.orfs() to find ORFs in a DNA sequence.

    :param seq: DNA sequence.
    :param min_len: Minimum ORF length in nucleotides.
    :param max_len: Maximum ORF length.
    :param strand: Strand option ('+', '-', or 'b' for both).
    :param starts: List of possible start codons.
    :param stops: List of possible stop codons.
    :param include_stop: Whether to include the stop codon.
    :param partial5: Allow partial ORFs at the 5' end.
    :param partial3: Allow partial ORFs at the 3' end.
    :return: List of ORFs as tuples (start, stop, strand, description).
    """
    if stops is None:
        stops = ["TAA", "TAG", "TGA"]
    if starts is None:
        starts = ["ATG", "TTG", "CTG"]
    return list(
        orfipy_core.orfs(
            seq,
            minlen=min_len,
            maxlen=max_len,
            strand=strand,
            starts=starts,
            stops=stops,
            include_stop=include_stop,
            partial3=partial3,
            partial5=partial5,
            between_stops=False,
        )
    )


def predict_orfs_in_haplotypes(results, min_len=30, orf_min_aa=100, required_prefix=None):
    """
    For each haplotype in results, find ORFs, translate them to peptides,
    and filter based on minimal peptide length and optional prefix.

    :param results: List of tuples (dna_seq, repeat_chain).
    :param min_len: Minimal ORF length in nucleotides (for orfipy).
    :param orf_min_aa: Minimal peptide length in amino acids.
    :param required_prefix: Optional prefix filter for peptides.
    :return: Dict mapping haplotype id to a list of ORF tuples.
    """
    haplotype_orfs = {}
    min_nt_length = orf_min_aa * 3

    for i, (dna_seq, _chain) in enumerate(results, start=1):
        hap_id = f"haplotype_{i}"
        orfs = find_orfs_in_memory(
            seq=dna_seq,
            min_len=min_nt_length,
            max_len=1000000,
            strand="b",
            starts=["ATG", "TTG", "CTG"],
            stops=["TAA", "TAG", "TGA"],
            include_stop=False,
        )
        orf_list = []
        orf_count = 0
        for start, stop, strand, desc in orfs:
            orf_count += 1
            if strand == "+":
                dna_sub = dna_seq[start:stop]
            else:
                dna_sub = reverse_complement(dna_seq[start:stop])
            peptide = dna_to_protein(dna_sub, include_stop=False)
            if len(peptide) < orf_min_aa:
                continue
            if required_prefix and not peptide.startswith(required_prefix):
                continue
            orf_id = f"{hap_id}_ORF{orf_count}"
            orf_list.append((orf_id, peptide, start, stop, strand, desc))
        haplotype_orfs[hap_id] = orf_list

    logging.info("Predicted ORFs for %d haplotypes.", len(results))
    return haplotype_orfs


def write_peptides_to_fasta(haplotype_orfs, output_pep):
    """
    Write predicted peptide sequences to a FASTA file.

    :param haplotype_orfs: Dict mapping haplotype_id to list of ORF tuples.
    :param output_pep: Output FASTA filename.
    """
    try:
        with Path(output_pep).open("w") as out_fh:
            for _hap_id, orf_list in haplotype_orfs.items():
                for orf_id, peptide, _start, _stop, strand, desc in orf_list:
                    out_fh.write(f">{orf_id} strand={strand} {desc}\n")
                    out_fh.write(peptide + "\n")
        logging.info("Peptide FASTA written to %s", output_pep)
    except Exception as e:
        logging.error("Error writing peptide FASTA: %s", e)
        raise


def run_orf_finder_in_memory(results, output_pep, min_len=30, orf_min_aa=100, required_prefix=None):
    """
    High-level function that predicts ORFs for each haplotype in memory,
    applies filters, and writes the peptides to a FASTA file.

    :param results: List of tuples (dna_seq, repeat_chain).
    :param output_pep: Output peptide FASTA filename.
    :param min_len: Minimum ORF length in nucleotides.
    :param orf_min_aa: Minimum peptide length in amino acids.
    :param required_prefix: Optional prefix filter for peptides.
    """
    haplotype_orfs = predict_orfs_in_haplotypes(
        results, min_len=min_len, orf_min_aa=orf_min_aa, required_prefix=required_prefix
    )
    write_peptides_to_fasta(haplotype_orfs, output_pep)
