# muc_one_up/translate.py

import orfipy_core
import sys

# A simple codon->amino_acid dict for standard table (table=1).
# In practice, you can expand this as needed or adapt from orfipy's internal table.
CODON_TABLE = {
    "ATA":"I", "ATC":"I", "ATT":"I", "ATG":"M",
    "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T",
    "AAC":"N", "AAT":"N", "AAA":"K", "AAG":"K",
    "AGC":"S", "AGT":"S", "AGA":"R", "AGG":"R",
    "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L",
    "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P",
    "CAC":"H", "CAT":"H", "CAA":"Q", "CAG":"Q",
    "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R",
    "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V",
    "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A",
    "GAC":"D", "GAT":"D", "GAA":"E", "GAG":"E",
    "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G",
    "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S",
    "TTC":"F", "TTT":"F", "TTA":"L", "TTG":"L",
    "TAC":"Y", "TAT":"Y", "TAA":"_", "TAG":"_",
    "TGC":"C", "TGT":"C", "TGA":"_", "TGG":"W",
}

def reverse_complement(seq):
    """Return reverse complement of a DNA string."""
    rc_map = {"A":"T","T":"A","G":"C","C":"G","N":"N",
              "a":"t","t":"a","g":"c","c":"g","n":"n"}
    return "".join(rc_map.get(base,"N") for base in reversed(seq))

def dna_to_protein(dna, codon_table=CODON_TABLE, include_stop=False):
    """Translate a DNA sequence in-frame to protein using a codon table."""
    protein = []
    # step 3 bases at a time
    for i in range(0, len(dna), 3):
        codon = dna[i:i+3].upper()
        if len(codon) < 3:
            break
        aa = codon_table.get(codon, 'X')  # unknown codons => 'X'
        # Optionally skip stop codons if we want truncated sequences
        if aa == '_' and not include_stop:
            break
        protein.append(aa)
    return "".join(protein)

def find_orfs_in_memory(seq, min_len=30, max_len=1000000, strand='b',
                        starts=['ATG','TTG','CTG'], stops=['TAA','TAG','TGA'],
                        include_stop=False, partial5=False, partial3=False):
    """
    Use the orfipy_core.orfs(...) function to find ORFs in a single DNA sequence.
    Returns a list of (start, stop, strand, description).
    """
    return list(orfipy_core.orfs(
        seq,
        minlen=min_len,
        maxlen=max_len,
        strand=strand,
        starts=starts,
        stops=stops,
        include_stop=include_stop,
        partial3=partial3,
        partial5=partial5,
        between_stops=False
    ))

def predict_orfs_in_haplotypes(
    results, 
    min_len=30, 
    orf_min_aa=100, 
    required_prefix=None
):
    """
    Given 'results' = list of (sequence, chain) for each haplotype,
    find ORFs for each haplotype and translate them to peptides,
    then apply filters:
      - minimal peptide length >= 'orf_min_aa'
      - optional prefix filter if 'required_prefix' is not None

    :param results: list of (dna_seq, repeat_chain)
    :param min_len: minimal ORF length in nucleotides (passed to orfipy)
    :param orf_min_aa: minimal peptide length in amino acids (default=100)
    :param required_prefix: if not None, only keep peptides starting with this string
    :return: dict haplotype_id -> list of (orf_id, peptide_seq, start, stop, strand, desc)
    """
    haplotype_orfs = {}
    # Convert 'orf_min_aa' to minimal nt length for orfipy
    min_nt_length = orf_min_aa * 3
    
    for i, (dna_seq, chain) in enumerate(results, start=1):
        hap_id = f"haplotype_{i}"
        # find ORFs, using min_nt_length for orfipy
        orfs = find_orfs_in_memory(
            seq=dna_seq,
            min_len=min_nt_length,
            max_len=1000000,
            strand='b',
            starts=['ATG','TTG','CTG'],  # or all possible starts
            stops=['TAA','TAG','TGA'],
            include_stop=False
        )
        # Translate and apply filters
        orf_list = []
        orf_count = 0
        for (start, stop, strand, desc) in orfs:
            orf_count += 1
            # Extract the subsequence
            if strand == '+':
                dna_sub = dna_seq[start:stop]
            else:
                dna_sub = reverse_complement(dna_seq[start:stop])

            peptide = dna_to_protein(dna_sub, include_stop=False)
            # Filter by actual AA length
            if len(peptide) < orf_min_aa:
                continue
            # If required_prefix is given, filter by prefix
            if required_prefix and not peptide.startswith(required_prefix):
                continue

            orf_id = f"{hap_id}_ORF{orf_count}"
            orf_list.append((orf_id, peptide, start, stop, strand, desc))

        haplotype_orfs[hap_id] = orf_list

    return haplotype_orfs

def write_peptides_to_fasta(haplotype_orfs, output_pep):
    """
    Write the predicted peptide sequences to a FASTA file.
    :param haplotype_orfs: dict haplotype_id -> list of (orf_id, peptide_seq, start, stop, strand, desc)
    :param output_pep: path to output FASTA
    """
    with open(output_pep, "w") as out_fh:
        for hap_id, orf_list in haplotype_orfs.items():
            for (orf_id, peptide, start, stop, strand, desc) in orf_list:
                out_fh.write(f">{orf_id} strand={strand} {desc}\n")
                out_fh.write(peptide + "\n")

def run_orf_finder_in_memory(
    results, 
    output_pep, 
    min_len=30, 
    orf_min_aa=100, 
    required_prefix=None
):
    """
    High-level function that:
     1) Predicts ORFs (in-memory) for each haplotype in 'results'
     2) Applies optional min AA length and prefix filters
     3) Writes them to 'output_pep' in FASTA format

    :param results: list of (dna_seq, chain)
    :param output_pep: filename for the predicted peptide FASTA
    :param min_len: initial minimal ORF length in nt for orfipy (overridden by orf_min_aa*3)
    :param orf_min_aa: minimal peptide length in amino acids (default=100)
    :param required_prefix: if not None, only keep peptides that start with this
    """
    haplotype_orfs = predict_orfs_in_haplotypes(
        results,
        min_len=min_len, 
        orf_min_aa=orf_min_aa,
        required_prefix=required_prefix
    )
    write_peptides_to_fasta(haplotype_orfs, output_pep)
