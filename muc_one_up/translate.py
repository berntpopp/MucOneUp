# muc_one_up/translate.py
"""ORF prediction and translation for simulated haplotypes.

This module provides Open Reading Frame (ORF) detection and translation
functionality for simulated MUC1 haplotype sequences. Uses orfipy_core for
efficient ORF finding and implements custom translation with the standard
genetic code.

Key Functions:
    predict_orfs_in_haplotypes: Find and translate ORFs in all haplotypes
    find_orfs_in_memory: Use orfipy to find ORFs in DNA sequence
    dna_to_protein: Translate DNA sequence to protein using codon table
    reverse_complement: Generate reverse complement of DNA sequence
    write_peptides_to_fasta: Write predicted peptides to FASTA file

Key Constants:
    CODON_TABLE: Standard genetic code mapping codons to amino acids

Workflow:
    1. Detect ORFs using orfipy_core (both strands)
    2. Filter by minimum amino acid length (orf_min_aa)
    3. Translate DNA to protein sequences
    4. Optional: Filter by required N-terminal prefix
    5. Write results to peptide FASTA file

Example:
    Predict ORFs with minimum 100 aa length::

        from muc_one_up.translate import run_orf_finder_in_memory

        results = [(seq1, chain1), (seq2, chain2)]
        run_orf_finder_in_memory(
            results=results,
            output_pep="output.pep.fa",
            orf_min_aa=100,
            required_prefix="MTRV"  # Optional N-terminal filter
        )

Notes:
    - Supports both forward and reverse strand ORF detection
    - Stop codons: TAA, TAG, TGA (standard genetic code)
    - Start codons: ATG, TTG, CTG (common bacterial starts)
    - Uses '_' character for stop codons in CODON_TABLE

See Also:
    - toxic_protein_detector.py: Analyze translated ORFs for toxic features
"""

import logging
from pathlib import Path

import orfipy_core

#: Standard genetic code codon table.
#:
#: Maps DNA codons (3-letter strings) to single-letter amino acid codes.
#: Stop codons (TAA, TAG, TGA) are represented by underscore ('_').
#:
#: Format:
#:     {codon: amino_acid, ...}
#:     Example: {"ATG": "M", "TAA": "_", "GCT": "A"}
#:
#: Properties:
#:     - 64 total codons (4^3 combinations)
#:     - 20 standard amino acids + stop
#:     - Multiple codons per amino acid (degeneracy)
#:     - Case-sensitive (uppercase DNA bases)
#:
#: Usage:
#:     Translate codon to amino acid::
#:
#:         from muc_one_up.translate import CODON_TABLE
#:         amino_acid = CODON_TABLE["ATG"]  # Returns "M" (methionine)
#:         amino_acid = CODON_TABLE["TAA"]  # Returns "_" (stop)
#:
#: See Also:
#:     - dna_to_protein: Translation function using this table
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


def reverse_complement(seq: str) -> str:
    """Generate reverse complement of DNA sequence.

    Handles both uppercase and lowercase DNA bases (A, T, G, C, N).
    Unknown bases are converted to 'N'.

    Args:
        seq: DNA sequence string (any case)

    Returns:
        Reverse complement DNA sequence

    Example:
        >>> reverse_complement("ATCG")
        'CGAT'
        >>> reverse_complement("atcgN")
        'Ncgat'
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


def dna_to_protein(
    dna: str, codon_table: dict[str, str] = CODON_TABLE, include_stop: bool = False
) -> str:
    """Translate DNA sequence to protein using codon table.

    Reads DNA in 3-base codons starting from position 0 (in-frame translation).
    Stops at first stop codon unless include_stop=True. Unknown codons are
    translated to 'X'.

    Args:
        dna: DNA sequence to translate
        codon_table: Mapping of codons to amino acids, defaults to standard code
        include_stop: If True, include stop codon character '_' in output

    Returns:
        Translated protein sequence as string of amino acid codes

    Example:
        >>> dna_to_protein("ATGGCTTAA")
        'MA'  # Stop codon TAA terminates translation
        >>> dna_to_protein("ATGGCTTAA", include_stop=True)
        'MA_'  # Stop codon included
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
    seq: str,
    min_len: int = 30,
    max_len: int = 1000000,
    strand: str = "b",
    starts: list[str] | None = None,
    stops: list[str] | None = None,
    include_stop: bool = False,
    partial5: bool = False,
    partial3: bool = False,
) -> list[tuple[int, int, str, str]]:
    """Find Open Reading Frames in DNA sequence using orfipy.

    Uses orfipy_core for efficient ORF detection on one or both strands.
    Default start codons: ATG, TTG, CTG (standard + alternative bacterial starts).
    Default stop codons: TAA, TAG, TGA (standard genetic code).

    Args:
        seq: DNA sequence to search
        min_len: Minimum ORF length in nucleotides
        max_len: Maximum ORF length in nucleotides
        strand: Strand option ('+' for forward, '-' for reverse, 'b' for both)
        starts: List of start codon sequences, defaults to ["ATG", "TTG", "CTG"]
        stops: List of stop codon sequences, defaults to ["TAA", "TAG", "TGA"]
        include_stop: If True, include stop codon in ORF coordinates
        partial5: Allow partial ORFs at 5' end (no start codon)
        partial3: Allow partial ORFs at 3' end (no stop codon)

    Returns:
        List of ORF tuples: (start_pos, stop_pos, strand, description)
        Positions are 0-based indices into the sequence
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


def predict_orfs_in_haplotypes(
    results: list[tuple[str, list[str]]],
    min_len: int = 30,
    orf_min_aa: int = 100,
    required_prefix: str | None = None,
) -> dict[str, list[tuple[str, str, int, int, str, str]]]:
    """Predict and translate ORFs for all haplotypes.

    Finds ORFs in each haplotype sequence, translates to protein, and filters
    by minimum amino acid length and optional N-terminal prefix requirement.

    Args:
        results: List of (dna_sequence, repeat_chain) tuples
        min_len: Minimum ORF length in nucleotides for orfipy detection
        orf_min_aa: Minimum peptide length in amino acids (post-translation filter)
        required_prefix: Optional N-terminal amino acid prefix filter
                        (e.g., "MTRV" to require specific signal peptide)

    Returns:
        Dictionary mapping haplotype IDs to ORF lists:
        {
            "haplotype_1": [
                (orf_id, peptide, start, stop, strand, description),
                ...
            ],
            ...
        }

    Example:
        >>> results = [(seq1, chain1), (seq2, chain2)]
        >>> orfs = predict_orfs_in_haplotypes(results, orf_min_aa=100)
        >>> orfs["haplotype_1"]
        [('haplotype_1_ORF1', 'MTRV...', 245, 1523, '+', 'len=1278')]
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


def write_peptides_to_fasta(haplotype_orfs: dict[str, list[tuple]], output_pep: str) -> None:
    """Write predicted peptide sequences to FASTA file.

    Args:
        haplotype_orfs: Dictionary mapping haplotype IDs to ORF lists
                       (from predict_orfs_in_haplotypes)
        output_pep: Output FASTA filename for peptide sequences

    Raises:
        Exception: If file writing fails (logged and re-raised)

    Example Output:
        >haplotype_1_ORF1 strand=+ len=1278
        MTRVPGTRPALLLLLVPLLLQTGALA...
        >haplotype_1_ORF2 strand=- len=456
        MKRELLSAAVPV...
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


def run_orf_finder_in_memory(
    results: list[tuple[str, list[str]]],
    output_pep: str,
    min_len: int = 30,
    orf_min_aa: int = 100,
    required_prefix: str | None = None,
) -> None:
    """Predict ORFs and write peptides to FASTA (complete pipeline).

    High-level orchestration function combining ORF detection, translation,
    filtering, and FASTA output. Convenience wrapper for common workflow.

    Args:
        results: List of (dna_sequence, repeat_chain) tuples from simulation
        output_pep: Output FASTA filename for translated peptides
        min_len: Minimum ORF length in nucleotides
        orf_min_aa: Minimum peptide length in amino acids
        required_prefix: Optional N-terminal prefix requirement

    Example:
        >>> from muc_one_up.simulate import simulate_diploid
        >>> from muc_one_up.translate import run_orf_finder_in_memory
        >>>
        >>> results = simulate_diploid(config)
        >>> run_orf_finder_in_memory(results, "output.pep.fa", orf_min_aa=100)
    """
    haplotype_orfs = predict_orfs_in_haplotypes(
        results, min_len=min_len, orf_min_aa=orf_min_aa, required_prefix=required_prefix
    )
    write_peptides_to_fasta(haplotype_orfs, output_pep)
