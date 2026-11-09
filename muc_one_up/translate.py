# muc_one_up/translate.py
"""ORF prediction and translation for simulated haplotypes.

This module provides Open Reading Frame (ORF) detection and translation
functionality for simulated MUC1 haplotype sequences. Uses orfipy_core for
efficient ORF finding and implements custom translation with the standard
genetic code.

Key Functions:
    find_orfs_in_memory: Use orfipy to find ORFs in DNA sequence
    dna_to_protein: Translate DNA sequence to protein using codon table

Key Constants:
    CODON_TABLE: Standard genetic code mapping codons to amino acids

Notes:
    - Supports both forward and reverse strand ORF detection
    - Stop codons: TAA, TAG, TGA (standard genetic code)
    - Default start codon: ATG
    - Uses '_' character for stop codons in CODON_TABLE
    - The ``analyze orfs`` command (cli/analysis.py) drives the full workflow

See Also:
    - toxic_protein_detector.py: Analyze translated ORFs for toxic features
"""

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
    Default start codon: ATG.
    Default stop codons: TAA, TAG, TGA (standard genetic code).

    Args:
        seq: DNA sequence to search
        min_len: Minimum ORF length in nucleotides
        max_len: Maximum ORF length in nucleotides
        strand: Strand option ('+' for forward, '-' for reverse, 'b' for both)
        starts: List of start codon sequences, defaults to ["ATG"]
        stops: List of stop codon sequences, defaults to ["TAA", "TAG", "TGA"]
        include_stop: If True, include stop codon in ORF coordinates
        partial5: Allow partial ORFs at 5' end (no start codon)
        partial3: Allow partial ORFs at 3' end (no stop codon)

    Returns:
        List of ORF tuples: (start_pos, stop_pos, strand, description)
        Positions are 0-based indices into the sequence
    """
    import orfipy_core

    if stops is None:
        stops = ["TAA", "TAG", "TGA"]
    if starts is None:
        starts = ["ATG"]
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
