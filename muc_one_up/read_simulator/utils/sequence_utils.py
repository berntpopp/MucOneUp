"""Shared sequence utilities for primer binding site detection.

Provides primer binding site detection used by both the SNaPshot
validation workflow and the amplicon extraction pipeline.
"""

import re

_COMPLEMENT = str.maketrans(
    "ACGTRYMKSWBDHVNacgtrymkswbdhvn",
    "TGCAYRKMSWVHDBNtgcayrkmswvhdbn",
)
_IUPAC_REGEX = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "[AGR]",
    "Y": "[CTY]",
    "M": "[ACM]",
    "K": "[GTK]",
    "S": "[CGS]",
    "W": "[ATW]",
    "B": "[CGTYSKB]",
    "D": "[AGTRWKD]",
    "H": "[ACTYWMH]",
    "V": "[ACGRMSV]",
    "N": "[ACGTRYMKSWBDHVN]",
}


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence.

    Args:
        seq: DNA sequence (IUPAC characters, case preserved).

    Returns:
        Reverse complement string.
    """
    return seq.translate(_COMPLEMENT)[::-1]


def find_primer_binding_sites(
    template: str,
    primer: str,
    search_reverse_complement: bool = False,
) -> list[int]:
    """Find all binding sites for a primer in a template sequence.

    Both template and primer are normalized to uppercase before matching.
    Searches support IUPAC degenerate bases and allow overlapping matches.

    Args:
        template: Template DNA sequence.
        primer: Primer sequence (5'->3').
        search_reverse_complement: If True, search for the reverse complement
            of the primer instead of the primer itself.

    Returns:
        List of 0-based start positions where the primer binds.
    """
    if not primer or len(primer) > len(template):
        return []

    template_upper = template.upper()
    search_seq = primer.upper()

    if search_reverse_complement:
        search_seq = reverse_complement(search_seq)

    pattern = "".join(_IUPAC_REGEX.get(ch, re.escape(ch)) for ch in search_seq)
    return [m.start() for m in re.finditer(f"(?={pattern})", template_upper)]
