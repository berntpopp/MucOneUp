"""Shared sequence utilities for primer binding site detection.

Provides primer binding site detection used by both the SNaPshot
validation workflow and the amplicon extraction pipeline.
"""

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence.

    Args:
        seq: DNA sequence (ACGT characters, case preserved).

    Returns:
        Reverse complement string.
    """
    return seq.translate(_COMPLEMENT)[::-1]


def find_primer_binding_sites(
    template: str,
    primer: str,
    reverse_complement: bool = False,
) -> list[int]:
    """Find all binding sites for a primer in a template sequence.

    Both template and primer are normalized to uppercase before matching.
    Searches use exact string matching with overlapping allowed.

    Args:
        template: Template DNA sequence.
        primer: Primer sequence (5'->3').
        reverse_complement: If True, search for the reverse complement
            of the primer instead of the primer itself.

    Returns:
        List of 0-based start positions where the primer binds.
    """
    template_upper = template.upper()
    search_seq = primer.upper()

    if reverse_complement:
        search_seq = search_seq.translate(_COMPLEMENT)[::-1]

    sites: list[int] = []
    start = 0
    while True:
        pos = template_upper.find(search_seq, start)
        if pos == -1:
            break
        sites.append(pos)
        start = pos + 1

    return sites
