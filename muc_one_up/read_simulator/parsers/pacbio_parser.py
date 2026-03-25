"""PacBio (pbsim3) MAF parser for read source tracking."""

from __future__ import annotations

import logging

from muc_one_up.read_simulator.source_tracking import ReadOrigin

logger = logging.getLogger(__name__)


def parse_pacbio_reads(
    maf_paths: list[str],
    haplotype_index: int,
    aligned_bam: str | None = None,
) -> list[ReadOrigin]:
    """Parse pbsim3 MAF files to extract read origins.

    Args:
        maf_paths: List of MAF file paths (one per haplotype sequence
            from pbsim3).
        haplotype_index: 1-based haplotype number for all reads from
            these MAF files.
        aligned_bam: Optional aligned BAM path for fallback position
            extraction.

    Returns:
        List of ReadOrigin entries.
    """
    origins: list[ReadOrigin] = []

    for maf_path in maf_paths:
        origins.extend(_parse_single_maf(maf_path, haplotype_index))

    return origins


def _parse_single_maf(maf_path: str, haplotype_index: int) -> list[ReadOrigin]:
    """Parse a single MAF file."""
    origins: list[ReadOrigin] = []

    try:
        with open(maf_path) as f:
            lines = f.readlines()
    except FileNotFoundError:
        logger.warning("MAF file not found: %s", maf_path)
        return origins

    # Parse MAF blocks: each block starts with 'a' line, followed by
    # 's' lines. First 's' line is reference, second is read.
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("a"):
            s_lines: list[str] = []
            i += 1
            while i < len(lines) and len(s_lines) < 2:
                sline = lines[i].strip()
                if sline.startswith("s"):
                    s_lines.append(sline)
                elif sline == "" or sline.startswith("a"):
                    break
                i += 1

            if len(s_lines) == 2:
                origin = _parse_maf_block(s_lines[0], s_lines[1], haplotype_index)
                if origin:
                    origins.append(origin)
        i += 1

    return origins


def _parse_maf_block(ref_line: str, read_line: str, haplotype_index: int) -> ReadOrigin | None:
    """Parse a MAF alignment block's s-lines.

    MAF s-line format: s name start size strand srcSize sequence
    """
    ref_parts = ref_line.split()
    read_parts = read_line.split()

    if len(ref_parts) < 7 or len(read_parts) < 7:
        return None

    try:
        ref_start = int(ref_parts[2])
        ref_size = int(ref_parts[3])
        ref_strand = ref_parts[4]

        read_name = read_parts[1]
    except (ValueError, IndexError):
        return None

    ref_end = ref_start + ref_size
    strand = "+" if ref_strand == "+" else "-"

    return ReadOrigin(
        read_id=read_name,
        haplotype=haplotype_index,
        ref_start=ref_start,
        ref_end=ref_end,
        strand=strand,
    )
