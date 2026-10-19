"""PacBio (pbsim3) MAF parser for read source tracking."""

from __future__ import annotations

import gzip
import logging
from collections.abc import Sequence
from pathlib import Path

from muc_one_up.exceptions import FileOperationError
from muc_one_up.read_simulator.source_tracking import ReadOrigin

logger = logging.getLogger(__name__)


def reconcile_ccs_reads(
    origins: Sequence[ReadOrigin],
    fastq_path: str | Path,
) -> list[ReadOrigin]:
    """Reconcile raw MAF read origins against surviving CCS reads in a FASTQ file.

    Raises:
        FileOperationError: If the FASTQ file cannot be read or any surviving CCS
            read cannot be reconciled with MAF origins.
    """
    fastq_file = Path(fastq_path)
    if not fastq_file.exists():
        raise FileOperationError(f"HiFi FASTQ file not found for CCS reconciliation: {fastq_path}")

    exact_map: dict[str, ReadOrigin] = {}
    zmw_map: dict[str, ReadOrigin] = {}

    for origin in origins:
        exact_map[origin.read_id] = origin
        if "/" in origin.read_id:
            zmw_map[origin.read_id.rsplit("/", 1)[0]] = origin
        else:
            zmw_map[origin.read_id] = origin

    surviving_origins: list[ReadOrigin] = []
    opener = gzip.open if str(fastq_file).endswith(".gz") else open

    try:
        with opener(fastq_file, "rt") as f:
            for i, line in enumerate(f):
                if i % 4 != 0:
                    continue
                line = line.strip()
                if not line.startswith("@"):
                    continue
                read_name = line[1:].split()[0]

                matched: ReadOrigin | None = None
                if read_name in exact_map:
                    matched = exact_map[read_name]
                elif "/" in read_name and read_name.rsplit("/", 1)[0] in zmw_map:
                    matched = zmw_map[read_name.rsplit("/", 1)[0]]
                elif read_name in zmw_map:
                    matched = zmw_map[read_name]

                if matched is None:
                    raise FileOperationError(
                        f"CCS read {read_name} could not be reconciled with MAF origins"
                    )

                surviving_origins.append(
                    ReadOrigin(
                        read_id=read_name,
                        haplotype=matched.haplotype,
                        ref_start=matched.ref_start,
                        ref_end=matched.ref_end,
                        strand=matched.strand,
                    )
                )
    except (OSError, gzip.BadGzipFile) as exc:
        raise FileOperationError(
            f"Could not read FASTQ file {fastq_path} for CCS reconciliation: {exc}"
        ) from exc

    return surviving_origins


def parse_pacbio_reads(
    maf_paths: list[str],
    haplotype_index: int,
    aligned_bam: str | None = None,
) -> list[ReadOrigin]:
    """Parse pbsim3 MAF files to extract read origins.

    Automatically handles both plain (.maf) and gzip-compressed (.maf.gz)
    files. When both exist for the same base name, only the uncompressed
    file is parsed to avoid double-counting.

    Args:
        maf_paths: List of MAF file paths (one per haplotype sequence
            from pbsim3). May include both .maf and .maf.gz variants.
        haplotype_index: 1-based haplotype number for all reads from
            these MAF files.
        aligned_bam: Optional aligned BAM path for fallback position
            extraction.

    Returns:
        List of ReadOrigin entries.
    """
    origins: list[ReadOrigin] = []

    dedup_map: dict[str, str] = {}
    for maf_path in maf_paths:
        p = Path(maf_path)
        if not p.exists():
            continue
        base = str(p).removesuffix(".gz")
        existing = dedup_map.get(base)
        if existing is None or (existing.endswith(".gz") and not maf_path.endswith(".gz")):
            dedup_map[base] = maf_path
    deduplicated = list(dedup_map.values())

    for maf_path in deduplicated:
        origins.extend(_parse_single_maf(maf_path, haplotype_index))

    return origins


def _parse_single_maf(maf_path: str, haplotype_index: int) -> list[ReadOrigin]:
    """Parse a single MAF file (plain or gzip-compressed), streaming line-by-line."""
    origins: list[ReadOrigin] = []

    opener = gzip.open if maf_path.endswith(".gz") else open
    try:
        with opener(maf_path, "rt") as f:
            s_lines: list[str] = []
            in_block = False

            for line in f:
                line = line.strip()
                if line.startswith("a"):
                    # Start of a new alignment block
                    in_block = True
                    s_lines = []
                elif in_block and line.startswith("s"):
                    s_lines.append(line)
                    if len(s_lines) == 2:
                        origin = _parse_maf_block(s_lines[0], s_lines[1], haplotype_index)
                        if origin:
                            origins.append(origin)
                        in_block = False
                elif line == "" or line.startswith("#"):
                    if in_block and len(s_lines) < 2:
                        in_block = False
                else:
                    # Non-s line inside block resets
                    if in_block and not line.startswith("s"):
                        in_block = False
    except FileNotFoundError:
        logger.warning("MAF file not found: %s", maf_path)
    except (OSError, gzip.BadGzipFile) as exc:
        logger.warning("Could not read MAF file %s: %s", maf_path, exc)

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

        read_name = read_parts[1]
        # In MAF, reference strand is typically always '+'.
        # The read s-line strand carries the actual read orientation.
        read_strand = read_parts[4]
    except (ValueError, IndexError):
        return None

    ref_end = ref_start + ref_size
    strand = "+" if read_strand == "+" else "-"

    return ReadOrigin(
        read_id=read_name,
        haplotype=haplotype_index,
        ref_start=ref_start,
        ref_end=ref_end,
        strand=strand,
    )
