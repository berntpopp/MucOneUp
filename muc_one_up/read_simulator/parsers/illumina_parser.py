"""Illumina fragment sidecar writer and read origin parser."""

from __future__ import annotations

import csv
import gzip
import logging

from muc_one_up.exceptions import ReadSimulationError
from muc_one_up.read_simulator.source_tracking import ReadOrigin

logger = logging.getLogger(__name__)


def write_fragment_origins(fragments: list[dict], output_path: str) -> None:
    """Write fragment origin data to a TSV sidecar file.

    Args:
        fragments: List of dicts with keys: fragment_index, chrom,
            fstart, fend, strand.
        output_path: Path to write the TSV file.
    """
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "fragment_index",
                "chrom",
                "fstart",
                "fend",
                "strand",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(fragments)


def parse_illumina_reads(
    sidecar_path: str,
    fastq_r1_path: str | None,
    seq_name_to_haplotype: dict[str, int],
) -> list[ReadOrigin]:
    """Parse Illumina reads using fragment origin sidecar and FASTQ names.

    reseq names each read pair after its fragment FASTA ID, which is the
    sidecar ``fragment_index``, so reads are matched to fragments by name.

    Args:
        sidecar_path: Path to fragment_origins.tsv.
        fastq_r1_path: Path to R1 FASTQ (gzipped). Used to get read
            names. If None or unreadable, uses fragment_index as read_id.
        seq_name_to_haplotype: Map from sequence name
            (e.g. "haplotype_1") to haplotype number (e.g. 1).

    Returns:
        List of ReadOrigin entries (one per read, or per fragment without FASTQ).

    Raises:
        ReadSimulationError: If a read name matches no sidecar fragment.
    """
    # Read sidecar
    fragments: list[dict] = []
    try:
        with open(sidecar_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                fragments.append(row)
    except FileNotFoundError:
        logger.warning("Fragment origins sidecar not found: %s", sidecar_path)
        return []

    if not fragments:
        return []

    read_names: list[str] = _extract_read_names(fastq_r1_path)
    if read_names:
        by_index = {frag["fragment_index"]: frag for frag in fragments}
        unknown = [n for n in read_names if n.removesuffix("/1") not in by_index]
        if unknown:
            raise ReadSimulationError(
                f"{len(unknown)} of {len(read_names)} reads in {fastq_r1_path} match no "
                f"fragment in {sidecar_path} (first: {unknown[0]!r})"
            )
        pairs = [(n, by_index[n.removesuffix("/1")]) for n in read_names]
    else:
        pairs = [(f"fragment_{frag['fragment_index']}", frag) for frag in fragments]

    origins: list[ReadOrigin] = []
    for read_id, frag in pairs:
        origins.append(
            ReadOrigin(
                read_id=read_id,
                haplotype=seq_name_to_haplotype.get(frag["chrom"], 1),
                ref_start=int(frag["fstart"]),
                ref_end=int(frag["fend"]),
                strand=frag["strand"],
            )
        )

    return origins


def _extract_read_names(fastq_path: str | None) -> list[str]:
    """Extract read names from a FASTQ file (plain or gzipped).

    Args:
        fastq_path: Path to FASTQ file, or None.

    Returns:
        List of read name strings (without @ prefix or /1 suffix).
    """
    if not fastq_path:
        return []

    read_names: list[str] = []
    try:
        opener = gzip.open if fastq_path.endswith(".gz") else open
        with opener(fastq_path, "rt") as f:
            line_num = 0
            for line in f:
                line_num += 1
                if line_num % 4 == 1:  # Header lines
                    name = line.strip().lstrip("@")
                    # Split on whitespace to get the read ID (drop comment fields)
                    name = name.split()[0] if name else name
                    read_names.append(name)
    except (FileNotFoundError, OSError) as e:
        logger.warning("Could not read FASTQ for read names: %s", e)

    return read_names
