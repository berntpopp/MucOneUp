"""Illumina fragment sidecar writer and read origin parser."""

from __future__ import annotations

import csv
import gzip
import itertools
import logging
from array import array
from collections.abc import Iterator

from muc_one_up.exceptions import ReadSimulationError
from muc_one_up.read_simulator.source_tracking import ReadOrigin

logger = logging.getLogger(__name__)

FRAGMENT_ORIGIN_FIELDS = ("fragment_index", "chrom", "fstart", "fend", "strand")


def write_fragment_origins(fragments: list[dict], output_path: str) -> None:
    """Write fragment origin data to a TSV sidecar file.

    Args:
        fragments: List of dicts with keys: fragment_index, chrom,
            fstart, fend, strand.
        output_path: Path to write the TSV file.
    """
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FRAGMENT_ORIGIN_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(fragments)


def parse_illumina_reads(
    sidecar_path: str,
    fastq_r1_path: str | None,
    seq_name_to_haplotype: dict[str, int],
) -> Iterator[ReadOrigin]:
    """Stream Illumina read origins from the fragment origin sidecar and FASTQ names.

    reseq names each read pair after its fragment FASTA ID, which is the
    sidecar ``fragment_index``, so reads are matched to fragments by name.

    Args:
        sidecar_path: Path to fragment_origins.tsv.
        fastq_r1_path: Path to R1 FASTQ (gzipped). Used to get read
            names. If None or unreadable, uses fragment_index as read_id.
        seq_name_to_haplotype: Map from sequence name
            (e.g. "haplotype_1") to haplotype number (e.g. 1).

    Yields:
        One ReadOrigin per read (or per fragment without FASTQ), in FASTQ order.

    Raises:
        ReadSimulationError: If a read name matches no sidecar fragment.
    """
    try:
        f = open(sidecar_path)
    except FileNotFoundError:
        logger.warning("Fragment origins sidecar not found: %s", sidecar_path)
        return

    names = _extract_read_names(fastq_r1_path)
    first = next(names, None)
    with f:
        rows = csv.DictReader(f, delimiter="\t")
        if first is None:
            for row in rows:
                yield ReadOrigin(
                    read_id=f"fragment_{row['fragment_index']}",
                    haplotype=seq_name_to_haplotype.get(row["chrom"], 1),
                    ref_start=int(row["fstart"]),
                    ref_end=int(row["fend"]),
                    strand=row["strand"],
                )
            return
        # Reads may come in any order, so index fragments by fragment_index in flat
        # arrays (~21 bytes per fragment) instead of one dict per row.
        haps, starts, ends, strands = array("i"), array("q"), array("q"), bytearray()
        for row in rows:
            i = int(row["fragment_index"])
            if i >= len(starts):
                grow = i + 1 - len(starts)
                haps.extend([0] * grow)
                starts.extend([-1] * grow)  # -1 marks an index with no fragment
                ends.extend([0] * grow)
                strands.extend(b"+" * grow)
            haps[i] = seq_name_to_haplotype.get(row["chrom"], 1)
            starts[i] = int(row["fstart"])
            ends[i] = int(row["fend"])
            strands[i] = ord(row["strand"])

    for name in itertools.chain([first], names):
        index = name.removesuffix("/1")
        i = int(index) if index.isdecimal() else -1
        if not 0 <= i < len(starts) or starts[i] < 0:
            raise ReadSimulationError(
                f"Reads in {fastq_r1_path} match no fragment in {sidecar_path} (first: {name!r})"
            )
        yield ReadOrigin(
            read_id=name,
            haplotype=haps[i],
            ref_start=starts[i],
            ref_end=ends[i],
            strand=chr(strands[i]),
        )


def _extract_read_names(fastq_path: str | None) -> Iterator[str]:
    """Yield read names from a FASTQ file (plain or gzipped).

    Args:
        fastq_path: Path to FASTQ file, or None.

    Yields:
        Read name strings (without @ prefix or comment fields).
    """
    if not fastq_path:
        return

    try:
        opener = gzip.open if fastq_path.endswith(".gz") else open
        with opener(fastq_path, "rt") as f:
            for line_num, line in enumerate(f):
                if line_num % 4 == 0:  # Header lines
                    name = line.strip().lstrip("@")
                    # Split on whitespace to get the read ID (drop comment fields)
                    yield name.split()[0] if name else name
    except (FileNotFoundError, OSError) as e:
        logger.warning("Could not read FASTQ for read names: %s", e)
