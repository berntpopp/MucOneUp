"""Illumina source manifest generation stage (pipeline stage 12a).

Extracts stage 12a of the Illumina read simulation pipeline: parsing
Illumina read origins from the fragment sidecar, filtering to reads that
survived alignment, annotating them via the source tracker, and writing
the gzip-compressed manifest file.
"""

from __future__ import annotations

import contextlib
import logging
import re
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

# pysam is an optional runtime dependency.  We import it at module scope so
# tests can patch ``muc_one_up.read_simulator.stages.source_manifest.pysam``
# cleanly.  The ImportError guard keeps the module importable in environments
# where pysam is not installed.
try:
    import pysam
except ImportError:  # pragma: no cover
    pysam = None  # type: ignore[assignment]

# parse_illumina_reads is imported at module scope (not inside the function)
# so that tests can patch it at the module level:
#   mocker.patch("muc_one_up.read_simulator.stages.source_manifest.parse_illumina_reads", ...)
# The relative import is safe here because this module is only ever loaded as
# part of the muc_one_up package.
from ..parsers.illumina_parser import parse_illumina_reads  # noqa: E402

_READ_SUFFIX_RE = re.compile(r"/[12]$")


def _normalize_read_id(read_id: str | None) -> str:
    """Strip a /1 or /2 mate suffix from a read identifier.

    Args:
        read_id: Raw read identifier, or None.

    Returns:
        Read identifier with /1 or /2 suffix removed, or empty string if
        *read_id* is None.
    """
    if read_id is None:
        return ""
    return _READ_SUFFIX_RE.sub("", read_id)


def generate_read_manifest(
    sidecar_path: str | None,
    final_bam: str,
    input_fa: str,
    source_tracker: Any | None,
    output_dir: Path,
    output_base: str,
) -> str | None:
    """Generate the Illumina read-source manifest (pipeline stage 12a).

    Parses Illumina read origins from the fragment sidecar file, optionally
    filters to reads that survived alignment (using pysam), annotates them
    via *source_tracker*, and writes a gzip-compressed TSV manifest.

    Args:
        sidecar_path: Path to the fragment origins TSV sidecar file, or None.
        final_bam: Path to the final aligned BAM file (used to filter
            surviving reads).
        input_fa: Path to the per-haplotype FASTA (used to build
            seq_name → haplotype index mapping from FASTA headers).
        source_tracker: A :class:`ReadSourceTracker` instance, or None.
            When None the function returns None immediately.
        output_dir: Directory in which the manifest file will be written.
        output_base: Base name used for output files.

    Returns:
        Absolute path string to the written manifest file, or None if no
        manifest was produced.
    """
    if source_tracker is None or sidecar_path is None:
        return None

    # ------------------------------------------------------------------
    # Build seq_name → haplotype-index mapping from FASTA headers
    # ------------------------------------------------------------------
    seq_name_to_haplotype: dict[str, int] = {}
    hap_idx = 1
    try:
        with open(input_fa) as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    seq_name_to_haplotype[name] = hap_idx
                    hap_idx += 1
    except OSError as exc:
        logger.warning("Could not read input FASTA %s: %s", input_fa, exc)

    # ------------------------------------------------------------------
    # Determine R1 FASTQ path
    # ------------------------------------------------------------------
    reads_fq1 = str(output_dir / f"{output_base}_R1.fastq.gz")

    # ------------------------------------------------------------------
    # Parse read origins from sidecar
    # ------------------------------------------------------------------
    origins = parse_illumina_reads(
        sidecar_path=sidecar_path,
        fastq_r1_path=reads_fq1,
        seq_name_to_haplotype=seq_name_to_haplotype,
    )

    # ------------------------------------------------------------------
    # Filter to reads that survived alignment (pysam optional)
    # ------------------------------------------------------------------
    if pysam is not None:
        try:
            surviving_read_ids: set[str] = set()
            with pysam.AlignmentFile(final_bam, "rb") as bam:
                for aln in bam:
                    normalized = _normalize_read_id(aln.query_name)
                    if normalized:
                        surviving_read_ids.add(normalized)
            origins = [o for o in origins if _normalize_read_id(o.read_id) in surviving_read_ids]
        except (OSError, AttributeError) as exc:
            logger.warning("pysam filtering skipped (%s); using all origins", exc)

    # ------------------------------------------------------------------
    # Annotate and write manifest
    # ------------------------------------------------------------------
    manifest_path = str(output_dir / f"{output_base}_read_manifest.tsv.gz")
    annotated = source_tracker.annotate_reads(origins)
    source_tracker.write_manifest(annotated, manifest_path)

    # ------------------------------------------------------------------
    # Clean up sidecar
    # ------------------------------------------------------------------
    with contextlib.suppress(OSError):
        Path(sidecar_path).unlink(missing_ok=True)

    return manifest_path
