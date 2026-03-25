"""ONT (NanoSim) read origin parser."""

from __future__ import annotations

import logging
import re

from muc_one_up.read_simulator.source_tracking import ReadOrigin

logger = logging.getLogger(__name__)

# NanoSim read name format:
# {ref_name}_{position}_{aligned|unaligned}_{index}_{F|R}_{head}_{middle}_{tail}
# For haplotype refs: haplotype_1_500_aligned_0_F_10_100_5
# The ref_name can contain underscores, so we parse from the right
_NANOSIM_PATTERN = re.compile(r"^(.+?)_(\d+)_(aligned|unaligned)_(\d+)_([FR])_(\d+)_(\d+)_(\d+)$")


def parse_nanosim_reads(
    fastq_path: str,
    haplotype_map: int | None = None,
) -> list[ReadOrigin]:
    """Parse NanoSim FASTQ read names to extract read origins.

    Args:
        fastq_path: Path to NanoSim output FASTQ file.
        haplotype_map: If provided, override haplotype for all reads
            (used in split-sim).

    Returns:
        List of ReadOrigin entries.
    """
    origins: list[ReadOrigin] = []

    with open(fastq_path) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # +
            f.readline()  # quality

            if not header.startswith("@"):
                continue

            read_id = header.lstrip("@")
            match = _NANOSIM_PATTERN.match(read_id)
            if not match:
                logger.warning(  # type: ignore[unreachable]
                    "Could not parse NanoSim read name: %s", read_id
                )
                continue

            ref_name = match.group(1)
            ref_start = int(match.group(2))
            strand_char = match.group(5)
            strand = "+" if strand_char == "F" else "-"

            if haplotype_map is not None:
                haplotype = haplotype_map
            else:
                hap_match = re.match(r"haplotype_(\d+)", ref_name)
                haplotype = int(hap_match.group(1)) if hap_match else 1

            ref_end = ref_start + len(seq)

            origins.append(
                ReadOrigin(
                    read_id=read_id,
                    haplotype=haplotype,
                    ref_start=ref_start,
                    ref_end=ref_end,
                    strand=strand,
                )
            )

    return origins
