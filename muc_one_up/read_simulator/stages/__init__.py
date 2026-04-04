"""Result dataclasses for Illumina pipeline stages.

These frozen dataclasses define the contracts between the orchestrator
(pipeline.py) and the extracted stage modules.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class FragmentResult:
    """Result of the fragment preparation stage (Illumina stages 1-8).

    Attributes:
        r1_fastq: Path to the R1 paired-end FASTQ (gzipped).
        r2_fastq: Path to the R2 paired-end FASTQ (gzipped).
        intermediate_files: Paths to intermediate files created during
            fragment preparation, suitable for cleanup.
    """

    r1_fastq: str
    r2_fastq: str
    intermediate_files: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class AlignmentResult:
    """Result of the alignment and refinement stage (Illumina stages 9-11).

    Attributes:
        final_bam: Path to the final BAM file (may be aligned, VNTR-biased,
            or downsampled depending on config).
        intermediate_bams: Paths to BAM files superseded during processing.
        intermediate_files: Paths to non-BAM intermediate files (depth files).
    """

    final_bam: str
    intermediate_bams: list[str] = field(default_factory=list)
    intermediate_files: list[str] = field(default_factory=list)
