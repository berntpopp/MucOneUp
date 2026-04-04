"""Samtools wrapper — re-exports from split modules.

All functions are available from this module for backward compatibility.
New code should import from specific submodules.
"""

from .samtools_convert import (
    FastqConversionOptions,
    _count_fastq_reads,
    convert_bam_to_fastq,
    convert_bam_to_paired_fastq,
    convert_sam_to_bam,
)
from .samtools_core import (
    extract_subset_reference,
    merge_bam_files,
    sort_and_index_bam,
)
from .samtools_coverage import (
    calculate_target_coverage,
    calculate_vntr_coverage,
    downsample_bam,
    downsample_entire_bam,
)

__all__ = [
    "FastqConversionOptions",
    "_count_fastq_reads",
    "calculate_target_coverage",
    "calculate_vntr_coverage",
    "convert_bam_to_fastq",
    "convert_bam_to_paired_fastq",
    "convert_sam_to_bam",
    "downsample_bam",
    "downsample_entire_bam",
    "extract_subset_reference",
    "merge_bam_files",
    "sort_and_index_bam",
]
