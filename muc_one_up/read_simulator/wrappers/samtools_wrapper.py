"""Samtools wrapper — re-exports from split modules.

All functions are available from this module for backward compatibility.
New code should import from specific submodules.
"""

from ..utils.samtools import (
    SamtoolsError,
    check_samtools_available,
    get_read_count,
    validate_bam,
)
from .samtools_convert import (
    FastqConversionOptions,
    convert_bam_to_fastq,
    convert_bam_to_paired_fastq,
    convert_sam_to_bam,
)
from .samtools_core import (
    extract_subset_reference,
    merge_bam_files,
)
from .samtools_coverage import (
    calculate_target_coverage,
    calculate_vntr_coverage,
    downsample_bam,
    downsample_entire_bam,
)

__all__ = [
    "FastqConversionOptions",
    "SamtoolsError",
    "calculate_target_coverage",
    "calculate_vntr_coverage",
    "check_samtools_available",
    "convert_bam_to_fastq",
    "convert_bam_to_paired_fastq",
    "convert_sam_to_bam",
    "downsample_bam",
    "downsample_entire_bam",
    "extract_subset_reference",
    "get_read_count",
    "merge_bam_files",
    "validate_bam",
]
