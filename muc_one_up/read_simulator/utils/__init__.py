"""
Utility modules for read simulation.

This package contains specialized utility modules for handling
reference files, FASTQ operations, and diploid-specific simulation logic.
"""

# Import from common_utils for backward compatibility
from .common_utils import (
    check_external_tools,
    cleanup_files,
    fix_field,
    run_command,
)
from .diploid_handler import (
    DiploidSimulationResult,
    calculate_corrected_coverage,
    prepare_diploid_simulation,
    run_split_simulation,
)
from .fastq_utils import (
    count_fastq_reads,
    merge_fastq_files,
    validate_fastq,
)

# Import from new specialized modules
from .reference_utils import (
    extract_haplotypes,
    get_reference_info,
    is_diploid_reference,
    validate_reference_compatibility,
)

__all__ = [
    # Diploid handler
    "DiploidSimulationResult",
    "calculate_corrected_coverage",
    # Common utilities (backward compatibility)
    "check_external_tools",
    "cleanup_files",
    # FASTQ utilities
    "count_fastq_reads",
    # Reference utilities
    "extract_haplotypes",
    "fix_field",
    "get_reference_info",
    "is_diploid_reference",
    "merge_fastq_files",
    "prepare_diploid_simulation",
    "run_command",
    "run_split_simulation",
    "validate_fastq",
    "validate_reference_compatibility",
]
