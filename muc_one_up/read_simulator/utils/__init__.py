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

# Import metadata and tool version tracking
from .metadata_writer import write_metadata_file

# Import from new specialized modules
from .reference_utils import (
    extract_haplotypes,
    get_reference_info,
    is_diploid_reference,
    validate_reference_compatibility,
)
from .tool_version import (
    capture_tool_versions,
    get_tool_version,
    log_tool_versions,
)

__all__ = [
    "DiploidSimulationResult",
    "calculate_corrected_coverage",
    "capture_tool_versions",
    "check_external_tools",
    "cleanup_files",
    "count_fastq_reads",
    "extract_haplotypes",
    "fix_field",
    "get_reference_info",
    "get_tool_version",
    "is_diploid_reference",
    "log_tool_versions",
    "merge_fastq_files",
    "prepare_diploid_simulation",
    "run_command",
    "run_split_simulation",
    "validate_fastq",
    "validate_reference_compatibility",
    "write_metadata_file",
]
