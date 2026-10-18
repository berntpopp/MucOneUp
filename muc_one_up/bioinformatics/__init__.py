"""Bioinformatics-specific reference validation and utilities for MucOneUp.

This package contains modules for:
- Reference genome validation
"""

from .reference_validation import (
    get_muc1_region_for_assembly,
    get_reference_path_for_assembly,
    validate_bam_file,
    validate_bed_file,
    validate_reference_for_assembly,
    validate_reference_genome,
)

__all__ = [
    "get_muc1_region_for_assembly",
    "get_reference_path_for_assembly",
    "validate_bam_file",
    "validate_bed_file",
    "validate_reference_for_assembly",
    "validate_reference_genome",
]
