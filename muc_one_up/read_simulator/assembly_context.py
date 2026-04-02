"""Centralized assembly resolution for read simulation pipelines.

Resolves assembly name, constants, and reference paths once at pipeline
entry. Eliminates repeated config.get("reference_assembly", "hg38")
calls scattered across pipeline stages.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any

logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class AssemblyContext:
    """Resolved assembly configuration for a pipeline run.

    Attributes:
        assembly_name: Reference assembly identifier (e.g., "hg38", "hg19").
        left_constant: Left flanking constant sequence.
        right_constant: Right flanking constant sequence.
        sample_bam: Path to sample BAM for the assembly (optional).
        vntr_region: VNTR region string for the assembly (optional).
        human_reference: Path to human reference FASTA (optional).
    """

    assembly_name: str
    left_constant: str
    right_constant: str
    sample_bam: str | None = None
    vntr_region: str | None = None
    human_reference: str | None = None

    @classmethod
    def from_configs(
        cls,
        config: dict[str, Any],
        rs_config: dict[str, Any] | None = None,
    ) -> AssemblyContext:
        """Construct from simulation config and read-simulator config.

        Args:
            config: Main simulation config with 'reference_assembly' and 'constants'.
            rs_config: Read-simulator config with assembly-keyed paths (optional).

        Returns:
            Resolved AssemblyContext.
        """
        rs_config = rs_config or {}
        assembly = config.get("reference_assembly", "hg38")
        constants = config.get("constants", {}).get(assembly, {})

        if not constants:
            logger.debug("No constants found for assembly '%s' in config", assembly)

        return cls(
            assembly_name=assembly,
            left_constant=constants.get("left", ""),
            right_constant=constants.get("right", ""),
            sample_bam=rs_config.get(f"sample_bam_{assembly}"),
            vntr_region=rs_config.get(f"vntr_region_{assembly}"),
            human_reference=rs_config.get(f"human_reference_{assembly}"),
        )
