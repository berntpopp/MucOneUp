"""Type definitions for MucOneUp.

This module provides type aliases and Protocol definitions to improve type safety
and IDE support throughout the codebase.

Following SOLID principles:
- Interface Segregation: Clear type interfaces document expectations
- Dependency Inversion: Protocol classes enable abstract interfaces
"""

from pathlib import Path
from typing import Any, Protocol

# Type aliases for haplotype representation
HaplotypeName = str
RepeatChain = list[str]
Haplotype = tuple[HaplotypeName, RepeatChain]
HaplotypeList = list[Haplotype]

# Configuration types
ConfigDict = dict[str, Any]
RepeatsDict = dict[str, str]
ProbabilitiesDict = dict[str, dict[str, float]]
ConstantsDict = dict[str, dict[str, str]]
LengthModelDict = dict[str, Any]

# Mutation types
MutationName = str
MutationTargets = list[tuple[int, int]]
MutatedUnits = dict[str, str]
MutationChange = dict[str, int | str]
MutationDefinition = dict[str, Any]

# SNP types
SNPRecord = tuple[int, int, str, str]
SNPList = list[SNPRecord]
SNPRegion = str

# File path types
FilePath = str | Path

# Sequence types
DNASequence = str
ProteinSequence = str
RepeatStructure = str

# Statistics types
SimulationStats = dict[str, Any]
HaplotypeStats = dict[str, Any]

# Read simulation types
ReadSimulationType = str
AlignmentFormat = str


# Protocol definitions for dependency injection
class ToolWrapper(Protocol):
    """Protocol for external tool wrappers."""

    def run_command(self, command: list[str], timeout: int | None = None) -> tuple[str, str]:
        """Run a command and return (stdout, stderr)."""
        ...


class ConfigLoader(Protocol):
    """Protocol for configuration loaders."""

    def load_config(self, config_path: str) -> ConfigDict:
        """Load and validate configuration."""
        ...


class SequenceValidator(Protocol):
    """Protocol for sequence validation."""

    def validate_dna_sequence(self, sequence: str, allow_ambiguous: bool = True) -> None:
        """Validate DNA sequence."""
        ...
