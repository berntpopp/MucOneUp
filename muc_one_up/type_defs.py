"""Type definitions for MucOneUp.

This module provides type aliases and Protocol definitions to improve type safety
and IDE support throughout the codebase.

Following SOLID principles:
- Interface Segregation: Clear type interfaces document expectations
- Dependency Inversion: Protocol classes enable abstract interfaces
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol


@dataclass(frozen=True, slots=True)
class RepeatUnit:
    """A single repeat unit in a VNTR chain.

    Replaces the convention of using plain strings with 'm' suffix
    for mutated repeats. Provides typed access to the symbol and
    mutation status.

    Attributes:
        symbol: The repeat type identifier (e.g., "1", "6p", "X").
        mutated: Whether this repeat has been mutated.
    """

    symbol: str
    mutated: bool = False

    def __str__(self) -> str:
        """Serialize to legacy string format (e.g., 'Xm' if mutated)."""
        return f"{self.symbol}m" if self.mutated else self.symbol

    @classmethod
    def from_str(cls, s: str) -> RepeatUnit:
        """Parse from legacy string format (e.g., 'Xm' -> RepeatUnit('X', True))."""
        if s.endswith("m"):
            return cls(symbol=s[:-1], mutated=True)
        return cls(symbol=s, mutated=False)


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
