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


@dataclass(slots=True)
class HaplotypeResult:
    """Result of simulating a single haplotype.

    Replaces the raw tuple[str, list[str]] convention with named fields.

    Attributes:
        sequence: The assembled DNA sequence.
        chain: The repeat chain as typed RepeatUnit objects.
    """

    sequence: str
    chain: list[RepeatUnit]

    def chain_strs(self) -> list[str]:
        """Return chain as legacy string list (e.g., ['1', 'Xm', '9'])."""
        return [str(ru) for ru in self.chain]

    def as_tuple(self) -> tuple[str, list[str]]:
        """Convert to legacy (sequence, chain_strs) tuple."""
        return (self.sequence, self.chain_strs())

    @classmethod
    def from_tuple(cls, t: tuple[str, list[str]]) -> HaplotypeResult:
        """Parse from legacy (sequence, chain_strs) tuple."""
        return cls(
            sequence=t[0],
            chain=[RepeatUnit.from_str(s) for s in t[1]],
        )


@dataclass(frozen=True, slots=True)
class MutationTarget:
    """A target position for mutation application.

    Replaces raw tuple[int, int] with named, validated fields.
    Uses 1-based indexing matching biological conventions.

    Attributes:
        haplotype_index: 1-based haplotype number.
        repeat_index: 1-based repeat position within the chain.
    """

    haplotype_index: int
    repeat_index: int

    def __post_init__(self) -> None:
        if self.haplotype_index < 1 or self.repeat_index < 1:
            raise ValueError(
                f"MutationTarget uses 1-based indexing, got "
                f"haplotype_index={self.haplotype_index}, "
                f"repeat_index={self.repeat_index}"
            )

    def as_tuple(self) -> tuple[int, int]:
        """Convert to legacy (haplotype_index, repeat_index) tuple."""
        return (self.haplotype_index, self.repeat_index)

    @classmethod
    def from_tuple(cls, t: tuple[int, int]) -> MutationTarget:
        """Parse from legacy (haplotype_index, repeat_index) tuple."""
        return cls(haplotype_index=t[0], repeat_index=t[1])


# Legacy type alias — used by thin wrappers accepting raw string chains
RepeatChain = list[str]

# Configuration types
ConfigDict = dict[str, Any]
RepeatsDict = dict[str, str]
ProbabilitiesDict = dict[str, dict[str, float]]
ConstantsDict = dict[str, dict[str, str]]
LengthModelDict = dict[str, Any]

# Mutation types
MutationName = str
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
