"""General validation functions for MucOneUp.

This module provides runtime validation for inputs, following the principle
of fail-fast: catch errors early with clear error messages.

Functions in this module:
- Validate data structures (haplotype indices, mutation targets)
- Validate file paths and existence
- Validate configuration values
"""

from pathlib import Path

from .exceptions import FileOperationError, ValidationError
from .type_defs import ConfigDict, MutationTargets


def validate_haplotype_index(index: int, num_haplotypes: int) -> None:
    """Validate haplotype index is within valid range.

    Args:
        index: 0-based haplotype index
        num_haplotypes: Total number of haplotypes

    Raises:
        ValidationError: If index is out of range
    """
    if not 0 <= index < num_haplotypes:
        raise ValidationError(
            f"Haplotype index {index} out of range. "
            f"Must be 0-{num_haplotypes - 1} (0-based indexing)."
        )


def validate_repeat_index(index: int, chain_length: int, context: str = "") -> None:
    """Validate repeat index is within chain bounds.

    Args:
        index: 0-based repeat index
        chain_length: Length of repeat chain
        context: Optional context for error message (e.g., "mutation application")

    Raises:
        ValidationError: If index is out of range
    """
    if not 0 <= index < chain_length:
        msg = f"Repeat index {index} out of range (chain length: {chain_length})."
        if context:
            msg = f"{context}: {msg}"
        raise ValidationError(msg)


def validate_file_exists(filepath: str | Path, description: str = "File") -> None:
    """Validate that a file exists.

    Args:
        filepath: Path to file
        description: Description for error message (e.g., "Config file", "FASTA file")

    Raises:
        FileOperationError: If file doesn't exist
    """
    path = Path(filepath)
    if not path.exists():
        raise FileOperationError(f"{description} not found: {filepath}")


def validate_directory_exists(dirpath: str | Path, description: str = "Directory") -> None:
    """Validate that a directory exists.

    Args:
        dirpath: Path to directory
        description: Description for error message

    Raises:
        FileOperationError: If directory doesn't exist
    """
    path = Path(dirpath)
    if not path.exists():
        raise FileOperationError(f"{description} not found: {dirpath}")
    if not path.is_dir():
        raise FileOperationError(f"{description} is not a directory: {dirpath}")


def validate_mutation_exists(mutation_name: str, config: ConfigDict) -> None:
    """Validate that a mutation is defined in config.

    Args:
        mutation_name: Name of mutation to check
        config: Configuration dictionary

    Raises:
        ValidationError: If mutation is not defined
    """
    mutations = config.get("mutations", {})
    if mutation_name not in mutations:
        available = list(mutations.keys())
        raise ValidationError(
            f"Mutation '{mutation_name}' not defined. " f"Available mutations: {available}"
        )


def validate_mutation_targets(
    targets: MutationTargets,
    num_haplotypes: int,
    chains: list[list[str]],
) -> None:
    """Validate mutation targets are within bounds.

    Args:
        targets: List of (haplotype_idx, repeat_idx) tuples (1-based)
        num_haplotypes: Number of haplotypes
        chains: List of repeat chains

    Raises:
        ValidationError: If any target is invalid
    """
    for hap_idx, repeat_idx in targets:
        # Validate haplotype index (1-based)
        if not 1 <= hap_idx <= num_haplotypes:
            raise ValidationError(
                f"Mutation target haplotype index {hap_idx} out of range. "
                f"Must be 1-{num_haplotypes} (1-based indexing)."
            )

        # Validate repeat index (1-based)
        chain = chains[hap_idx - 1]  # Convert to 0-based
        if not 1 <= repeat_idx <= len(chain):
            raise ValidationError(
                f"Mutation target repeat index {repeat_idx} out of range "
                f"for haplotype {hap_idx} (chain length: {len(chain)}). "
                f"Must be 1-{len(chain)} (1-based indexing)."
            )


def validate_repeat_symbol(symbol: str, valid_symbols: set[str]) -> None:
    """Validate repeat symbol is in valid set.

    Args:
        symbol: Repeat symbol to validate (e.g., '1', '2', 'X')
        valid_symbols: Set of valid repeat symbols

    Raises:
        ValidationError: If symbol is invalid
    """
    # Remove mutation marker if present
    clean_symbol = symbol.rstrip("m")

    if clean_symbol not in valid_symbols:
        raise ValidationError(
            f"Invalid repeat symbol: '{symbol}'. " f"Valid symbols: {sorted(valid_symbols)}"
        )


def validate_positive_integer(value: int, name: str) -> None:
    """Validate value is a positive integer.

    Args:
        value: Value to validate
        name: Name of value for error message

    Raises:
        ValidationError: If value is not positive
    """
    if not isinstance(value, int) or value <= 0:
        raise ValidationError(f"{name} must be a positive integer, got: {value}")


def validate_non_negative_integer(value: int, name: str) -> None:
    """Validate value is a non-negative integer.

    Args:
        value: Value to validate
        name: Name of value for error message

    Raises:
        ValidationError: If value is negative
    """
    if not isinstance(value, int) or value < 0:
        raise ValidationError(f"{name} must be non-negative, got: {value}")


def validate_probability(value: float, name: str) -> None:
    """Validate value is a valid probability (0.0-1.0).

    Args:
        value: Value to validate
        name: Name of value for error message

    Raises:
        ValidationError: If value is not in range [0.0, 1.0]
    """
    if not isinstance(value, int | float) or not 0.0 <= value <= 1.0:
        raise ValidationError(f"{name} must be between 0.0 and 1.0, got: {value}")


def validate_reference_assembly(assembly: str) -> None:
    """Validate reference assembly is supported.

    Args:
        assembly: Reference assembly name (e.g., 'hg19', 'hg38')

    Raises:
        ValidationError: If assembly is not supported
    """
    valid_assemblies = {"hg19", "hg38"}
    if assembly not in valid_assemblies:
        raise ValidationError(
            f"Invalid reference assembly: '{assembly}'. "
            f"Supported assemblies: {sorted(valid_assemblies)}"
        )
