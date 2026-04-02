"""Centralized sequence assembly from repeat chains.

Single source of truth for chain -> DNA sequence conversion.
Replaces assemble_haplotype_from_chain() in simulate.py and
rebuild_haplotype_sequence() in mutate.py.
"""

from __future__ import annotations

import logging

from .type_defs import ConfigDict, DNASequence, RepeatUnit

logger = logging.getLogger(__name__)


def assemble_sequence(chain: list[RepeatUnit], config: ConfigDict) -> DNASequence:
    """Assemble a DNA sequence from a repeat chain and flanking constants.

    Concatenates: left_constant + repeat_units + right_constant.
    The right constant is only appended if the chain ends with repeat "9"
    (the canonical terminal repeat in MUC1 VNTR).

    Each RepeatUnit carries its symbol directly, so no string manipulation
    is needed to resolve the base repeat type.

    Args:
        chain: List of RepeatUnit objects describing the repeat chain.
        config: Configuration dict with 'repeats', 'constants', and
                optionally 'reference_assembly' keys.

    Returns:
        Assembled DNA sequence string.

    Raises:
        KeyError: If a repeat symbol is not in config repeats,
            or if required constants for the reference assembly are missing.
    """
    repeats_dict = config["repeats"]
    ref_assembly = config.get("reference_assembly", "hg38")
    constants = config["constants"][ref_assembly]
    left_const = constants["left"]
    right_const = constants["right"]

    # Build repeat region
    parts: list[str] = [left_const]
    for unit in chain:
        if unit.symbol not in repeats_dict:
            raise KeyError(f"Repeat symbol '{unit.symbol}' not found in config repeats")
        parts.append(repeats_dict[unit.symbol])

    # Right constant only if chain ends with canonical terminal repeat "9"
    if chain and chain[-1].symbol == "9":
        parts.append(right_const)
    elif chain:
        logger.debug("Chain does not end with '9'; right constant omitted.")

    return "".join(parts)
