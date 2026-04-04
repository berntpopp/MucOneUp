"""Centralized sequence assembly from repeat chains.

Single source of truth for chain -> DNA sequence conversion.
Replaces assemble_haplotype_from_chain() in simulate.py and
rebuild_haplotype_sequence() in mutate.py.
"""

from __future__ import annotations

import logging

from .type_defs import AssemblyConstants, DNASequence, RepeatUnit

logger = logging.getLogger(__name__)


def assemble_sequence(
    chain: list[RepeatUnit],
    repeats: dict[str, str],
    constants: AssemblyConstants,
) -> DNASequence:
    """Assemble a DNA sequence from a repeat chain and flanking constants.

    Concatenates: left_constant + repeat_units + right_constant.
    The right constant is only appended if the chain ends with repeat "9"
    (the canonical terminal repeat in MUC1 VNTR).

    Each RepeatUnit carries its symbol directly, so no string manipulation
    is needed to resolve the base repeat type.

    Args:
        chain: List of RepeatUnit objects describing the repeat chain.
        repeats: Mapping of repeat symbol to DNA sequence.
        constants: Assembly constants with 'left' and 'right' flanking
                   sequences.

    Returns:
        Assembled DNA sequence string.

    Raises:
        KeyError: If a repeat symbol is not found in repeats.
    """
    left_const = constants["left"]
    right_const = constants["right"]

    # Build repeat region
    parts: list[str] = [left_const]
    for unit in chain:
        if unit.symbol not in repeats:
            raise KeyError(f"Repeat symbol '{unit.symbol}' not found in repeats")
        parts.append(repeats[unit.symbol])

    # Right constant only if chain ends with canonical terminal repeat "9"
    if chain and chain[-1].symbol == "9":
        parts.append(right_const)
    elif chain:
        logger.debug("Chain does not end with '9'; right constant omitted.")

    return "".join(parts)
