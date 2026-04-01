"""Centralized sequence assembly from repeat chains.

Single source of truth for chain -> DNA sequence conversion.
Replaces assemble_haplotype_from_chain() in simulate.py and
rebuild_haplotype_sequence() in mutate.py.
"""

from __future__ import annotations

import logging

from .type_defs import ConfigDict, DNASequence, RepeatChain

logger = logging.getLogger(__name__)


def assemble_sequence(chain: RepeatChain, config: ConfigDict) -> DNASequence:
    """Assemble a DNA sequence from a repeat chain and flanking constants.

    Concatenates: left_constant + repeat_units + right_constant.
    The right constant is only appended if the chain ends with repeat "9"
    (the canonical terminal repeat in MUC1 VNTR).

    Mutation markers (trailing 'm' on symbols like 'Xm') are stripped
    before looking up repeat sequences, so mutated and unmutated repeats
    use the same base sequence.

    Args:
        chain: List of repeat symbols (e.g., ["1", "2", "Xm", "7", "8", "9"]).
        config: Configuration dict with 'repeats', 'constants', and
                optionally 'reference_assembly' keys.

    Returns:
        Assembled DNA sequence string.

    Raises:
        KeyError: If a repeat symbol (after stripping 'm') is not in config.
    """
    repeats_dict = config["repeats"]
    ref_assembly = config.get("reference_assembly", "hg38")
    constants = config.get("constants", {}).get(ref_assembly, {})
    left_const = constants.get("left", "")
    right_const = constants.get("right", "")

    # Build repeat region
    parts: list[str] = [left_const]
    for symbol in chain:
        base_symbol = symbol.rstrip("m")
        if base_symbol not in repeats_dict:
            raise KeyError(
                f"Repeat symbol '{base_symbol}' not found in config repeats"
            )
        parts.append(repeats_dict[base_symbol])

    # Right constant only if chain ends with canonical terminal repeat "9"
    if chain and chain[-1].rstrip("m") == "9":
        parts.append(right_const)
    elif chain:
        logger.debug("Chain does not end with '9'; right constant omitted.")

    return "".join(parts)
