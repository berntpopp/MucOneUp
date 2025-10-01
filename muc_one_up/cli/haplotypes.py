"""
Haplotype generation functions for MucOneUp CLI.

Single Responsibility: Generate diploid haplotypes from configuration.
"""

import logging
import sys

from ..simulate import simulate_diploid, simulate_from_chains


def generate_haplotypes(
    args, config, fixed_conf, predefined_chains
) -> list[tuple[str, list[str]]]:
    """
    Generate haplotypes based on simulation mode.

    Single Responsibility: Haplotype generation.

    Args:
        args: Parsed command-line arguments
        config: Configuration dictionary
        fixed_conf: Fixed length configuration or "from_structure"
        predefined_chains: Predefined chains from structure file (if any)

    Returns:
        List of (sequence, chain) tuples
    """
    try:
        if fixed_conf == "from_structure":
            results = simulate_from_chains(predefined_chains=predefined_chains, config=config)
        else:
            results = simulate_diploid(
                config=config,
                num_haplotypes=args.num_haplotypes,
                fixed_lengths=fixed_conf,
                seed=args.seed,
            )
        return results
    except Exception as e:
        logging.error("Simulation failed: %s", e)
        sys.exit(1)
