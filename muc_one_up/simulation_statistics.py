#!/usr/bin/env python3
"""
simulation_statistics.py

This module implements the generation of comprehensive simulation statistics
for each simulation run. Statistics include:

    - Runtime for the simulation.
    - Haplotype-level summaries:
        * Number of repeats (the VNTR chain length).
        * VNTR region length (in nucleotides).
        * GC content of the simulated sequence.
        * Repeat lengths (each repeat’s nucleotide length, plus summary: min, max, average).
        * Count and breakdown of repeat types.
        * Mutation details (positions and count of mutated repeats).
    - Overall aggregated statistics across haplotypes.
    - Mutation information (if any mutation was applied).
    - VNTR coverage statistics extracted from BAM files (if available).

The report is output as a JSON file and can be integrated into the pipeline’s logging
and reporting systems.
"""

import json
import logging
from typing import Any, Dict, List, Optional

# ---------------------------------------------------------------------------
# Utility functions


def compute_gc_content(seq: str) -> float:
    """
    Compute the GC content percentage of a DNA sequence.

    Args:
        seq (str): DNA sequence.

    Returns:
        float: GC content as a percentage.
    """
    if not seq:
        return 0.0
    seq_upper = seq.upper()
    gc_count = seq_upper.count("G") + seq_upper.count("C")
    return 100.0 * gc_count / len(seq)


def extract_vntr_region(seq: str, config: Dict[str, Any]) -> str:
    """
    Extract the VNTR region from a simulated haplotype sequence by removing
    the left and right constant flanks.

    Args:
        seq (str): The full simulated haplotype sequence.
        config (Dict[str, Any]): Configuration dict containing 'constants'.

    Returns:
        str: The VNTR region sequence.
    """
    reference_assembly = config.get("reference_assembly", "hg38")
    left = config.get("constants", {}).get(reference_assembly, {}).get("left", "")
    right = config.get("constants", {}).get(reference_assembly, {}).get("right", "")
    if seq.startswith(left) and seq.endswith(right):
        return seq[len(left) : -len(right)]
    return seq


def get_repeat_lengths(chain: List[str], config: Dict[str, Any]) -> List[int]:
    """
    Calculate the nucleotide lengths of each repeat in the chain.

    Args:
        chain (List[str]): List of repeat symbols (possibly with a mutation marker).
        config (Dict[str, Any]): Configuration dict containing the 'repeats' sequences.

    Returns:
        List[int]: List of nucleotide lengths for each repeat.
    """
    lengths = []
    repeats_dict = config.get("repeats", {})
    for symbol in chain:
        pure_symbol = symbol.rstrip("m")
        repeat_seq = repeats_dict.get(pure_symbol, "")
        lengths.append(len(repeat_seq))
    return lengths


def count_repeat_types(chain: List[str]) -> Dict[str, int]:
    """
    Count the frequency of each repeat type in the chain.

    Args:
        chain (List[str]): List of repeat symbols (with or without a mutation marker).

    Returns:
        Dict[str, int]: A dictionary mapping each repeat type (without the 'm' marker)
                        to its count.
    """
    counts = {}
    for symbol in chain:
        pure_symbol = symbol.rstrip("m")
        counts[pure_symbol] = counts.get(pure_symbol, 0) + 1
    return counts


def get_mutation_details(chain: List[str]) -> List[Dict[str, Any]]:
    """
    Identify the positions in the chain where mutations have been applied.

    A repeat is considered mutated if its symbol ends with "m".

    Args:
        chain (List[str]): List of repeat symbols.

    Returns:
        List[Dict[str, Any]]: List of mutation details with 1-based positions and the
                              underlying repeat type.
    """
    mutations = []
    for idx, symbol in enumerate(chain, start=1):
        if symbol.endswith("m"):
            mutations.append({"position": idx, "repeat": symbol.rstrip("m")})
    return mutations


# ---------------------------------------------------------------------------
# Haplotype and overall statistics functions


def generate_haplotype_stats(
    simulation_results: List[tuple], config: Dict[str, Any]
) -> List[Dict[str, Any]]:
    """
    Generate statistics for each haplotype simulation result.

    For each haplotype (a tuple of (sequence, chain)), the function calculates:
        - The VNTR region (by removing constant flanks).
        - The number of repeats (chain length).
        - The overall GC content of the full simulated sequence.
        - The individual repeat lengths and their summary (min, max, average).
        - The count of each repeat type.
        - The number of mutated repeats and their positions.

    Args:
        simulation_results (List[tuple]): List of (sequence, chain) for each haplotype.
        config (Dict[str, Any]): Configuration dictionary.

    Returns:
        List[Dict[str, Any]]: A list of dictionaries, one per haplotype, containing
                              statistics.
    """
    hap_stats = []
    for seq, chain in simulation_results:
        vntr_region = extract_vntr_region(seq, config)
        repeat_lengths = get_repeat_lengths(chain, config)
        hap_stat = {
            "repeat_count": len(chain),
            "vntr_length": len(vntr_region),
            "gc_content": compute_gc_content(seq),
            "repeat_lengths": repeat_lengths,
            "repeat_lengths_summary": {
                "min": min(repeat_lengths) if repeat_lengths else 0,
                "max": max(repeat_lengths) if repeat_lengths else 0,
                "average": (
                    sum(repeat_lengths) / len(repeat_lengths) if repeat_lengths else 0
                ),
            },
            "repeat_type_counts": count_repeat_types(chain),
            "mutant_repeat_count": sum(1 for r in chain if r.endswith("m")),
            "mutation_details": get_mutation_details(chain),
            # Initialize SNP fields that will be populated later if SNPs were applied
            "snp_count": 0,
            "applied_snps": [],
        }
        hap_stats.append(hap_stat)
    return hap_stats


def generate_overall_stats(hap_stats: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Generate overall aggregated statistics from haplotype-level data.

    Aggregates include:
        - Repeat count (average, min, max).
        - VNTR region length (average, min, max).
        - GC content (average, min, max).

    Args:
        hap_stats (List[Dict[str, Any]]): List of haplotype statistics dictionaries.

    Returns:
        Dict[str, Any]: Dictionary containing overall statistics.
    """
    repeat_counts = [hs["repeat_count"] for hs in hap_stats]
    vntr_lengths = [hs["vntr_length"] for hs in hap_stats]
    gc_contents = [hs["gc_content"] for hs in hap_stats]

    overall = {
        "repeat_count": {
            "average": sum(repeat_counts) / len(repeat_counts) if repeat_counts else 0,
            "min": min(repeat_counts) if repeat_counts else 0,
            "max": max(repeat_counts) if repeat_counts else 0,
        },
        "vntr_length": {
            "average": sum(vntr_lengths) / len(vntr_lengths) if vntr_lengths else 0,
            "min": min(vntr_lengths) if vntr_lengths else 0,
            "max": max(vntr_lengths) if vntr_lengths else 0,
        },
        "gc_content": {
            "average": sum(gc_contents) / len(gc_contents) if gc_contents else 0,
            "min": min(gc_contents) if gc_contents else 0,
            "max": max(gc_contents) if gc_contents else 0,
        },
    }
    return overall


# ---------------------------------------------------------------------------
# Main statistics report functions


def generate_simulation_statistics(
    start_time: float,
    end_time: float,
    simulation_results: List[tuple],
    config: Dict[str, Any],
    mutation_info: Optional[Dict[str, Any]] = None,
    vntr_coverage: Optional[Dict[str, Any]] = None,
    applied_snp_info: Optional[Dict[int, List[Dict[str, Any]]]] = None,
) -> Dict[str, Any]:
    """
    Generate a comprehensive simulation statistics report.

    The report includes:
      - Simulation runtime (in seconds).
      - Haplotype-level statistics (repeat counts, VNTR lengths, GC content,
        repeat lengths and summaries, repeat type counts, and mutation details).
      - Overall aggregated statistics.
      - Mutation information (if any).
      - VNTR coverage statistics (if available).

    Args:
        start_time (float): Simulation start time (timestamp).
        end_time (float): Simulation end time (timestamp).
        simulation_results (List[tuple]): List of (sequence, chain) for each haplotype.
        config (Dict[str, Any]): Configuration dictionary.
        mutation_info (Optional[Dict[str, Any]]): Mutation details (e.g., mutation name,
            target positions). Defaults to None.
        vntr_coverage (Optional[Dict[str, Any]]): VNTR coverage statistics (e.g., from BAM files).
            Defaults to None.
        applied_snp_info (Optional[Dict[int, List[Dict[str, Any]]]]): Dictionary mapping haplotype
            index (0-based) to list of successfully applied SNPs. Defaults to None.

    Returns:
        Dict[str, Any]: Comprehensive simulation statistics.
    """
    runtime = end_time - start_time
    haplotype_stats = generate_haplotype_stats(simulation_results, config)
    overall_stats = generate_overall_stats(haplotype_stats)

    # Get the reference assembly to include in the report
    reference_assembly = config.get("reference_assembly", "hg38")

    # Add SNP information to haplotype statistics if available
    if applied_snp_info:
        # Convert to 1-indexed for consistency with other parts of the report
        snp_info_1indexed = {}
        for hap_idx, snps in applied_snp_info.items():
            # Convert each SNP dictionary to a serializable form
            # (Removing any potential complex objects like objects with __dict__ etc.)
            serializable_snps = []
            for snp in snps:
                serializable_snp = {
                    "position": snp.get("position"),
                    "ref_base": snp.get("ref_base"),
                    "alt_base": snp.get("alt_base"),
                }
                serializable_snps.append(serializable_snp)

            # Store with 1-indexed haplotype
            snp_info_1indexed[hap_idx + 1] = serializable_snps

            # Add SNP count to the corresponding haplotype statistics
            if 0 <= hap_idx < len(haplotype_stats):
                haplotype_stats[hap_idx]["snp_count"] = len(snps)
                haplotype_stats[hap_idx]["applied_snps"] = serializable_snps

    report = {
        "runtime_seconds": runtime,
        "reference_assembly": reference_assembly,
        "haplotype_statistics": haplotype_stats,
        "overall_statistics": overall_stats,
        "mutation_info": mutation_info if mutation_info is not None else {},
        "vntr_coverage": vntr_coverage if vntr_coverage is not None else {},
        "snp_info": snp_info_1indexed if applied_snp_info else {},
    }
    return report


def write_statistics_report(report: Dict[str, Any], output_path: str) -> None:
    """
    Write the simulation statistics report to a JSON file.

    Args:
        report (Dict[str, Any]): Simulation statistics report.
        output_path (str): File path for the JSON report.
    """
    try:
        with open(output_path, "w") as fh:
            json.dump(report, fh, indent=4)
        logging.info("Simulation statistics report written to %s", output_path)
    except Exception as exc:
        logging.error("Failed to write statistics report: %s", exc)
        raise


# ---------------------------------------------------------------------------
# Example usage when run as a script

if __name__ == "__main__":
    import time

    # Dummy configuration for demonstration.
    dummy_config = {
        "constants": {"left": "TTTT", "right": "AAAA"},
        "repeats": {"1": "AAA", "2": "CCC", "9": "GGG"},
    }
    # Dummy simulation results: two haplotypes (sequence, chain).
    dummy_results = [
        ("TTTTAAACCCGGGAAAA", ["1", "2", "9"]),
        ("TTTTAAACCCGGGAAAA", ["1", "2", "9m"]),
    ]

    # Simulate a run by waiting a bit.
    start = time.time()
    time.sleep(1.2)
    end = time.time()

    # Dummy mutation information and VNTR coverage statistics.
    dummy_mutation = {"mutation_name": "dupC", "targets": [(2, 3)]}
    dummy_vntr_coverage = {"mean_coverage": 30.5, "region": "chr1:1000-2000"}

    # Generate and print the report.
    stats_report = generate_simulation_statistics(
        start_time=start,
        end_time=end,
        simulation_results=dummy_results,
        config=dummy_config,
        mutation_info=dummy_mutation,
        vntr_coverage=dummy_vntr_coverage,
    )
    print(json.dumps(stats_report, indent=4))
