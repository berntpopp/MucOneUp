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
        * Repeat lengths (each repeat`s nucleotide length, plus summary: min, max, average).
        * Count and breakdown of repeat types.
        * Mutation details (positions and count of mutated repeats).
    - Overall aggregated statistics across haplotypes.
    - Mutation information (if any mutation was applied).
    - VNTR coverage statistics extracted from BAM files (if available).

The report is output as a JSON file and can be integrated into the pipeline`s logging
and reporting systems.
"""

import json
import logging
from pathlib import Path
from typing import Any

from .type_defs import HaplotypeResult, RepeatUnit

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


def get_haplotype_vntr_bounds(
    seq: str,
    config: dict[str, Any],
    chain: list[RepeatUnit] | None = None,
    repeat_lengths: list[int] | None = None,
) -> tuple[int, int]:
    """Calculate the [start, end) nucleotide bounds of the VNTR region in seq.

    Args:
        seq: Haplotype DNA sequence.
        config: Configuration dict with constants and optional reference assembly.
        chain: Optional repeat chain.
        repeat_lengths: Optional nucleotide lengths for each repeat in chain.

    Returns:
        tuple[int, int]: (start, end) 0-based boundary indices within seq.
    """
    reference_assembly = config.get("reference_assembly", "hg38")
    constants = config.get("constants", {})
    if isinstance(constants, dict):
        if reference_assembly in constants and isinstance(constants[reference_assembly], dict):
            ref_const = constants[reference_assembly]
        else:
            ref_const = constants
    else:
        ref_const = {}
    left = ref_const.get("left", "") if isinstance(ref_const, dict) else ""
    right = ref_const.get("right", "") if isinstance(ref_const, dict) else ""
    left = str(left) if left is not None else ""
    right = str(right) if right is not None else ""

    if not left and not right:
        return (0, len(seq))

    if chain is None:
        if len(seq) < len(left) + len(right):
            return (0, len(seq))
        return (len(left), len(seq) - len(right))

    if len(seq) < len(left):
        return (0, len(seq))

    vntr_start = len(left)
    if len(right) == 0:
        right_len = 0
    elif not chain or chain[-1].symbol == "9":
        has_right = len(seq) >= len(left) + len(right)
        right_len = len(right) if has_right else 0
    else:
        repeats_dict = config.get("repeats", {})
        expected_vntr_len = (
            sum(repeat_lengths)
            if repeat_lengths is not None
            else sum(len(repeats_dict.get(u.symbol, "")) for u in chain)
        )
        has_right = seq.endswith(right) or len(seq) >= len(left) + expected_vntr_len + len(right)
        right_len = len(right) if has_right else 0

    vntr_end = max(vntr_start, len(seq) - right_len)
    return (vntr_start, vntr_end)


def extract_vntr_region(
    seq: str,
    config: dict[str, Any],
    chain: list[RepeatUnit] | None = None,
    repeat_lengths: list[int] | None = None,
) -> str:
    """Extract the VNTR region from a simulated haplotype sequence.

    Args:
        seq: The full simulated haplotype sequence.
        config: Configuration dict containing 'constants'.
        chain: Repeat chain if available.
        repeat_lengths: Optional precomputed repeat lengths.

    Returns:
        The VNTR region sequence.
    """
    start, end = get_haplotype_vntr_bounds(seq, config, chain=chain, repeat_lengths=repeat_lengths)
    return seq[start:end]


def get_repeat_lengths(
    chain: list[RepeatUnit],
    config: dict[str, Any],
    hap_mutated_units: list[tuple[int, str]] | dict[int, str] | None = None,
    vntr_length: int | None = None,
) -> list[int]:
    """Calculate the nucleotide lengths of each repeat in the chain.

    Args:
        chain: List of RepeatUnit objects.
        config (Dict[str, Any]): Configuration dict containing the 'repeats' sequences.
        hap_mutated_units: Mutated units for this haplotype as list of (rep_idx, seq) or dict.
        vntr_length: Extracted VNTR region length to infer single-indel delta.

    Returns:
        List[int]: List of nucleotide lengths for each repeat.
    """
    repeats_dict = config.get("repeats", {})
    lengths = [len(repeats_dict.get(unit.symbol, "")) for unit in chain]

    if hap_mutated_units:
        items = (
            hap_mutated_units.items() if isinstance(hap_mutated_units, dict) else hap_mutated_units
        )
        for rep_idx, seq in items:
            if 1 <= rep_idx <= len(lengths):
                lengths[rep_idx - 1] = len(seq)
    elif vntr_length is not None:
        mutated_indices = [i for i, u in enumerate(chain) if u.mutated]
        if len(mutated_indices) == 1 and vntr_length >= 0:
            delta = vntr_length - sum(lengths)
            if lengths[mutated_indices[0]] + delta >= 0:
                lengths[mutated_indices[0]] += delta

    return lengths


def count_repeat_types(chain: list[RepeatUnit]) -> dict[str, int]:
    """
    Count the frequency of each repeat type in the chain.

    Args:
        chain: List of RepeatUnit objects.

    Returns:
        Dict[str, int]: A dictionary mapping each repeat type to its count.
    """
    counts: dict[str, int] = {}
    for unit in chain:
        counts[unit.symbol] = counts.get(unit.symbol, 0) + 1
    return counts


def get_mutation_details(chain: list[RepeatUnit]) -> list[dict[str, Any]]:
    """
    Identify the positions in the chain where mutations have been applied.

    Args:
        chain: List of RepeatUnit objects.

    Returns:
        List[Dict[str, Any]]: List of mutation details with 1-based positions and the
                              underlying repeat type.
    """
    mutations = []
    for idx, unit in enumerate(chain, start=1):
        if unit.mutated:
            mutations.append({"position": idx, "repeat": unit.symbol})
    return mutations


# ---------------------------------------------------------------------------
# Haplotype and overall statistics functions


def generate_haplotype_stats(
    simulation_results: list[HaplotypeResult],
    config: dict[str, Any],
    mutated_units: dict[int, list[tuple[int, str]]] | None = None,
) -> list[dict[str, Any]]:
    """
    Generate statistics for each haplotype simulation result.

    For each HaplotypeResult, the function calculates:
        - The VNTR region (by removing constant flanks).
        - The number of repeats (chain length).
        - The overall GC content of the full simulated sequence.
        - The individual repeat lengths and their summary (min, max, average).
        - The count of each repeat type.
        - The number of mutated repeats and their positions.

    Args:
        simulation_results: List of HaplotypeResult objects.
        config (Dict[str, Any]): Configuration dictionary.
        mutated_units: Mapping of haplotype index to mutated unit tuples.

    Returns:
        List[Dict[str, Any]]: A list of dictionaries, one per haplotype, containing
                              statistics.
    """
    hap_stats = []
    for idx, hr in enumerate(simulation_results):
        hr_lengths = getattr(hr, "_repeat_lengths", None)
        vntr_region = extract_vntr_region(
            hr.sequence, config, chain=hr.chain, repeat_lengths=hr_lengths
        )
        hap_mut = (mutated_units or {}).get(idx + 1)
        if hap_mut is None and hr_lengths is not None:
            repeat_lengths = list(hr_lengths)
        else:
            repeat_lengths = get_repeat_lengths(
                hr.chain,
                config,
                hap_mutated_units=hap_mut,
                vntr_length=len(vntr_region),
            )
        hap_stat = {
            "repeat_count": len(hr.chain),
            "vntr_length": len(vntr_region),
            "gc_content": compute_gc_content(hr.sequence),
            "repeat_lengths": repeat_lengths,
            "repeat_lengths_summary": {
                "min": min(repeat_lengths) if repeat_lengths else 0,
                "max": max(repeat_lengths) if repeat_lengths else 0,
                "average": (sum(repeat_lengths) / len(repeat_lengths) if repeat_lengths else 0),
            },
            "repeat_type_counts": count_repeat_types(hr.chain),
            "mutant_repeat_count": sum(1 for ru in hr.chain if ru.mutated),
            "mutation_details": get_mutation_details(hr.chain),
            # Initialize SNP fields that will be populated later if SNPs were applied
            "snp_count": 0,
            "applied_snps": [],
        }
        hap_stats.append(hap_stat)
    return hap_stats


def generate_overall_stats(hap_stats: list[dict[str, Any]]) -> dict[str, Any]:
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
    simulation_results: list[HaplotypeResult],
    config: dict[str, Any],
    mutation_info: dict[str, Any] | None = None,
    vntr_coverage: dict[str, Any] | None = None,
    applied_snp_info: dict[int, list[dict[str, Any]]] | None = None,
    provenance_info: dict[str, Any] | None = None,
    mutated_units: dict[int, list[tuple[int, str]]] | None = None,
) -> dict[str, Any]:
    """
    Generate a comprehensive simulation statistics report.

    The report includes:
      - Simulation runtime (in seconds).
      - Haplotype-level statistics (repeat counts, VNTR lengths, GC content,
        repeat lengths and summaries, repeat type counts, and mutation details).
      - Overall aggregated statistics.
      - Mutation information (if any).
      - VNTR coverage statistics (if available).
      - Provenance metadata (software version, config fingerprint, seed, timestamps).

    Args:
        start_time (float): Simulation start time (timestamp).
        end_time (float): Simulation end time (timestamp).
        simulation_results: List of HaplotypeResult objects.
        config (Dict[str, Any]): Configuration dictionary.
        mutation_info (Optional[Dict[str, Any]]): Mutation details (e.g., mutation name,
            target positions). Defaults to None.
        vntr_coverage (Optional[Dict[str, Any]]): VNTR coverage statistics (e.g., from BAM files).
            Defaults to None.
        applied_snp_info (Optional[Dict[int, List[Dict[str, Any]]]]): Dictionary mapping haplotype
            index (0-based) to list of successfully applied SNPs. Defaults to None.
        provenance_info (Optional[Dict[str, Any]]): Provenance metadata (version, fingerprint,
            seed, timestamps). Defaults to None. (New in v0.28.0)
        mutated_units (Optional[Dict[int, List[Tuple[int, str]]]]): Mapping of haplotype index
            to mutated units. Defaults to None.

    Returns:
        Dict[str, Any]: Comprehensive simulation statistics.

    Note:
        The 'runtime_seconds' field is computed from provenance.duration_seconds
        if provenance is provided (single source of truth). For backward compatibility,
        it falls back to end_time - start_time if provenance is missing.
    """
    # Compute runtime from provenance if available (single source of truth)
    # Otherwise fall back to direct calculation (backward compatibility)
    if provenance_info and "duration_seconds" in provenance_info:
        runtime = provenance_info["duration_seconds"]
    else:
        runtime = end_time - start_time

    if mutated_units is None and mutation_info and "mutated_units" in mutation_info:
        mutated_units = mutation_info["mutated_units"]

    haplotype_stats = generate_haplotype_stats(
        simulation_results, config, mutated_units=mutated_units
    )
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
        "provenance": provenance_info if provenance_info is not None else {},
    }
    return report


def write_statistics_report(report: dict[str, Any], output_path: str) -> None:
    """
    Write the simulation statistics report to a JSON file.

    Args:
        report (Dict[str, Any]): Simulation statistics report.
        output_path (str): File path for the JSON report.
    """
    try:
        with Path(output_path).open("w") as fh:
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
    # Dummy simulation results: two haplotypes.
    dummy_results = [
        HaplotypeResult(
            sequence="TTTTAAACCCGGGAAAA",
            chain=[RepeatUnit("1"), RepeatUnit("2"), RepeatUnit("9")],
        ),
        HaplotypeResult(
            sequence="TTTTAAACCCGGGAAAA",
            chain=[RepeatUnit("1"), RepeatUnit("2"), RepeatUnit("9", mutated=True)],
        ),
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
