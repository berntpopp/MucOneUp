"""Source tracker construction for simulation orchestration.

Extracted from orchestration.py to separate domain logic (tracker creation,
coordinate map writing) from orchestration sequencing.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from ..read_simulator.source_tracking import ReadSourceTracker
from ..type_defs import HaplotypeResult


def build_source_trackers(
    *,
    config: dict[str, Any],
    results: list[HaplotypeResult],
    mutated_results: list[HaplotypeResult] | None,
    mutation_positions: Any,
    mutation_name: str | None,
    applied_snp_info_normal: Any,
    applied_snp_info_mut: Any,
    reference_assembly: str | None,
    out_dir: str,
    out_base: str,
    sim_index: int,
    dual_mutation_mode: bool,
) -> tuple[ReadSourceTracker | None, ReadSourceTracker | None]:
    """Build read source trackers and write coordinate map files.

    Resolves the reference assembly (defaults to config value, then hg38),
    constructs ReadSourceTracker objects via from_simulation_results(),
    writes coordinate map TSV files, and logs paths.

    Returns:
        (normal_tracker, mutated_tracker). In non-dual mode,
        mutated_tracker is None.
    """
    ref_assembly = reference_assembly or config.get("reference_assembly", "hg38")

    if dual_mutation_mode:
        # Normal tracker: no mutations
        source_tracker = ReadSourceTracker.from_simulation_results(
            results=results,
            config=config,
            applied_snp_info=applied_snp_info_normal,
            reference_assembly=ref_assembly,
        )

        # Mutated tracker: with mutations, using mutated chains
        source_tracker_mut = ReadSourceTracker.from_simulation_results(
            results=mutated_results if mutated_results else results,
            config=config,
            mutation_positions=mutation_positions,
            mutation_name=mutation_name,
            applied_snp_info=applied_snp_info_mut,
            reference_assembly=ref_assembly,
        )

        # Write coordinate maps for both variants
        normal_coord_path = str(
            Path(out_dir) / f"{out_base}.{sim_index:03d}.normal.repeat_coordinates.tsv"
        )
        source_tracker.write_coordinate_map(normal_coord_path)
        logging.info("Normal repeat coordinate map written: %s", normal_coord_path)

        mut_coord_path = str(
            Path(out_dir) / f"{out_base}.{sim_index:03d}.mut.repeat_coordinates.tsv"
        )
        source_tracker_mut.write_coordinate_map(mut_coord_path)
        logging.info("Mutated repeat coordinate map written: %s", mut_coord_path)

        return source_tracker, source_tracker_mut
    else:
        source_tracker = ReadSourceTracker.from_simulation_results(
            results=results,
            config=config,
            mutation_positions=mutation_positions,
            mutation_name=mutation_name,
            applied_snp_info=applied_snp_info_normal,
            reference_assembly=ref_assembly,
        )

        coord_map_path = str(Path(out_dir) / f"{out_base}.{sim_index:03d}.repeat_coordinates.tsv")
        source_tracker.write_coordinate_map(coord_map_path)
        logging.info("Repeat coordinate map written: %s", coord_map_path)

        return source_tracker, None
