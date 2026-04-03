"""Read source tracking for simulated reads.

This module provides data models and logic for annotating simulated reads
with their source haplotype, genomic position, repeat unit overlap, mutation
status, and SNP overlap. Annotations are written as gzip-compressed TSV
manifests suitable for benchmarking variant callers, aligners, and
VNTR-aware pipelines.

Key Components:
    build_coordinate_map: Build a repeat coordinate map from a haplotype chain
    ReadSourceTracker: Annotate reads and write manifests/coordinate maps

Data Models:
    RepeatRegion: A single repeat unit with coordinates and mutation status
    SNPPosition: A SNP with its position and optional repeat unit assignment
    RepeatCoordinateMap: Full coordinate map for one haplotype
    ReadOrigin: Raw read origin (haplotype, position, strand)
    AnnotatedRead: Read annotated with VNTR overlap, repeat units, mutations, SNPs
"""

from __future__ import annotations

import gzip
import json
import logging
from collections.abc import Iterable, Iterator
from dataclasses import dataclass

from ..type_defs import RepeatUnit

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


@dataclass
class RepeatRegion:
    """A single repeat unit within a haplotype VNTR region."""

    index: int
    """1-based repeat unit index within the chain."""

    repeat_type: str
    """Repeat symbol (e.g. '1', '2', 'X', 'Xm')."""

    start: int
    """0-based start coordinate (inclusive) within the haplotype sequence."""

    end: int
    """0-based end coordinate (exclusive) within the haplotype sequence."""

    is_mutated: bool
    """Whether this repeat unit carries a mutation."""

    mutation_name: str | None
    """Name of the applied mutation, or None if not mutated."""

    @property
    def length(self) -> int:
        """Length of this repeat region in bases."""
        return self.end - self.start


@dataclass
class SNPPosition:
    """A single-nucleotide polymorphism within a haplotype sequence."""

    position: int
    """0-based position within the haplotype sequence."""

    ref_base: str
    """Reference base at this position."""

    alt_base: str
    """Alternate (applied) base at this position."""

    repeat_index: int | None
    """1-based repeat unit index this SNP falls within, or None if outside VNTR."""


@dataclass
class RepeatCoordinateMap:
    """Complete coordinate map for one haplotype."""

    haplotype: int
    """1-based haplotype number."""

    regions: list[RepeatRegion]
    """Ordered list of repeat regions."""

    vntr_start: int
    """0-based start of the VNTR region (first repeat unit start)."""

    vntr_end: int
    """0-based end of the VNTR region (last repeat unit end, exclusive)."""

    snp_positions: list[SNPPosition]
    """SNPs mapped to this haplotype."""


@dataclass
class ReadOrigin:
    """Raw read origin extracted from simulator output."""

    read_id: str
    """Unique identifier for this read."""

    haplotype: int
    """1-based haplotype number the read was simulated from."""

    ref_start: int
    """0-based start position on the haplotype reference (inclusive)."""

    ref_end: int
    """0-based end position on the haplotype reference (exclusive)."""

    strand: str
    """Strand orientation: '+' or '-'."""


@dataclass
class AnnotatedRead:
    """A read annotated with source tracking information."""

    read_id: str
    haplotype: int
    ref_start: int
    ref_end: int
    strand: str
    overlaps_vntr: bool
    repeat_units: str
    overlaps_mutation: bool
    mutation_name: str
    overlaps_snp: bool


# ---------------------------------------------------------------------------
# Coordinate map builder
# ---------------------------------------------------------------------------


def build_coordinate_map(
    haplotype: int,
    chain: list[RepeatUnit],
    repeats_dict: dict[str, str],
    left_const_len: int,
    mutation_positions: list[tuple[int, int]],
    mutation_name: str | None,
    snp_info: list[dict[str, object]] | None,
) -> RepeatCoordinateMap:
    """Build a coordinate map from a haplotype repeat chain.

    Args:
        haplotype: 1-based haplotype number.
        chain: List of RepeatUnit objects.
        repeats_dict: Mapping of repeat symbol to DNA sequence.
        left_const_len: Length of the left constant region in bases.
        mutation_positions: List of (haplotype, repeat_index) tuples, both 1-based.
        mutation_name: Name of the mutation applied, or None.
        snp_info: List of SNP dicts with 'position', 'ref_base', 'alt_base' keys,
            or None if no SNPs.

    Returns:
        RepeatCoordinateMap for this haplotype.
    """
    regions: list[RepeatRegion] = []
    current_pos = left_const_len

    # Build set of mutated repeat indices for this haplotype (1-based)
    mutated_indices: set[int] = {rep_idx for hap, rep_idx in mutation_positions if hap == haplotype}

    for i, unit in enumerate(chain):
        repeat_index = i + 1  # 1-based

        if unit.symbol not in repeats_dict:
            raise ValueError(
                f"Unknown repeat symbol '{unit.symbol}' (from '{unit}') for haplotype "
                f"{haplotype} at index {repeat_index} in build_coordinate_map. "
                f"Available symbols: {sorted(repeats_dict.keys())}"
            )
        seq = repeats_dict[unit.symbol]
        region_len = len(seq)

        is_mutated = repeat_index in mutated_indices or unit.mutated
        region = RepeatRegion(
            index=repeat_index,
            repeat_type=str(unit),
            start=current_pos,
            end=current_pos + region_len,
            is_mutated=is_mutated,
            mutation_name=mutation_name if is_mutated else None,
        )
        regions.append(region)
        current_pos += region_len

    # VNTR boundaries
    if regions:
        vntr_start = regions[0].start
        vntr_end = regions[-1].end
    else:
        vntr_start = left_const_len
        vntr_end = left_const_len

    # Map SNP positions to repeat indices
    snp_positions: list[SNPPosition] = []
    if snp_info:
        for snp in snp_info:
            pos = int(snp.get("position", 0))  # type: ignore[call-overload]
            ref_base = str(snp.get("ref_base", ""))
            alt_base = str(snp.get("alt_base", ""))

            # Determine which repeat region this SNP falls within
            repeat_idx: int | None = None
            for region in regions:
                if region.start <= pos < region.end:
                    repeat_idx = region.index
                    break

            snp_positions.append(
                SNPPosition(
                    position=pos,
                    ref_base=ref_base,
                    alt_base=alt_base,
                    repeat_index=repeat_idx,
                )
            )

    return RepeatCoordinateMap(
        haplotype=haplotype,
        regions=regions,
        vntr_start=vntr_start,
        vntr_end=vntr_end,
        snp_positions=snp_positions,
    )


# ---------------------------------------------------------------------------
# ReadSourceTracker
# ---------------------------------------------------------------------------


class ReadSourceTracker:
    """Annotate simulated reads with source haplotype and VNTR overlap.

    Builds coordinate maps from repeat chains and provides methods to
    annotate reads, write manifests, and write coordinate map files.
    """

    def __init__(
        self,
        repeat_chains: dict[int, list[RepeatUnit]],
        repeats_dict: dict[str, str],
        left_const_len: int,
        mutation_positions: list[tuple[int, int]] | None = None,
        mutation_name: str | None = None,
        snp_info: dict[int, list[dict[str, object]]] | None = None,
    ) -> None:
        """Initialise the tracker and build coordinate maps.

        Args:
            repeat_chains: Mapping of 1-based haplotype number to RepeatUnit chain.
            repeats_dict: Mapping of repeat symbol to DNA sequence.
            left_const_len: Length of the left constant region in bases.
            mutation_positions: List of (haplotype, repeat_index) tuples, both
                1-based. Defaults to empty list.
            mutation_name: Name of the applied mutation, or None.
            snp_info: Mapping of 0-based haplotype index to list of SNP dicts.
                Keys are converted to 1-based internally.
        """
        if mutation_positions is None:
            mutation_positions = []
        if snp_info is None:
            snp_info = {}

        # Convert snp_info from 0-based haplotype index to 1-based
        snp_info_1based: dict[int, list[dict[str, object]]] = {
            hap_idx + 1: snps for hap_idx, snps in snp_info.items()
        }

        self.coordinate_maps: dict[int, RepeatCoordinateMap] = {}
        for hap_num, chain in repeat_chains.items():
            hap_snps = snp_info_1based.get(hap_num)
            self.coordinate_maps[hap_num] = build_coordinate_map(
                haplotype=hap_num,
                chain=chain,
                repeats_dict=repeats_dict,
                left_const_len=left_const_len,
                mutation_positions=mutation_positions,
                mutation_name=mutation_name,
                snp_info=hap_snps,
            )

    @classmethod
    def from_simulation_results(
        cls,
        results: list,
        config: dict,
        mutation_positions: list | None = None,
        mutation_name: str | None = None,
        applied_snp_info: list | None = None,
        reference_assembly: str | None = None,
    ) -> ReadSourceTracker:
        """Build a tracker from simulation results and config.

        Extracts repeat chains, constants, and SNP info from simulation
        outputs so callers don't need to reshape the data themselves.

        Args:
            results: List of HaplotypeResult objects.
            config: Configuration dict with 'repeats', 'constants', 'reference_assembly'.
            mutation_positions: List of MutationTarget objects (or None).
            mutation_name: Mutation name string (e.g. "dupC"), or None.
            applied_snp_info: List of SNP info lists, indexed by haplotype (0-based).
            reference_assembly: Override for reference assembly (defaults to config value).

        Returns:
            Configured ReadSourceTracker instance.
        """
        # Build repeat chains dict (1-based haplotype keys)
        repeat_chains: dict[int, list[RepeatUnit]] = {}
        for i, hr in enumerate(results):
            repeat_chains[i + 1] = hr.chain

        # Get left constant length
        ref_assembly = reference_assembly or config.get("reference_assembly", "hg38")
        left_const = config.get("constants", {}).get(ref_assembly, {}).get("left", "")

        # Get repeats dict
        repeats_dict = config.get("repeats", {})

        # Convert MutationTarget objects to tuples
        mut_positions: list[tuple[int, int]] = []
        mut_name: str | None = None
        if mutation_positions:
            mut_positions = [(mt.haplotype_index, mt.repeat_index) for mt in mutation_positions]
            if mutation_name and mutation_name != "normal":
                # Extract non-"normal" mutation name from comma-separated string
                for part in mutation_name.split(","):
                    if part != "normal":
                        mut_name = part
                        break

        # Build SNP info dict (0-based haplotype indices)
        snp_info_dict: dict[int, list[dict[str, object]]] = {}
        if applied_snp_info:
            for i, snp_list in enumerate(applied_snp_info):
                if snp_list:
                    snp_info_dict[i] = snp_list

        return cls(
            repeat_chains=repeat_chains,
            repeats_dict=repeats_dict,
            left_const_len=len(left_const),
            mutation_positions=mut_positions,
            mutation_name=mut_name,
            snp_info=snp_info_dict if snp_info_dict else None,
        )

    def annotate_reads(self, origins: Iterable[ReadOrigin]) -> Iterator[AnnotatedRead]:
        """Annotate reads with VNTR overlap, repeat units, and mutation status.

        Args:
            origins: Iterable of ReadOrigin objects to annotate.

        Yields:
            AnnotatedRead for each origin.
        """
        for origin in origins:
            coord_map = self.coordinate_maps.get(origin.haplotype)
            if coord_map is None:
                yield AnnotatedRead(
                    read_id=origin.read_id,
                    haplotype=origin.haplotype,
                    ref_start=origin.ref_start,
                    ref_end=origin.ref_end,
                    strand=origin.strand,
                    overlaps_vntr=False,
                    repeat_units="",
                    overlaps_mutation=False,
                    mutation_name="",
                    overlaps_snp=False,
                )
                continue

            # Check VNTR overlap: read [ref_start, ref_end) vs VNTR [vntr_start, vntr_end)
            overlaps_vntr = (
                origin.ref_start < coord_map.vntr_end and origin.ref_end > coord_map.vntr_start
            )

            # Find overlapping repeat regions
            overlapping_indices: list[int] = []
            overlaps_mutation = False
            mutation_name_found = ""

            for region in coord_map.regions:
                # Interval overlap check
                if origin.ref_start < region.end and origin.ref_end > region.start:
                    overlapping_indices.append(region.index)
                    if region.is_mutated:
                        overlaps_mutation = True
                        if region.mutation_name:
                            mutation_name_found = region.mutation_name

            repeat_units = ",".join(str(idx) for idx in overlapping_indices)

            # Check SNP overlap
            overlaps_snp = any(
                origin.ref_start <= snp.position < origin.ref_end for snp in coord_map.snp_positions
            )

            yield AnnotatedRead(
                read_id=origin.read_id,
                haplotype=origin.haplotype,
                ref_start=origin.ref_start,
                ref_end=origin.ref_end,
                strand=origin.strand,
                overlaps_vntr=overlaps_vntr,
                repeat_units=repeat_units,
                overlaps_mutation=overlaps_mutation,
                mutation_name=mutation_name_found,
                overlaps_snp=overlaps_snp,
            )

    def write_manifest(self, annotated_reads: Iterable[AnnotatedRead], path: str) -> None:
        """Write a gzip-compressed TSV manifest of annotated reads.

        Args:
            annotated_reads: Iterable of AnnotatedRead objects.
            path: Output file path (will be gzip-compressed).
        """
        header = (
            "read_id\thaplotype\tref_start\tref_end\tstrand\t"
            "overlaps_vntr\trepeat_units\toverlaps_mutation\t"
            "mutation_name\toverlaps_snp\n"
        )

        with gzip.open(path, "wt", encoding="utf-8") as f:
            f.write(header)
            for read in annotated_reads:
                repeat_units = read.repeat_units if read.repeat_units else "."
                f.write(
                    f"{read.read_id}\t"
                    f"{read.haplotype}\t"
                    f"{read.ref_start}\t"
                    f"{read.ref_end}\t"
                    f"{read.strand}\t"
                    f"{str(read.overlaps_vntr).lower()}\t"
                    f"{repeat_units}\t"
                    f"{str(read.overlaps_mutation).lower()}\t"
                    f"{read.mutation_name if read.mutation_name else '.'}\t"
                    f"{str(read.overlaps_snp).lower()}\n"
                )

    def write_coordinate_map(self, path: str) -> None:
        """Write a plain TSV coordinate map of repeat regions.

        Args:
            path: Output file path.
        """
        header = (
            "haplotype\tindex\trepeat_type\tstart\tend\t"
            "is_mutated\tmutation_name\tsnp_count\tsnp_positions\n"
        )

        with open(path, "w", encoding="utf-8") as f:
            f.write(header)
            for hap_num in sorted(self.coordinate_maps):
                coord_map = self.coordinate_maps[hap_num]

                # Pre-index SNPs by region for efficient lookup
                snps_by_region: dict[int, list[int]] = {}
                for snp in coord_map.snp_positions:
                    if snp.repeat_index is not None:
                        snps_by_region.setdefault(snp.repeat_index, []).append(snp.position)

                for region in coord_map.regions:
                    mut_name = region.mutation_name if region.is_mutated else "."
                    region_snps = snps_by_region.get(region.index, [])
                    snp_count = len(region_snps)
                    snp_pos_str = ",".join(str(p) for p in region_snps) if region_snps else "."
                    f.write(
                        f"{hap_num}\t"
                        f"{region.index}\t"
                        f"{region.repeat_type}\t"
                        f"{region.start}\t"
                        f"{region.end}\t"
                        f"{str(region.is_mutated).lower()}\t"
                        f"{mut_name}\t"
                        f"{snp_count}\t"
                        f"{snp_pos_str}\n"
                    )

    @classmethod
    def from_companion_files(cls, stats_path: str) -> ReadSourceTracker | None:
        """Reconstruct a ReadSourceTracker from a simulation_stats.json file.

        Extracts repeat chains, repeats dictionary, mutation details, SNP info,
        and left constant length from the stats JSON to build coordinate maps.

        Args:
            stats_path: Path to a simulation_stats.json file.

        Returns:
            A new ReadSourceTracker, or None if required fields are missing.
        """
        try:
            with open(stats_path, encoding="utf-8") as f:
                stats = json.load(f)
        except (OSError, json.JSONDecodeError) as exc:
            logger.warning("Could not load stats file %s: %s", stats_path, exc)
            return None

        # Extract repeat chains: stats["haplotypes"]["haplotype_N"]["repeat_chain"]
        haplotypes_section = stats.get("haplotypes")
        if not haplotypes_section or not isinstance(haplotypes_section, dict):
            logger.warning("Missing or invalid 'haplotypes' section in %s", stats_path)
            return None

        repeat_chains: dict[int, list[RepeatUnit]] = {}
        for key, hap_data in haplotypes_section.items():
            # Parse "haplotype_N" -> N
            if not key.startswith("haplotype_"):
                continue
            try:
                hap_num = int(key.split("_", 1)[1])
            except (ValueError, IndexError):
                continue
            chain_str = hap_data.get("repeat_chain")
            if not chain_str or not isinstance(chain_str, str):
                logger.warning("Missing repeat_chain for %s in %s", key, stats_path)
                return None
            repeat_chains[hap_num] = [RepeatUnit.from_str(s) for s in chain_str.split("-")]

        if not repeat_chains:
            logger.warning("No haplotype chains found in %s", stats_path)
            return None

        # Extract repeats dictionary: stats["config"]["repeats"]
        config_section = stats.get("config", {})
        repeats_dict = config_section.get("repeats")
        if not repeats_dict or not isinstance(repeats_dict, dict):
            logger.warning("Missing or invalid 'config.repeats' in %s", stats_path)
            return None

        # Extract left constant length: stats["config"]["constants"][assembly]["left"]
        constants = config_section.get("constants", {})
        assembly = stats.get("reference_assembly", "hg38")
        assembly_constants = constants.get(assembly, {})
        left_const = assembly_constants.get("left", "")
        left_const_len = len(left_const) if isinstance(left_const, str) else 0

        # Extract mutation details
        mutation_details = stats.get("mutation_details", {})
        mutation_name = mutation_details.get("mutation_name")
        raw_targets = mutation_details.get("targets", [])
        mutation_positions: list[tuple[int, int]] = []
        for target in raw_targets:
            if isinstance(target, (list, tuple)) and len(target) >= 2:
                mutation_positions.append((int(target[0]), int(target[1])))

        # Extract SNP info: stats["snp_info"] keyed by 1-based haplotype string
        # Convert to 0-based haplotype index for the constructor
        raw_snp_info = stats.get("snp_info", {})
        snp_info: dict[int, list[dict[str, object]]] = {}
        for hap_key, snps in raw_snp_info.items():
            try:
                hap_1based = int(hap_key)
                # Constructor expects 0-based keys
                snp_info[hap_1based - 1] = snps
            except (ValueError, TypeError):
                continue

        return cls(
            repeat_chains=repeat_chains,
            repeats_dict=repeats_dict,
            left_const_len=left_const_len,
            mutation_positions=mutation_positions,
            mutation_name=mutation_name,
            snp_info=snp_info if snp_info else None,
        )
