# Read Source Tracking Design Spec

**Issue:** [#53 — Read source tracking: ground-truth annotations for simulated reads](https://github.com/berntpopp/MucOneUp/issues/53)
**Date:** 2026-03-25
**Status:** Draft

## Problem

MucOneUp tracks rich mutation metadata throughout the simulation pipeline (FASTA headers, structure files, `simulation_stats.json`), but this metadata is entirely severed at the read simulation boundary. External simulators encode read origin information — NanoSim in read names, pbsim3 in MAF files, reseq/Wessim2 in fragment coordinates — but MucOneUp discards it all. This prevents users from using MucOneUp as a benchmarking platform for VNTR variant callers, aligners, and analysis pipelines.

## Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Scope | Full manifest (benchmarking + aligner evaluation) | Both use cases need the same position + annotation data |
| BAM tags (`XH`, `XO`, `XM`, `XR`) | Deferred to follow-up issue | Requires `pysam` dependency; manifest provides same data |
| Efficiency bias handling | Track only surviving reads | Manifest stays 1:1 with final output; simpler validation |
| Flag placement | Both `reads` and `simulate` subcommands | Orchestration path is the primary workflow |
| Illumina fragment mapping | Sidecar file + order-based correlation | Zero risk to existing pipeline; header modification could break reseq |
| Platform scope | All three (Illumina, ONT, PacBio) in first deliverable | Shared infrastructure benefits from validating across all platforms |
| Repeat coordinate map | Persisted as TSV | Useful for debugging, external tools, and independent verification |
| Architecture | Centralized `ReadSourceTracker` class | Single place for coordinate math and annotation; parsers stay simple |

## Architecture

### Approach: Centralized Manifest Builder

A single `ReadSourceTracker` class receives simulation metadata, builds coordinate maps, accepts raw read origins from platform-specific parsers, annotates them, and writes the manifest. Platform parsers are simple extractors that return a standardized `ReadOrigin` format.

```
Simulation Phase (orchestration.py)
    │
    ├── repeat_chains, mutation_positions, snp_info
    │
    ▼
ReadSourceTracker (constructed with simulation metadata)
    │
    ├── Builds RepeatCoordinateMap per haplotype
    ├── Persists coordinate map to TSV
    │
    ▼
Passed to read simulation pipeline
    │
    ├── Platform parser extracts list[ReadOrigin]
    ├── tracker.annotate_reads(origins) → list[AnnotatedRead]
    └── tracker.write_manifest(path) → compressed TSV
```

## Component Design

### 1. Repeat Coordinate Map

**Module:** `muc_one_up/read_simulator/source_tracking.py`

Translates repeat chain indices into genomic coordinate ranges on the simulated reference.

**Input:**
- Repeat chain per haplotype (e.g., `['1', '2', 'X', 'Xm', '7', '8', '9']`)
- Repeats dictionary from config (symbol to DNA sequence)
- Left constant sequence length (offset for VNTR start)
- Mutation positions and SNP positions from simulation

**Data model — `RepeatRegion`:**
- `index`: 1-based repeat unit index
- `repeat_type`: symbol (e.g., `X`)
- `start` / `end`: 0-based genomic coordinates on the haplotype reference
- `is_mutated`: bool
- `mutation_name`: str or None

**Data model — `RepeatCoordinateMap`:**
- `haplotype`: haplotype number (1 or 2)
- `regions`: list of `RepeatRegion`
- `vntr_start` / `vntr_end`: derived from first/last repeat boundaries
- `snp_positions`: SNP positions mapped to the repeat units they fall within

**Data model — `SNPPosition`** (new dataclass; no `SNPInfo` exists in the codebase — raw SNP data from `integrate_snps_unified()` is a plain dict with keys `position`, `ref_base`, `alt_base`):
- `position`: 0-based genomic position on the haplotype reference
- `ref_base`: reference base
- `alt_base`: alternate base
- `repeat_index`: 1-based repeat unit index this SNP falls within (or None if outside VNTR)

**Persisted as:** `{output_base}_repeat_coordinates.tsv`

| Column | Type | Description |
|--------|------|-------------|
| `haplotype` | int | Haplotype number (1 or 2) |
| `index` | int | 1-based repeat unit index |
| `repeat_type` | str | Repeat symbol (e.g., `X`) |
| `start` | int | 0-based start coordinate |
| `end` | int | 0-based end coordinate (exclusive) |
| `is_mutated` | bool | Whether this repeat was mutated |
| `mutation_name` | str | Mutation name if mutated, `.` otherwise |
| `snp_count` | int | Number of SNPs within this repeat unit |
| `snp_positions` | str | Comma-separated 0-based SNP positions, `.` if none |

One map per haplotype, since haplotypes have different chains and lengths.

### 2. ReadSourceTracker Class

**Module:** `muc_one_up/read_simulator/source_tracking.py`

**Construction** (in `orchestration.py`, after simulation completes):

Parameters:
- `repeat_chains: dict[int, list[str]]` — per-haplotype chains
- `repeats_dict: dict[str, str]` — symbol to DNA sequence
- `left_const_len: int` — offset to VNTR start
- `mutation_positions: list[tuple[int, int]]` — (haplotype, repeat_idx) pairs
- `mutation_name: str` — mutation name
- `snp_info: dict[int, list[dict]]` — per-haplotype SNP list (keys are 0-based haplotype indices, matching `apply_snps_to_sequences()` output; converted to 1-based internally. Each dict has keys: `position`, `ref_base`, `alt_base`)

**Responsibilities:**
1. Build `RepeatCoordinateMap` per haplotype at construction
2. Persist coordinate maps to TSV via `write_coordinate_map(path)`
3. Annotate reads via `annotate_reads(origins: Iterable[ReadOrigin]) -> Iterator[AnnotatedRead]` (iterator-based to support streaming for large read sets; V1 may buffer internally but the API supports future streaming)
4. For each read, intersect `[ref_start, ref_end]` against coordinate map to determine:
   - `overlaps_vntr`: does the read span any part of the VNTR region?
   - `repeat_units`: which repeat indices are overlapped?
   - `overlaps_mutation` / `mutation_name`: does it overlap a mutated repeat?
   - `overlaps_snp`: does it overlap an applied SNP position?
5. Write compressed manifest via `write_manifest(path)` producing `{output_base}_read_manifest.tsv.gz`

**Data types:**
- `ReadOrigin`: platform-parser output — `read_id`, `haplotype`, `ref_start`, `ref_end`, `strand`
- `AnnotatedRead`: `ReadOrigin` plus all annotation columns from manifest schema

**Threading:** Passed from `orchestration.py` → `run_read_simulation()` → platform pipeline as an optional parameter (None when `--track-read-source` is not set). Pipelines check for its presence; if None, skip origin extraction entirely. Zero overhead when disabled.

### 3. Platform-Specific Parsers

Each parser extracts `list[ReadOrigin]` from the platform's native output. The `ReadSourceTracker` handles all annotation logic.

#### ONT (NanoSim) Parser

**Source data:** NanoSim read names encode origin: `{haplotype}_{position}_{status}_{index}_{strand}_{head}_{middle}_{tail}`

**Approach:**
- After NanoSim generates reads, parse FASTQ read names with regex
- Extract `haplotype` (1 or 2), `ref_start` (position field), `strand` (F/R)
- For `ref_end`: use aligned read span from the minimap2 BAM (CIGAR-derived reference consumption) rather than raw read length, since NanoSim introduces insertions/deletions that make read length differ from reference span. If BAM is not yet available at parsing time, use read length as initial approximation and refine after alignment.
- For diploid split-sim: haplotype assignment comes from which split file the read originated from

**Complexity:** Low.

#### PacBio (pbsim3) Parser

**Source data:** MAF alignment files generated alongside reads.

**PacBio pipeline specifics:** Unlike ONT, PacBio does NOT use split-simulation. pbsim3 simulates against the combined diploid FASTA and produces separate output files per sequence (e.g., `sd_0001.bam`, `sd_0002.bam`) with corresponding MAF files per sequence. Haplotype assignment comes from which per-sequence output file the reads belong to (sequence index maps to haplotype).

**HiFi consensus identity mapping:** The PacBio pipeline includes a CCS (Circular Consensus Sequence) stage that converts CLR subreads into HiFi consensus reads. CCS creates new read identifiers that don't correspond 1:1 to CLR reads. The MAF files describe CLR read positions, not HiFi read positions. To bridge this gap:
1. Parse MAF files to get CLR read origin positions
2. Group CLR reads by ZMW (zero-mode waveguide) — CLR reads from the same ZMW produce one HiFi read
3. For each HiFi read, use the consensus reference position from its constituent CLR reads (all CLR subreads from the same ZMW map to the same reference region)
4. If ZMW-level tracing proves too complex, fall back to: after CCS + alignment, use the HiFi read's aligned position from the BAM (minimap2 alignment against the simulated reference) as `ref_start`/`ref_end`, since the reference IS the simulated sequence. This fallback is less pure (alignment-derived vs. ground-truth) but still accurate for VNTR overlap detection since the reference is known.

**Approach:**
- When tracking is enabled, preserve MAF files during cleanup (currently deleted at `pacbio_pipeline.py` lines 292-293)
- Parse MAF `s` (sequence) lines: extract reference name, start position, strand, alignment length
- Extract haplotype from per-sequence output file index (sequence 1 = haplotype 1, etc.)
- Build CLR-to-position map, then correlate with HiFi reads via ZMW grouping or alignment fallback
- Compute `ref_start`, `ref_end` from MAF alignment blocks (CLR path) or BAM alignment (fallback path)
- Delete MAF files after parsing

**Complexity:** Highest among long-read platforms due to CLR-to-HiFi mapping.

#### Illumina (reseq/Wessim2) Parser

**Source data:** Sidecar file written during fragment simulation.

**Approach:**
1. Modify `fragment_simulation.py` to write `{output_base}_fragment_origins.tsv` alongside fragment FASTA — columns: `fragment_index`, `chrom`, `fstart`, `fend`, `strand`
2. After `reseq seqToIllumina` produces reads, correlate by output order: fragment N maps to read pair N
3. Map `chrom` to haplotype number using FASTA sequence names

**Reseq order preservation — blocking prerequisite:**
- Write an integration test that runs a small fragment set through reseq and verifies output read count matches input fragment count and order is preserved
- Test with varying fragment counts (1, 10, 100+) to confirm consistency
- If order is NOT preserved: fallback to embedding a unique ID in the fragment header's unused comment field and grepping for it in reseq output
- This must be verified before the Illumina parser implementation begins

**Post efficiency-bias handling:** After VNTR efficiency bias downsampling, filter the manifest to only include surviving read IDs by matching against the final BAM/FASTQ.

**Complexity:** Highest.

### 4. Metadata Threading and CLI Integration

#### Flag Propagation

`--track-read-source` (boolean, default False) added to:
- `reads` subcommand (Click option)
- `simulate` subcommand (propagated through orchestration)

When enabled:
1. `orchestration.py` constructs `ReadSourceTracker` after simulation completes, using already-available metadata (repeat chains, mutation positions, SNPs)
2. Passes tracker to `run_read_simulation()` as new optional parameter
3. `run_read_simulation()` passes tracker to the platform pipeline
4. Pipeline calls the appropriate parser, feeds results to `tracker.annotate_reads()`, then `tracker.write_manifest()`

When disabled (default): tracker is None, pipelines skip all tracking logic.

#### Dual Mutation Mode

When `--mutation-name normal,dupC` is used, orchestration produces two separate FASTA files (normal and mutated) and runs `simulate_reads_pipeline` twice. Each invocation requires its own `ReadSourceTracker` constructed with the appropriate metadata:
- **Normal tracker:** repeat chains without mutation markers, no mutation positions
- **Mutated tracker:** repeat chains with mutation markers, mutation positions populated

Each tracker produces its own manifest and coordinate map, named with the corresponding output base (e.g., `{base}_normal_read_manifest.tsv.gz`, `{base}_dupC_read_manifest.tsv.gz`).

#### Standalone `reads` Invocation

When `reads` is called directly (not through `simulate`), simulation metadata is not in memory. Resolution:

- **Read metadata from companion files:** Look for `simulation_stats.json` and structure file alongside the input FASTA. Required fields from companion files:
  - From structure file (`*.structure.txt`): repeat chains per haplotype (parsed via `io.py` `parse_structure_file()`)
  - From `simulation_stats.json`: `mutation_details` (name, targets), `snp_info` (position, ref/alt per haplotype), `repeats_dict` (symbol to sequence), `left_constant_length`
  - If stats JSON schema changes across versions, validate presence of required keys and warn on missing fields rather than failing
- **If companion files are missing:** Warn and produce a partial manifest with positions only (no VNTR/mutation/SNP annotations).

### 5. Manifest Schema

**File:** `{output_base}_read_manifest.tsv.gz` (gzip-compressed TSV)

| Column | Type | Description |
|--------|------|-------------|
| `read_id` | str | Read name from FASTQ/BAM |
| `haplotype` | int | Source haplotype (1 or 2) |
| `ref_start` | int | 0-based start on simulated reference |
| `ref_end` | int | 0-based end on simulated reference |
| `strand` | str | `+` or `-` |
| `overlaps_vntr` | bool | Whether read overlaps the VNTR region |
| `repeat_units` | str | Comma-separated repeat unit indices (e.g., `12,13,14`) or `.` |
| `overlaps_mutation` | bool | Whether read overlaps a mutated repeat |
| `mutation_name` | str | Mutation name if overlapping, `.` otherwise |
| `overlaps_snp` | bool | Whether read overlaps an applied SNP |

**Boolean serialization:** `true` / `false` (lowercase, for interoperability with standard TSV parsers and pandas).

**Paired-end Illumina reads:** Each mate gets its own row. Both mates share the same fragment origin coordinates (`ref_start`/`ref_end` reflect the full fragment span for both mates). Read IDs use `/1` and `/2` suffixes (or as determined by reseq output format). This ensures correct VNTR overlap detection — a mate that only partially covers the fragment still inherits the fragment's full coordinate span for annotation purposes.

### 6. Output Files

When `--track-read-source` is enabled:

| File | Platforms | Description |
|------|-----------|-------------|
| `{base}_repeat_coordinates.tsv` | All | Repeat unit coordinate map per haplotype |
| `{base}_read_manifest.tsv.gz` | All | Compressed read-level annotations |
| `{base}_fragment_origins.tsv` | Illumina only | Intermediate; deleted only after successful manifest generation (preserved on failure for retry) |

## Testing Strategy

### Unit Tests

- **Repeat coordinate map builder:** Given known chains + repeats dict, verify coordinates are correct. Edge cases: single repeat, empty chain, mutated repeats, varying repeat lengths.
- **ReadSourceTracker annotation:** Given a coordinate map and synthetic `ReadOrigin` entries, verify overlap detection for VNTR, mutation, and SNP regions. Edge cases: read spanning multiple repeats, read partially overlapping VNTR boundary, read entirely outside VNTR.
- **ONT parser:** Synthetic NanoSim-style read names produce correct `ReadOrigin` extraction.
- **PacBio parser:** Synthetic MAF file produces correct `ReadOrigin` extraction.
- **Illumina sidecar writer:** Verify fragment origins TSV is written with correct coordinates during fragment simulation.
- **Manifest writer:** Verify gzip-compressed TSV output with correct schema, headers, and data types.

### Integration Tests

- **Illumina order preservation (blocking prerequisite):** Small fragment set through reseq, verify output count and order match input. Test with varying fragment counts. Document fallback approach if assumption fails.
- **End-to-end per platform:** Run a minimal simulation + read generation with `--track-read-source`, verify:
  - Manifest file exists and is valid gzipped TSV
  - Row count matches read count in output FASTQ/BAM
  - Coordinate map file exists with expected columns
  - Every `read_id` in manifest exists in output
  - Annotation columns are populated (not all empty)
- **Standalone `reads` with companion files:** Verify manifest generation works when calling `reads` directly with existing simulation outputs.
- **Standalone `reads` without companion files:** Verify warning is emitted and partial manifest is produced.
- **Flag disabled:** Verify no tracking files are produced when `--track-read-source` is omitted.

### Documentation Updates

- Update CLI help text for `--track-read-source` flag on both `reads` and `simulate` subcommands
- Add read source tracking section to user documentation covering:
  - Feature overview and use cases
  - Manifest TSV schema (column descriptions, types, example values)
  - Repeat coordinates TSV schema
  - Usage examples for each platform
  - How to use manifests for benchmarking workflows
- Update CLAUDE.md if new make targets or workflow steps are added

## Scope Exclusions

- **BAM tags** (`XH`, `XO`, `XM`, `XR`): Deferred to follow-up issue. Requires `pysam` dependency.
- **Filtered read tracking:** Only surviving reads (post efficiency bias) are included. Pre-filter tracking is a separate concern.
- **Config file changes:** No new config fields anticipated. All behavior controlled via CLI flag.

## Risks

| Risk | Impact | Mitigation |
|------|--------|------------|
| reseq does not preserve fragment order | Illumina parser breaks | Blocking integration test before implementation; fallback to header-embedded ID |
| NanoSim read name format changes across versions | ONT parser breaks | Pin expected format in parser; validate with regex; fail loudly on mismatch |
| pbsim3 MAF format varies between versions | PacBio parser breaks | Validate MAF structure during parsing; test with known pbsim3 output |
| Large read counts cause memory pressure in manifest generation | Slow or OOM | Stream annotations — process reads in batches rather than loading all into memory |
| Diploid split-sim haplotype assignment | Incorrect haplotype labels | Use split file origin (already tracked in diploid_handler.py) rather than parsing read names |
| PacBio HiFi CCS loses CLR read identity | Cannot map HiFi reads to MAF positions | ZMW-based grouping as primary approach; alignment-based fallback using minimap2 BAM positions against known simulated reference |
