# Pipeline Decomposition Design Spec

**Date:** 2026-04-04
**Issue:** #71 — Decompose Illumina pipeline (expanded to all 4 pipelines)
**Status:** Approved

## Goal

Decompose the 4 read simulation pipeline monoliths (Illumina 711 LOC, PacBio 495 LOC, ONT 382 LOC, Amplicon 335 LOC) into well-defined modules with shared utilities, and add true E2E tests that run all workflows from input FASTA to output BAM with real tools locally.

## Approach: Composition with Shared Utility Functions

Based on research into Python bioinformatics best practices (Snakemake, Biopython, scikit-bio patterns) and guidance from Effective Python and Hynek Schlawack:

- **No base class or protocol.** The 4 pipelines share structural similarity (setup → simulate → align → track → finalize) but not behavioral similarity — each simulation stage is deeply different.
- **Shared utility functions** for cross-cutting concerns (output resolution, metadata, cleanup).
- **Module-level separation** for the Illumina pipeline (largest, most complex).
- **Private helper extraction** for ONT/PacBio/Amplicon (smaller, don't warrant separate modules).

## Section 1: Shared Pipeline Utilities

**New module: `muc_one_up/read_simulator/pipeline_utils.py`**

| Function | Responsibility | Replaces duplication in |
|----------|---------------|------------------------|
| `resolve_pipeline_outputs(input_fa, rs_config, output_config, suffix) -> tuple[Path, str, Path]` | Output dir, base name, and path resolution | All 4 pipelines (~20 LOC each) |
| `resolve_human_reference(config, assembly_ctx) -> str` | Reference resolution with 3-path fallback + validation | Illumina, ONT |
| `create_pipeline_metadata(output_dir, base, config, start, end, platform, tools) -> None` | Timing capture + `write_metadata_file()` call | All 4 pipelines |
| `track_and_cleanup(intermediate_files, intermediate_bams, final_bam, keep_intermediates) -> None` | Safe cleanup that never deletes final outputs | All 4 pipelines |

**Not shared (intentionally):**

- Source tracking — each pipeline's parser is genuinely different (Illumina: fragment origins sidecar + BAM filtering; ONT: split-sim haplotype parsing; PacBio: MAF file parsing; Amplicon: not supported).
- Core simulation stages — no behavioral overlap.
- Assembly context construction — only Illumina uses it fully.

## Section 2: Illumina Pipeline Decomposition

**Current state:** Single 615-line `simulate_reads_pipeline()` function in `pipeline.py`.

**Proposed split:** `pipeline.py` becomes a thin orchestrator (~80-100 LOC) that calls into extracted stage modules.

### `muc_one_up/read_simulator/stages/fragment_preparation.py`

- Function: `prepare_fragments(tools, rs_config, input_fa, output_dir, base, assembly_ctx) -> FragmentResult`
- `tools`: dict of tool paths; `rs_config`: read_simulation subsection only (not full config)
- Stages 1-8: Replace Ns → systematic errors → 2bit → extract subset ref → pblat → simulate fragments → create reads → split reads
- Returns a `FragmentResult` dataclass with paths to R1/R2 FASTQ and list of intermediate files
- ~120 LOC, fully deterministic, no optional branches

### `muc_one_up/read_simulator/stages/alignment.py`

- Function: `align_and_refine(tools, rs_config, r1, r2, human_ref, output_dir, base, assembly_ctx) -> AlignmentResult`
- `tools`: dict of tool paths; `rs_config`: read_simulation subsection only
- Stages 9-11: bwa align → optional VNTR capture efficiency bias → optional downsampling
- Returns `AlignmentResult` with final BAM path, intermediate BAMs list
- Fixes the mutable `output_bam` problem: uses explicit `aligned_bam`, `vntr_biased_bam`, `downsampled_bam` variable names; assigns `final_bam` once at the end
- ~220 LOC; VNTR bias and downsampling branches stay here (tightly coupled to alignment output)

### `muc_one_up/read_simulator/stages/source_manifest.py`

- Function: `generate_read_manifest(sidecar_path, final_bam, input_fa, source_tracker) -> Path | None`
- Stage 12a: Parse fragment origins, match to surviving BAM reads, annotate haplotypes, write manifest TSV
- ~70 LOC, only called when `source_tracker is not None`

### Result Dataclasses (`stages/__init__.py`)

```python
@dataclass(frozen=True)
class FragmentResult:
    r1_fastq: str
    r2_fastq: str
    intermediate_files: list[str]

@dataclass(frozen=True)
class AlignmentResult:
    final_bam: str
    intermediate_bams: list[str]
```

### Orchestrator (`pipeline.py`)

Becomes:

```
setup & validate → prepare_fragments() → align_and_refine() → generate_read_manifest() → create_pipeline_metadata() → track_and_cleanup()
```

Each stage function receives only what it needs — no passing the entire config dict through. The orchestrator unpacks config and passes specific parameters.

**Public API preserved:** `simulate_reads_pipeline(config, input_fa, source_tracker, output_config) -> str` signature unchanged.

## Section 3: ONT, PacBio, and Amplicon Pipeline Decomposition

Lighter decomposition — extract shared utilities and break up the single large function, but not into separate files per stage.

### `ont_pipeline.py` (382 LOC → ~250 LOC)

- Extract output resolution and metadata to `pipeline_utils.py`.
- Extract diploid split-simulation decision logic into helper: `_resolve_simulation_mode(input_fa, config) -> tuple[bool, list[str]]` — returns whether to split-sim and the haplotype FASTA paths.
- Keep the rest inline — only 2 core stages, not worth separate modules.
- Move late imports to function top with clear feature guards.

### `pacbio_pipeline.py` (495 LOC → ~350 LOC)

- Extract output resolution, metadata, and cleanup to shared utils.
- Extract per-haplotype CLR+CCS loop into: `_simulate_haplotype_hifi(hap_fa, config, seed, hap_idx) -> list[str]` — returns list of HiFi BAM paths for one haplotype.
- Keep alignment and merge inline.

### `amplicon_pipeline.py` (335 LOC → ~220 LOC)

- Extract output resolution, metadata, and cleanup to shared utils.
- Amplicon-specific stages (primer extraction, PCR bias, template generation) stay inline.
- Extract per-haplotype template+CCS loop into: `_simulate_haplotype_amplicon(hap_fa, templates, config, seed, hap_idx) -> list[str]`.
- Pattern matches PacBio's `_simulate_haplotype_hifi` for consistency.

**Key principle:** These pipelines don't get a `stages/` directory. Private helper functions within the same file are the right level of abstraction.

## Section 4: CLI Dispatch Consolidation

**Current state:** `reads.py` already has `_run_batch_simulation()` and `_setup_read_config()`. The remaining duplication is PacBio/Amplicon `pacbio_params` setup.

**New helper:**

```python
def _apply_pacbio_params(
    config: dict,
    model_type: str | None,
    model_file: str | None,
    seed: int | None,
    threads: int | None = None,
    pass_num: int | None = None,
    min_passes: int | None = None,
    min_rq: float | None = None,
) -> None:
```

- Initializes `pacbio_params` section.
- Applies CLI overrides (only non-None values).
- Validates required params (`model_type`, `model_file`).
- Used by both `pacbio()` and `amplicon()` commands.

**PacBio command** shrinks from ~45 lines of config setup to ~5 lines. **Amplicon command** similarly shrinks, keeping its amplicon-specific `pcr_bias` and `amplicon_params` setup.

## Section 5: E2E Testing Strategy

### Approach

True end-to-end tests that simulate a VNTR repeat and produce BAM files from all 4 pipeline workflows, running locally with real bioinformatics tools.

### Test Fixture: Minimal Simulated VNTR Repeat

- Shared fixture generates a small diploid FASTA (~500bp per haplotype) using MucOneUp's own `simulate_from_chains()` with a short chain (e.g., `1-2-X-B-6-7-8-9`).
- Real MUC1 repeat — not a synthetic toy sequence — so we test realistic data flow.
- Fixture writes companion files (structure, simulation stats) needed for source tracking.

### 4 E2E Test Functions

| Test | Tools required | Validates |
|------|---------------|-----------|
| `test_illumina_e2e` | `reseq`, `bwa`, `samtools`, `pblat`, `faToTwoBit` | R1/R2 FASTQ + sorted BAM + BAM index exist; BAM has aligned reads; metadata JSON written |
| `test_ont_e2e` | `nanosim`, `minimap2`, `samtools` | BAM exists with long reads; split-sim produces balanced haplotype coverage |
| `test_pacbio_e2e` | `pbsim3`, `ccs`, `minimap2`, `samtools` | HiFi BAM exists; reads pass CCS quality threshold |
| `test_amplicon_e2e` | `pbsim3`, `ccs`, `minimap2`, `samtools` | Amplicon BAM exists; reads span expected amplicon region |

### Gating

- All tests marked `@pytest.mark.e2e` and `@pytest.mark.requires_tools(...)`.
- Skipped automatically if any required tool is missing (existing `conftest.py` infra handles this).
- Not run in CI — add `pyproject.toml` config: default pytest excludes `e2e` marker.

### Assertions Per Test

1. Output files exist (BAM, BAM index, metadata JSON).
2. BAM is valid (pysam can open it, has reads, has header).
3. Read count > 0 (simulation actually produced something).
4. Intermediate files cleaned up (unless `keep_intermediate_files` is set).

### Config Fixtures

- Each test builds a minimal config with real tool paths discovered via `shutil.which()`.
- Uses small coverage (e.g., 5x) and minimal parameters for speed.
- Amplicon test needs primer sequences in `amplicon_params` — use known MUC1 primer pair from existing config examples.

### Estimated Local Run Time

~2-5 minutes total (small inputs, low coverage).

## Section 6: File Layout

### New Files

```
muc_one_up/read_simulator/pipeline_utils.py              # Shared utilities
muc_one_up/read_simulator/stages/__init__.py              # FragmentResult, AlignmentResult dataclasses
muc_one_up/read_simulator/stages/fragment_preparation.py  # Illumina stages 1-8
muc_one_up/read_simulator/stages/alignment.py             # Illumina stages 9-11
muc_one_up/read_simulator/stages/source_manifest.py       # Illumina stage 12a
tests/e2e/__init__.py
tests/e2e/conftest.py                                     # Shared E2E fixtures
tests/e2e/test_pipeline_e2e.py                            # 4 E2E tests
```

### Modified Files

```
muc_one_up/read_simulator/pipeline.py         # Slim orchestrator (~80-100 LOC)
muc_one_up/read_simulator/ont_pipeline.py     # Use shared utils, extract helpers
muc_one_up/read_simulator/pacbio_pipeline.py  # Use shared utils, extract haplotype helper
muc_one_up/read_simulator/amplicon_pipeline.py  # Use shared utils, extract haplotype helper
muc_one_up/cli/commands/reads.py              # Extract _apply_pacbio_params
pyproject.toml                                # Add e2e marker config
```

### Not Touched

- Wrappers (`bwa_wrapper.py`, `samtools_wrapper.py`, etc.) — interfaces unchanged.
- `read_simulation.py` strategy dispatcher — unchanged, same function signatures.
- Existing unit/integration tests — should continue passing since public API doesn't change.
- `cli/analysis.py` — calls `simulate_reads_pipeline` which keeps same signature.

## Architectural Issues Fixed

1. **Mutable `output_bam` references** — Illumina's `alignment.py` uses explicit variable names (`aligned_bam`, `vntr_biased_bam`, `downsampled_bam`) and assigns `final_bam` once at the end.
2. **Scattered late imports** — moved to function top with feature guards in all pipelines.
3. **Triple-path reference resolution** — consolidated into `resolve_human_reference()` in `pipeline_utils.py`.
4. **CLI param duplication** — PacBio/Amplicon share `_apply_pacbio_params()`.

## Out of Scope

- Changing the strategy dispatcher in `read_simulation.py`.
- Modifying wrapper interfaces.
- Adding new simulation features.
- CI integration for E2E tests (requires bioinformatics tools not available in CI).
