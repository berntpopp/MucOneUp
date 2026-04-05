# Systematic Codebase Hardening — v0.39.0

**Date:** 2026-04-04
**Status:** DRAFT
**Baseline:** v0.38.0 (all quality gates green: mypy 0 errors, ruff clean, 1240 tests passing, 87% coverage)
**Release strategy:** Single release (v0.39.0) batching all four waves.

## Context

MucOneUp completed a major refactoring cycle (v0.31.0–v0.38.0) that resolved the CLI god module, introduced typed domain objects, centralized subprocess handling, and added config TypedDicts. The codebase is in good shape (87% coverage, clean linting, clean mypy) but has accumulated specific debt:

- CLI commands duplicate option declarations and config setup logic
- `samtools_wrapper.py` at 1018 LOC exceeds maintainability thresholds
- Section-level TypedDicts exist but aren't used in domain function signatures
- CI coverage threshold (70% in `.github/workflows/test.yml:104`) has room to tighten toward the actual 87% baseline
- Error path tests are thin in subprocess wrappers
- Three modules still use global randomness

This spec defines a four-wave hardening plan that addresses all remaining items from `.planning/remaining-work-v0.37.md` and `.planning/codebase-review-report-v2.md`.

## Constraints

- CI must stay green after every wave.
- No behavioral changes to the CLI — all refactors are internal.
- One version bump at the end (v0.39.0).
- All work on a single feature branch, merged as one PR.

---

## Wave 1: Quick Wins — Dev Loop & CI

**Goal:** Improve the development feedback loop and prevent coverage regressions.

### 1a. Raise CI coverage threshold: 70% → 75%

The CI workflow (`.github/workflows/test.yml:104`) already enforces `--cov-fail-under=70`. Raise to 75% — still conservative against the actual 87%, but a tighter safety net. This is a small, safe change.

### 1b. Enable branch coverage

Add `branch = true` under `[tool.coverage.run]` in `pyproject.toml`. Branch coverage catches missed conditional paths (e.g., untested `if` branches) that line coverage misses. Verify the overall branch coverage percentage after enabling and adjust the threshold if needed.

### 1c. Add mypy override for rfc8785

Add a `[[tool.mypy.overrides]]` block for `rfc8785` with `ignore_missing_imports = true`. This matches the existing pattern for `orfipy` and `jsonschema`.

### 1d. Clean up pytest markers

Remove `slow` and `requires_tools` from `pyproject.toml` marker declarations. Wave 4 will re-add them when tests actually use them.

### 1e. Fix stale simulate command comment

`muc_one_up/cli/commands/simulate.py:150`: change `# IMPORTANT: Disable pipeline options (simulate is PURE)` to `# simulate command generates haplotypes only — no read simulation or ORF output`.

### Verification

- `uv run mypy muc_one_up/` reports 0 errors.
- `uv run pytest` passes with branch coverage enabled.
- CI threshold set to 75%.

---

## Wave 2: Structural — DRY & Module Split

**Goal:** Reduce duplication in CLI commands and split the largest module.

### 2a. reads.py deduplication

**Problem:** The illumina, ont, and pacbio subcommands in `cli/commands/reads.py` repeat 6 identical Click option declarations and identical config initialization blocks.

**Solution:**

1. **Shared decorator stack in `cli/options.py`:** Create a `shared_read_options()` function that applies the 6 common Click options (`input_fastas` argument, `--out-dir`, `--out-base`, `--coverage`, `--seed`, `--track-read-source`). Each subcommand applies `@shared_read_options` instead of declaring them inline.

2. **Config setup helper:** Extract `_setup_read_config(config, simulator_type, coverage, seed, seed_key)` that handles:
   - Ensuring `config["read_simulation"]` exists
   - Setting simulator type
   - Applying coverage with default fallback
   - Setting seed and logging it
   
   Each subcommand calls this instead of duplicating the setup block.

3. **Result:** illumina ~30 LOC, ont ~35 LOC, pacbio ~80 LOC (down from 93, 92, 180).

**CLI help-text preservation:** The refactor must produce identical `--help` output for all three subcommands (option names, order, help strings, defaults). Decorator stacking order in Click determines help order — verify with snapshot comparison before/after.

**Files modified:**
- `muc_one_up/cli/options.py` — add `shared_read_options()` decorator
- `muc_one_up/cli/commands/reads.py` — refactor 3 subcommands to use shared decorator and helper

### 2b. samtools_wrapper.py split (1018 LOC → 3 modules)

**Problem:** Single module with 11 functions spanning reference extraction, coverage, downsampling, sort/index, format conversion, and merging.

**Solution:** Split into 3 focused modules grouped by responsibility:

| New Module | Functions | ~LOC |
|---|---|---|
| `samtools_core.py` | `extract_subset_reference()`, `sort_and_index_bam()`, `merge_bam_files()` | ~190 |
| `samtools_coverage.py` | `calculate_vntr_coverage()`, `calculate_target_coverage()`, `downsample_bam()`, `downsample_entire_bam()` | ~280 |
| `samtools_convert.py` | `convert_sam_to_bam()`, `convert_bam_to_fastq()`, `convert_bam_to_paired_fastq()`, `FastqConversionOptions`, `_count_fastq_reads()` | ~480 |

`samtools_wrapper.py` becomes a thin re-export barrel:
```python
from .samtools_core import extract_subset_reference, sort_and_index_bam, merge_bam_files
from .samtools_coverage import calculate_vntr_coverage, calculate_target_coverage, downsample_bam, downsample_entire_bam
from .samtools_convert import convert_sam_to_bam, convert_bam_to_fastq, convert_bam_to_paired_fastq, FastqConversionOptions
```

All existing imports (`from ...wrappers.samtools_wrapper import X`) keep working. New code imports from specific modules.

**Test migration:** Existing `test_samtools_wrapper.py` tests continue working through the re-exports. New tests target specific modules.

**Files created:**
- `muc_one_up/read_simulator/wrappers/samtools_core.py`
- `muc_one_up/read_simulator/wrappers/samtools_coverage.py`
- `muc_one_up/read_simulator/wrappers/samtools_convert.py`

**Files modified:**
- `muc_one_up/read_simulator/wrappers/samtools_wrapper.py` — gutted to re-export barrel

### Verification

- All 1240+ tests pass.
- `reads illumina --help`, `reads ont --help`, `reads pacbio --help` output identical to pre-refactor (option names, order, help strings, defaults). Capture before/after snapshots to verify.
- `grep -c "def " muc_one_up/read_simulator/wrappers/samtools_wrapper.py` returns 0 (only re-exports).
- No module exceeds 500 LOC.

---

## Wave 3: Type Hardening — Domain Signatures

**Goal:** Narrow `config: ConfigDict` parameters to specific section TypedDicts where the dependency boundary is clean, and use small composite typed inputs where a function genuinely crosses multiple config sections.

**Strategy:** Audit each domain function's actual config access. Where a function reads from one section, replace `config: ConfigDict` with the section TypedDict. Where a function crosses sections (e.g., `apply_mutations` reads `config["mutations"]`, `config["repeats"]`, and calls `assemble_sequence` which needs constants), either pass multiple section params or introduce a small composite TypedDict. Do not force artificial narrowing that would just push complexity to every call site.

### Already narrow (no changes needed)

These functions already accept section-level types — the migration from PR #67 already landed here:

| Module | Function | Current Param | Status |
|---|---|---|---|
| `distribution.py` | `sample_repeat_count()` | `length_model: LengthModelDict` | Already narrow |
| `probabilities.py` | `pick_next_repeat()` | `probabilities: ProbabilitiesDict` | Already narrow |

### Clean single-section migrations

These functions access one config section and can be narrowed directly:

| Module | Function | Current Param | New Param | Config Access |
|---|---|---|---|---|
| `simulate.py:91` | `assemble_haplotype_from_chain()` | `config: ConfigDict` | `constants: AssemblyConstants` | Only reads `config["constants"]` via `assemble_sequence()` |
| `read_simulator/pipeline.py:76` | `simulate_reads_pipeline()` | `config: dict[str, Any]` | `config: ReadSimulationConfig` | Reads `config["read_simulation"]` and `config["tools"]` — may need a `ReadPipelineConfig` composite |

### Multi-section functions (composite input needed)

These functions genuinely cross config sections:

| Module | Function | Sections Accessed | Approach |
|---|---|---|---|
| `simulate.py:170` | `simulate_diploid()` | `config["length_model"]`, `config["probabilities"]`, `config["repeats"]`, `config["constants"]` | Keep `config: ConfigDict` but add a `SimulationSections` TypedDict that declares the expected keys, or pass individual sections as params |
| `simulate.py:214` | `simulate_single_haplotype()` | `config["probabilities"]`, `config["repeats"]`, `config["constants"]` | Same approach as `simulate_diploid` |
| `mutate.py:88` | `apply_mutations()` | `config["mutations"]`, `config["repeats"]` + calls `assemble_sequence` (needs constants) | Pass `mutations_section`, `repeats_section`, `constants` as separate params |
| `mutate.py:60` | `validate_allowed_repeats()` | `config: ConfigDict` only for `config["repeats"]` | Narrow to `repeats: RepeatsDict` (keys are the valid symbols) |
| `mutate.py:216` | `rebuild_haplotype_sequence()` | `config: ConfigDict` only for `assemble_sequence` call | Narrow to `constants: AssemblyConstants` |
| `simulate.py:108` | `simulate_from_chains()` | `config["repeats"]`, `config["constants"]` | Narrow to `repeats` + `constants` params |

For multi-section functions, the preferred approach is: pass individual typed sections as parameters rather than a monolithic config. Where this would result in >3 section params on one function, define a lightweight composite TypedDict (e.g., `SimulationSections`) that groups just the needed keys.

### What stays ConfigDict

- `config.py` — infrastructure (load/validate/normalize)
- `cli/commands/*.py` — top-level config variable passed through Click context
- `cli/orchestration.py` — orchestrates across sections, unpacks at boundary
- `read_simulator/pipeline.py` — may stay `dict[str, Any]` if the pipeline config shape is too intertwined with tool paths; evaluate during implementation

### TypedDict additions

If existing TypedDicts from PR #67 are missing fields needed by migrated functions, add them. New composite TypedDicts are allowed only when a function genuinely needs 3+ sections and passing individual params would be unwieldy. Each composite must be documented with which sections it groups and why.

### Verification

- `uv run mypy muc_one_up/` passes with 0 errors.
- `config: ConfigDict` eliminated from: `simulate.py` (all functions except where composite is used), `mutate.py` (all functions), `distribution.py`, `probabilities.py` (already done).
- All tests pass without modification (TypedDicts are structurally compatible with plain dicts).
- No function accepts `ConfigDict` when it only reads one section.

---

## Wave 4: Coverage & Polish

**Goal:** Fill error path test gaps, fix remaining global randomness, operationalize pytest markers, raise CI threshold.

### 4a. Error path tests for subprocess wrappers

Add parametrized tests for error scenarios:

| Target Module | Test Scenarios |
|---|---|
| `common_utils.py` (60% coverage) | Timeout, non-zero exit, missing binary, `run_pipeline` mid-chain failure |
| `samtools_core.py` (post-split) | `sort_and_index_bam` subprocess failure, `merge_bam_files` empty input list |
| `samtools_convert.py` (post-split) | Paired FASTQ validation failure, truncated BAM output |
| `pipeline.py` | Missing reference config, tool-not-found error paths |

Pattern per function:
```python
@pytest.mark.parametrize("returncode,stderr,expected_exc", [
    (1, "error: file not found", ExternalToolError),
    (137, "", ExternalToolError),  # SIGKILL / timeout
])
def test_sort_and_index_bam_errors(mocker, returncode, stderr, expected_exc):
    mocker.patch("...run_command", return_value=RunResult(returncode, None, stderr, "cmd"))
    with pytest.raises(expected_exc):
        sort_and_index_bam(...)
```

### 4b. Fix remaining global randomness

Thread `rng: random.Random` through:
- `cli/mutations.py:85` — `select_random_mutation()` or similar
- `snp_integrator.py:249` — SNP placement logic
- `read_simulator/fragment_simulation.py:330` — fragment sampling

Follow the existing pattern from `simulate.py` and `mutate.py` where `rng` is already threaded.

### 4c. Operationalize pytest markers

**Slow tests:**
- Run `pytest --durations=20` to identify tests > 1s.
- Add `@pytest.mark.slow` to those tests.
- Re-add `slow` marker to `pyproject.toml`.

**External tool tests (`requires_tools`):**

The pipeline depends on multiple external binaries: `samtools`, `bwa`, `reseq`, `pblat`, `faToTwoBit`, and simulator-specific tools (`nanosim`, `pbsim3`, `ccs`). A fixture-only approach won't work unless every test explicitly requests the fixture; a marker-driven approach is more robust.

Implementation: add a `pytest_collection_modifyitems` hook in root `conftest.py` that inspects `@pytest.mark.requires_tools("samtools", "bwa")` and skips the test if any named binary is absent via `shutil.which()`:

```python
def pytest_collection_modifyitems(config, items):
    for item in items:
        for marker in item.iter_markers("requires_tools"):
            for tool in marker.args:
                if not shutil.which(tool):
                    item.add_marker(pytest.mark.skip(
                        reason=f"requires external tool: {tool}"
                    ))
                    break
```

- Apply `@pytest.mark.requires_tools("samtools")`, `@pytest.mark.requires_tools("bwa", "samtools")`, etc. to tests that invoke external binaries.
- Re-add `requires_tools` marker to `pyproject.toml` with description noting it accepts tool names as arguments.
- Verify: `pytest -m "not requires_tools"` runs the fast/pure-Python subset.

### 4d. Raise CI coverage threshold to 80%

After all error path tests land, raise `--cov-fail-under` in `.github/workflows/test.yml` from 75% (Wave 1) to 80%.

### Verification

- `uv run pytest` passes with branch coverage.
- Coverage >= 85% with branch coverage enabled.
- `grep -rn "random\.(choice|choices|seed|randint)" muc_one_up/{cli/mutations,snp_integrator,read_simulator/fragment_simulation}.py` returns 0 hits.
- `pytest -m "not slow and not requires_tools"` runs fast subset successfully.

---

## Execution Order

```
Wave 1 → Wave 2 → Wave 3 → Wave 4 → version bump → PR
```

Waves are sequential — each builds on the previous. All work lands on a single feature branch and ships as one PR with version bump to v0.39.0.

## Success Criteria (final state)

- mypy: 0 errors
- ruff: clean
- pytest: all pass, branch coverage enabled, >= 85%
- CI threshold: 80%
- No module > 500 LOC (samtools_wrapper split)
- No duplicated Click option blocks in reads.py
- No function accepts `ConfigDict` when it only reads one config section (simulate, mutate, distribution, probabilities)
- No global `random` usage outside config/CLI boundary
- pytest markers actively used and meaningful
