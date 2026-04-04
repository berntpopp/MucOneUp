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
- CI coverage threshold (30%) doesn't protect the actual 87% baseline
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

### 1a. Raise CI coverage threshold: 30% → 75%

Edit `pyproject.toml` pytest addopts or the CI workflow to set `--cov-fail-under=75`. This is conservative — current coverage is 87% — but leaves headroom for new modules that start with lower coverage.

### 1b. Enable branch coverage

Add `branch = true` under `[tool.coverage.run]` in `pyproject.toml`. Branch coverage catches missed conditional paths (e.g., untested `if` branches) that line coverage misses. Verify the overall branch coverage percentage after enabling and adjust the threshold if needed.

### 1c. Add mypy override for rfc8785

Add a `[[tool.mypy.overrides]]` block for `rfc8785` with `ignore_missing_imports = true`. This matches the existing pattern for `orfipy` and `jsonschema`.

### 1d. Clean up pytest markers

Remove `slow` and `requires_tools` from `pyproject.toml` marker declarations. Wave 4 will re-add them when tests actually use them.

### 1e. Fix stale simulate.py comment

Line 150: change `# IMPORTANT: Disable pipeline options (simulate is PURE)` to `# simulate command generates haplotypes only — no read simulation or ORF output`.

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
- `reads --help` output unchanged for all 3 subcommands.
- `grep -c "def " muc_one_up/read_simulator/wrappers/samtools_wrapper.py` returns 0 (only re-exports).
- No module exceeds 500 LOC.

---

## Wave 3: Type Hardening — Domain Signatures

**Goal:** Replace `config: ConfigDict` with specific section TypedDicts in domain function signatures.

**Strategy:** Domain functions receive narrow config sections, not the full dict. The CLI orchestration layer unpacks the full config into sections at the boundary.

### Target migrations

| Module | Function(s) | Current Param | New Param |
|---|---|---|---|
| `distribution.py` | `select_repeat_count()`, `generate_repeat_chain()` | `config: ConfigDict` | `length_config: LengthModelConfig` |
| `probabilities.py` | `compute_repeat_probabilities()` | `config: ConfigDict` | `length_config: LengthModelConfig` |
| `mutate.py` | `apply_mutations()`, `select_mutation_targets()` | `config: ConfigDict` | `mutation_config: MutationConfig` |
| `simulate.py` | `simulate_haplotype()` | `config: ConfigDict` | `constants: ConstantsConfig`, `length_config: LengthModelConfig` (unpack at call site) |
| `read_simulation.py` | `simulate_reads_pipeline()` | `config: ConfigDict` | `read_config: ReadSimulationConfig` |

### Caller updates

`cli/orchestration.py` and `cli/commands/*.py` unpack the full config:
```python
# Before
result = simulate_haplotype(config, ...)
# After
result = simulate_haplotype(config["constants"][assembly], ...)
```

### What stays ConfigDict

- `config.py` (load/validate/normalize) — infrastructure, not domain
- `cli/commands/*.py` top-level config variable — the full config dict is still passed through Click context
- `orchestration.py` config parameter — it orchestrates across sections

### TypedDict additions

If any existing TypedDicts from PR #67 are missing fields needed by the migrated functions, add them. No new top-level TypedDicts — use the existing `LengthModelConfig`, `ReadSimulationConfig`, `NanosimConfig`, `PacbioConfig`, `MutationConfig`, `ConstantsConfig`.

### Verification

- `uv run mypy muc_one_up/` passes with 0 errors.
- `grep -rn "config: ConfigDict" muc_one_up/{simulate,mutate,distribution,probabilities,read_simulation}.py` returns 0 hits.
- All tests pass without modification (TypedDicts are structurally compatible with plain dicts).

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

- Run `pytest --durations=20` to identify slow tests (> 1s).
- Add `@pytest.mark.slow` to those tests.
- Add a `requires_tools` fixture in `conftest.py`:
  ```python
  @pytest.fixture
  def requires_samtools():
      if not shutil.which("samtools"):
          pytest.skip("samtools not found")
  ```
- Apply `@pytest.mark.requires_tools` to tests that need external binaries.
- Re-add `slow` and `requires_tools` markers to `pyproject.toml`.

### 4d. Raise CI coverage threshold to 80%

After all error path tests land, raise `--cov-fail-under` from 75% (Wave 1) to 80%.

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
- No `config: ConfigDict` in domain function signatures (simulate, mutate, distribution, probabilities, read_simulation)
- No global `random` usage outside config/CLI boundary
- pytest markers actively used and meaningful
