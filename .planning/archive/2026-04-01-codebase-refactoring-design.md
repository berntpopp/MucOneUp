# MucOneUp Codebase Refactoring Design

**Date:** 2026-04-01
**Status:** COMPLETE. All phases (1, 2A, 2B, 3A, 3B, 4A, 4B) merged. v0.36.0 released 2026-04-02.
**Based on:** `.planning/codebase-review-report.md` (overall rating: 5/10)
**Research sources:** Click Complex Applications docs, Python subprocess best practices, OpenStack shell-avoidance guide, Simon Willison CLI design principles

## Context

MucOneUp is a Python CLI tool for simulating MUC1 VNTR diploid references. A codebase review identified 13 findings across architecture, modularity, type safety, and operational robustness. The main issues are:

1. A 1586-line CLI god module with a fake argparse namespace bridging Click to the backend.
2. Weakly typed domain objects relying on raw `dict[str, Any]` and string conventions (`"m"` suffix, positional tuples).
3. Inconsistent subprocess handling — 44+ direct `subprocess` calls scattered across wrappers, bypassing the existing `run_command()` helper.
4. An import-time `rfc8785` dependency that breaks test collection and `--help` in minimal environments.
5. Config normalization that is not consistently applied — some CLI paths bypass `load_config()` and pass raw JSON downstream, causing flat-vs-nested constant key mismatches.
6. CLI options (`--out-dir`, `--out-base`) that don't reach the backend — the read-simulation API has no output location parameter.

## Constraints

- **No backward-compatibility requirements.** Pre-1.0, primary-user project. All internal and CLI APIs may break.
- **CI must stay green per phase.** GitHub Actions: ruff, mypy (continue-on-error), pytest across Python 3.10-3.12, 30% coverage threshold.
- **Strict phase isolation.** Each phase is self-contained. Phase N merged and CI-green before Phase N+1 starts.
- **One seed per CLI run** is the reproducibility model. No need for per-process parallel RNG isolation.
- **`rfc8785` stays as a dependency** but must be lazy-loaded. Config fingerprinting is metadata-only (attestation, not control) but valuable for reproducibility verification.

---

## Phase 1: Stabilize Basic Operability

**Goal:** Make `import muc_one_up.cli`, `muconeup --help`, `mypy`, and `pytest` collection work in all environments. Enforce a single config-loading boundary so all downstream code sees the same normalized shape.

### 1.1 Lazy-load `rfc8785` in `config_fingerprint.py`

Move the `import rfc8785` from module level into `compute_config_fingerprint()`. The module is only called from `collect_provenance_metadata()` at runtime, never at import time. This breaks the eager import chain:

```
CLI entry -> cli module -> analysis -> provenance -> config_fingerprint -> rfc8785 (FAILS)
```

### 1.2 Guard `config_fingerprint` import in `provenance.py`

Use a try/import pattern so provenance collection degrades gracefully if `rfc8785` is missing, returning an error sentinel fingerprint rather than crashing.

### 1.3 Enforce a single config-loading boundary

**This is the highest-value change in Phase 1.** Several CLI paths bypass `load_config()` and read raw JSON directly:

- `click_main.py:505` — Illumina reads command loads config but skips normalization
- `click_main.py:655` — ONT reads command, same issue
- `click_main.py:877` — PacBio reads command, same issue
- `click_main.py:1046` — ORF analysis reads raw config

`load_config()` in `config.py:476` normalizes flat legacy constants into the nested assembly-aware shape (`config["constants"]["hg38"]["left"]` vs `config["constants"]["left"]`). When CLI paths bypass this, downstream code sees inconsistent shapes. This is why ORF code still uses wrong flat-key access (`click_main.py:1049`, `analysis.py:60`).

**Fix:** Audit every CLI path that reads a config file. Replace all raw `json.load()` or direct dict access with calls to `load_config()`. Add a test that imports each command module and verifies it calls `load_config()` rather than parsing JSON directly. This is a prerequisite for Phase 3's `SimulationConfig` dataclass — if source dicts have inconsistent shapes, the dataclass constructor will produce wrong results.

### 1.4 Add smoke tests

- `test_cli_import.py`: Verify `from muc_one_up.cli import main` succeeds.
- `test_cli_help.py`: Verify `muconeup --help` exits 0 with expected output.
- Verify `mypy muc_one_up/` passes (fix missing stubs as needed).

### 1.5 Verify CI stays green

**Scope:** ~5-6 files modified, ~2 test files added.

**Targets:**
- `muc_one_up/config_fingerprint.py` (lazy import)
- `muc_one_up/provenance.py` (guarded import)
- `muc_one_up/cli/click_main.py` (replace raw JSON reads with `load_config()` calls)
- `muc_one_up/cli/__init__.py` (if import chain needs adjustment)

**Behavioral verification:**
- `muconeup --help` exits 0 in a virtualenv without `rfc8785`
- `pytest --collect-only` succeeds in a minimal environment
- `mypy muc_one_up/` passes
- All existing tests pass

---

## Phase 2: Fix the CLI Boundary

**Goal:** Split the god module, replace the fake argparse namespace, fix the read-simulation output contract, consolidate ORF analysis, and centralize error handling. The priority is making CLI options actually control backend behavior, not just producing a nicer file layout.

### 2.1 Fix the read-simulation output contract

**This is the highest-value change in Phase 2.** The `--out-dir` and `--out-base` CLI options currently don't control where read simulation output lands. The root cause is below Click:

- The CLI computes `actual_out_base` but never passes it downstream (`click_main.py:539`, `click_main.py:691`, `click_main.py:931`).
- The read-simulation API (`read_simulation.py:95`) accepts only `config` and `input_fa` — no output location parameter.
- Each pipeline derives output names from the input file path: `pipeline.py:138` (Illumina), `ont_pipeline.py:148` (ONT), `pacbio_pipeline.py:239` (PacBio).

**Fix:** Introduce an `OutputConfig` dataclass:

```python
@dataclass
class OutputConfig:
    out_dir: Path
    out_base: str
```

Change `simulate_reads_pipeline()` and the three pipeline entry points to accept `OutputConfig`. Each pipeline uses `output_config.out_dir / output_config.out_base` as the stem for all generated files instead of deriving names from `input_fa`. The CLI constructs `OutputConfig` from `--out-dir` and `--out-base` and passes it through.

**Behavioral verification:** Add integration tests that invoke each read command with `--out-dir /tmp/test_output --out-base my_reads` and assert that all generated files (FASTQ, BAM, index) land in the specified directory with the specified base name.

### 2.2 Consolidate ORF/toxic-protein analysis

ORF analysis logic exists in both the CLI god module (`click_main.py:1019`) and the CLI helper (`analysis.py:22`). These implementations have already diverged — `click_main.py:1049` uses flat constant keys while `analysis.py:60` uses assembly-keyed constants (a consequence of the normalization gap fixed in Phase 1).

**Fix:** Choose `analysis.py:22` as the single implementation (it's already partially extracted). Remove the inline ORF logic from `click_main.py`. The `analyze orfs` command and the orchestrator both call into `analysis.run_orf_prediction()`. Verify with a test that both paths produce identical results for the same input.

### 2.3 Split `click_main.py` into per-command modules

Following Click's [Complex Applications](https://click.palletsprojects.com/en/stable/complex/) pattern:

```
muc_one_up/cli/
    __init__.py          # exports main()
    click_main.py        # ~50 lines: @click.group + group.add_command() registration
    _common.py           # CONTEXT_SETTINGS, configure_logging, require_config, shared decorators
    commands/
        simulate.py      # simulate command
        reads.py         # @reads group + illumina, ont, pacbio subcommands
        analyze.py       # @analyze group + orfs, stats, vntr-stats, snapshot-validate
    options.py           # typed option dataclasses
    error_handling.py    # centralized exception-to-exit-code mapping
    config.py            # (existing) setup_configuration, process_mutation_config
    orchestration.py     # (existing) run_single_simulation_iteration
    analysis.py          # (existing, now sole ORF implementation)
    outputs.py           # (existing, extended with _generate_output_base)
```

Commands register via `cli.add_command()` in `click_main.py`. Shared config is passed through `@click.pass_context` with a typed dataclass on `ctx.obj`.

### 2.4 Replace `_make_args_namespace()` with typed dataclasses

The fake `argparse.Namespace` at `click_main.py:1410` is replaced with:

- `SimulationOptions` — seed, iterations, output paths, mutation config, mode flags.
- `ReadSimOptions` — coverage, threads, seed, `OutputConfig`, simulator-specific params.
- `AnalysisOptions` — ORF/stats/validation parameters.

Defined in `cli/options.py`. Each command constructs the appropriate dataclass from its Click parameters. Backend functions accept these typed objects instead of `Namespace`.

### 2.5 Centralize error handling

Create `cli/error_handling.py` with a decorator that maps domain exceptions (`FileOperationError`, `ReadSimulationError`, `ExternalToolError`) to consistent Click exit codes and user-facing messages. Remove scattered try/except blocks in individual commands that collapse errors into generic exits.

### 2.6 Move `_generate_output_base()` to `cli/outputs.py`

**Scope:** ~12 files modified/created, 1 large file broken up.

**Targets:**
- `muc_one_up/cli/click_main.py` (broken up into ~50-line root)
- `muc_one_up/cli/_common.py` (new, shared utilities)
- `muc_one_up/cli/commands/` (new directory with 3 command modules)
- `muc_one_up/cli/options.py` (new, typed dataclasses)
- `muc_one_up/cli/error_handling.py` (new)
- `muc_one_up/cli/outputs.py` (extended)
- `muc_one_up/cli/analysis.py` (sole ORF implementation)
- `muc_one_up/read_simulation.py` (API changed to accept `OutputConfig`)
- `muc_one_up/read_simulator/pipeline.py` (propagate `OutputConfig`)
- `muc_one_up/read_simulator/ont_pipeline.py` (propagate `OutputConfig`)
- `muc_one_up/read_simulator/pacbio_pipeline.py` (propagate `OutputConfig`)

**Behavioral verification:**
- All existing CLI tests pass with updated import paths.
- `--out-dir` and `--out-base` produce output in the specified location for all three read pipelines (new integration tests).
- `analyze orfs` and orchestrator ORF paths produce identical output for same input.
- Error exit codes are consistent and documented.

---

## Phase 3: Repair Domain Modeling

**Goal:** Replace implicit string conventions and raw dicts with explicit typed domain objects. Centralize assembly logic. Clean up RNG usage.

### 3.1 Typed domain models

New module `muc_one_up/models.py`:

- **`SimulationConfig`** dataclass: typed view of the config dict's simulation-relevant fields (repeats, probabilities, length_model, constants, mutations). Constructed from `ConfigDict` at the CLI boundary via a `SimulationConfig.from_config_dict()` classmethod. All domain functions accept this instead of raw dicts.
- **`HaplotypeResult`** dataclass: replaces `tuple[str, list[str]]` with named fields `sequence: str`, `chain: list[RepeatUnit]`.
- **`MutationTarget`** dataclass: replaces `tuple[int, int]` with `haplotype_index: int`, `repeat_index: int` (1-based, documented). Includes validation.
- **`MutationResult`** dataclass: holds the mutated haplotype, list of changes applied, and mutated units.

### 3.2 Replace the "m" suffix convention

Currently mutated repeats are marked by appending `"m"` to the symbol string, then stripped with `.rstrip("m")` for lookup. This is fragile (breaks if a repeat symbol legitimately ends in "m").

Introduce **`RepeatUnit`** dataclass: `symbol: str`, `mutated: bool`. The chain becomes `list[RepeatUnit]` instead of `list[str]`. Serialization to structure files preserves the `"Xm"` format for output compatibility, but internal logic uses the typed representation.

### 3.3 Centralize sequence/chain assembly

`assemble_haplotype_from_chain()` in `simulate.py` and `rebuild_haplotype_sequence()` in `mutate.py` do the same thing — concatenate constants + repeat sequences. Consolidate into a single `assemble_sequence(chain: list[RepeatUnit], constants: Constants) -> str` function in `muc_one_up/assembly.py`. Both `simulate.py` and `mutate.py` call this instead of maintaining their own assembly logic. This becomes the single source of truth for chain -> sequence conversion.

### 3.4 Pass RNG instances instead of using global `random`

Functions in `simulate.py`, `mutate.py`, `distribution.py`, `probabilities.py` currently call `random.choices()`, `random.choice()`, etc. on the global RNG.

Change signatures to accept `rng: random.Random` parameter. CLI entry point creates `rng = random.Random(seed)` and threads it through. This doesn't change behavior for single-seed CLI runs but makes the code explicit and testable.

### 3.5 Update `type_defs.py`

- Keep as canonical location for type aliases.
- Replace raw dict aliases with references to new dataclasses where appropriate.
- Remove `MutationTargets = list[tuple[int, int]]` in favor of `list[MutationTarget]`.

**Scope:** ~1 new module, ~6-8 modules modified.

**Targets:**
- `muc_one_up/models.py` (new)
- `muc_one_up/assembly.py` (new)
- `muc_one_up/simulate.py`
- `muc_one_up/mutate.py`
- `muc_one_up/distribution.py`
- `muc_one_up/probabilities.py`
- `muc_one_up/type_defs.py`
- `muc_one_up/cli/orchestration.py`
- `muc_one_up/cli/mutations.py`

**Behavioral verification:**
- All simulation outputs are byte-identical given the same seed (regression test).
- `mypy --strict muc_one_up/models.py muc_one_up/assembly.py` passes with no errors.
- No `dict[str, Any]` remains in domain function signatures (grep check).
- No `.rstrip("m")` or `+ "m"` in domain code (grep check).

---

## Phase 4: Consolidate Read-Simulation Infrastructure

**Goal:** Unify assembly resolution, replace all direct subprocess usage with a richer execution abstraction, eliminate wrapper duplication, and fix temp-file lifetime issues.

### 4.1 Unify assembly resolution

Introduce `AssemblyContext` dataclass (or extend `SimulationConfig` from Phase 3): holds resolved assembly name, left/right constants, VNTR region coordinates, human reference path. Constructed once at pipeline entry, passed through all steps. No function re-derives assembly from config.

Apply consistently across `pipeline.py` (Illumina), `ont_pipeline.py`, and `pacbio_pipeline.py`.

### 4.2 Build a richer execution abstraction

The existing `run_command()` in `common_utils.py:20` only returns an exit code and streams logs. This is why 44+ call sites bypass it — they need captured output, pipe chains, or redirected stderr. Before migrating callers, extend the abstraction.

**New API shape** (based on subprocess best practices research):

```python
@dataclass
class RunResult:
    returncode: int
    stdout: str | None   # None when streaming mode
    stderr: str | None
    command: str          # for error messages

def run_command(
    cmd: list[str],
    *,
    capture: bool = False,       # False = stream logs (current behavior)
    timeout: int | None = None,
    cwd: Path | None = None,
) -> RunResult: ...

def run_pipeline(
    cmds: list[list[str]],       # e.g., [["bwa", "mem", ...], ["samtools", "sort", ...]]
    *,
    capture: bool = True,
    timeout: int | None = None,
) -> RunResult: ...
```

Key design decisions:
- **No `shell=True`, ever.** Command construction stays separate from execution.
- **Pipe chains via `Popen` chaining.** Connect `stdout=PIPE` of process N to `stdin` of N+1. All processes in the same process group for cleanup.
- **Process-group cleanup.** Use `process_group=0` (Python 3.9+) instead of `preexec_fn=os.setsid` for fork safety. On timeout: SIGTERM, wait 5s, SIGKILL.
- **`RunResult` replaces bare `int` return.** Callers that need output use `capture=True` and read `result.stdout`.
- **Deadlock prevention.** Use `communicate()` for captured mode. Keep threaded line readers for streaming mode.

### 4.3 Migrate all subprocess callers

Direct `subprocess.run`/`Popen` calls exist in:
- `bwa_wrapper.py`: 13 calls (all ad-hoc, several are pipe chains)
- `samtools_wrapper.py`: 27 calls (mix of captured output and fire-and-forget)
- `nanosim_wrapper.py`: 4 calls (alignment output redirection)
- `read_simulator/utils/samtools.py`: direct subprocess calls
- `read_simulator/utils/bed.py`: direct subprocess calls
- `read_simulator/utils/tool_version.py`: direct subprocess calls

Migrate all to `run_command()` or `run_pipeline()`. The migration is mechanical once the abstraction supports capture and pipe modes.

### 4.4 Remove duplicated alignment paths

Audit `minimap2_wrapper.py` and `nanosim_wrapper.py` for overlapping alignment logic. Consolidate into a shared alignment function that both call.

### 4.5 Fix temp-file lifetime in ONT diploid handling

`diploid_handler.py`'s `run_split_simulation()` already uses a correct copy-out pattern internally — the `TemporaryDirectory` context manager owns the temp dir, work happens inside, and the merged FASTQ is written to a caller-specified `output_fastq` path before cleanup. The actual bug is narrower:

When `keep_intermediate=False`, the `DiploidSimulationResult` still stores temp-dir paths in `hap1_fastq`, `hap2_fastq`, `hap1_reference`, `hap2_reference` (lines 310-311, 326-327). These paths are invalid after the `with` block exits. The caller in `ont_pipeline.py:329` may read stale paths.

**Fix:** When `keep_intermediate=False`, set `hap1_fastq`, `hap2_fastq`, `hap1_reference`, `hap2_reference` to `None` in the returned `DiploidSimulationResult`. Update `ont_pipeline.py` and tests to handle `None` fields. This is a minimal, safe fix that prevents accidental use of stale paths without requiring an API migration to context managers.

### 4.6 Consistent cleanup and observability

Once all wrappers use the extended runner, timeout behavior, error reporting, and log output become uniform. Add a wrapper-level summary log line on completion (tool name, elapsed time, exit code).

**Scope:** ~12-14 files modified, `common_utils.py` substantially extended.

**Targets:**
- `muc_one_up/read_simulator/utils/common_utils.py` (new `RunResult`, capture mode, `run_pipeline()`)
- `muc_one_up/read_simulator/pipeline.py`
- `muc_one_up/read_simulator/ont_pipeline.py`
- `muc_one_up/read_simulator/wrappers/bwa_wrapper.py`
- `muc_one_up/read_simulator/wrappers/samtools_wrapper.py`
- `muc_one_up/read_simulator/wrappers/nanosim_wrapper.py`
- `muc_one_up/read_simulator/wrappers/minimap2_wrapper.py`
- `muc_one_up/read_simulator/utils/diploid_handler.py`
- `muc_one_up/read_simulator/utils/samtools.py`
- `muc_one_up/read_simulator/utils/bed.py`
- `muc_one_up/read_simulator/utils/tool_version.py`

**Behavioral verification:**
- `grep -r "subprocess\." muc_one_up/` returns only `common_utils.py` (grep check).
- All existing wrapper tests pass with `RunResult` returns.
- ONT diploid tests updated: `result.hap1_fastq is None` when `keep_intermediate=False`.
- Pipe-chain commands (bwa | samtools) produce identical output to current ad-hoc implementation.

---

## Phase Dependencies

```
Phase 1 -> Phase 2 -> Phase 3 -> Phase 4
```

- Phase 1 unblocks Phase 2: tests can collect, config is consistently normalized, CI baseline is reliable.
- Phase 2 creates clean CLI boundaries and fixes the output contract that Phase 3 domain models slot into.
- Phase 3 introduces `SimulationConfig`/`AssemblyContext` that Phase 4 propagates through pipelines.

---

## Success Criteria

Each phase defines **behavioral checks** (not just structural goals) that serve as regression guards:

### Phase 1
- `muconeup --help` exits 0 in a virtualenv without `rfc8785` installed.
- `pytest --collect-only` succeeds in a minimal environment.
- `mypy muc_one_up/` passes.
- `grep -r "json.load" muc_one_up/cli/` returns zero hits outside `config.py` (all config reads go through normalized loader).
- All existing tests pass.

### Phase 2
- `click_main.py` is under 100 lines (registration only).
- No `argparse.Namespace` usage anywhere in the codebase.
- Integration tests: `--out-dir /tmp/X --out-base Y` produces all output files in `/tmp/X/Y.*` for Illumina, ONT, and PacBio read commands.
- Single ORF implementation: `grep -rn "run_orf_finder" muc_one_up/cli/` returns only `analysis.py`.
- Error exit codes: each domain exception type maps to a documented, distinct exit code.
- All existing tests pass (with updated imports).

### Phase 3
- `grep -rn "dict\[str, Any\]" muc_one_up/{simulate,mutate,distribution,probabilities}.py` returns zero hits.
- `grep -rn '\.rstrip("m")' muc_one_up/` returns zero hits.
- `grep -rn "random\.(choice|choices|seed|randint)" muc_one_up/{simulate,mutate,distribution,probabilities}.py` returns zero hits (all use injected `rng`).
- Simulation output is byte-identical given the same seed (regression test comparing before/after).
- `mypy --strict muc_one_up/models.py muc_one_up/assembly.py` passes.

### Phase 4
- `grep -r "subprocess\." muc_one_up/` returns only `common_utils.py`.
- `RunResult` is the return type for all execution functions (type check).
- ONT diploid: `result.hap1_fastq is None` when `keep_intermediate=False` (test assertion).
- Pipe-chain outputs are byte-identical to current ad-hoc implementation (regression test).
- Assembly is resolved exactly once per pipeline run (single construction of `AssemblyContext`).

### All Phases
- CI green (ruff, mypy, pytest across 3.10-3.12).
- Coverage >= 30% threshold.
- No new ruff warnings.

---

## Completion Summary

**All phases completed and merged as of 2026-04-02.**

| Phase | Version | PR | Description | Key Deliverables |
|-------|---------|-----|-------------|-----------------|
| 1 | v0.31.0 | #55 | Stabilize operability | Lazy rfc8785, config normalization, smoke tests |
| 2A | v0.31.0 | #55 | Output contract & ORF | Single ORF implementation, output contract fix |
| 2B | v0.32.0 | #57 | CLI split & typed options | click_main.py split, typed option dataclasses |
| 3A | v0.33.0 | #58 | RNG threading & assembly | Centralized assembly.py, explicit RNG instances |
| 3B | v0.34.0 | #59 | Typed domain models | RepeatUnit, HaplotypeResult, MutationTarget dataclasses |
| 4A | v0.35.0 | #60 | Execution abstraction | RunResult, run_command capture/stdout_path, run_pipeline |
| 4B | v0.36.0 | #61 | Assembly context & fixes | AssemblyContext, stale temp-file fix, alignment dedup |

### Success criteria verification (final state at v0.36.0)

**Phase 1:** All criteria met.
**Phase 2:** All criteria met. click_main.py is ~75 lines (registration only). No argparse.Namespace usage.
**Phase 3:** Partially met.
- `.rstrip("m")` eliminated from domain code.
- Global random replaced with explicit rng instances.
- RepeatUnit, HaplotypeResult, MutationTarget replace raw types.
- `SimulationConfig` dataclass deferred (would touch every function signature — diminishing returns for pre-1.0 project).
- `dict[str, Any]` remains in domain signatures via `ConfigDict` (deferred with SimulationConfig).
**Phase 4:** All criteria met.
- subprocess centralized to common_utils.py (+ 1 justified conda exception).
- AssemblyContext constructed once per pipeline run.
- DiploidSimulationResult returns None for stale temp paths.
- nanosim alignment reuses samtools_wrapper sort_and_index_bam.

### Deferred items
- **SimulationConfig dataclass** (spec 3.1): Typed wrapper for ConfigDict. Would replace `config: ConfigDict` with `config: SimulationConfig` in all domain functions. Deferred as diminishing returns — the config dict pattern is deeply embedded and the project is pre-1.0.
- **MutationResult dataclass** (spec 3.1): The current `apply_mutations` return type is already descriptive. Adding another dataclass would be overengineering.
