# Remaining Work — Post v0.37.0

Date: 2026-04-03 (revised)
Baseline: v0.37.0 (boundary hardening complete)

## Work Groups (by theme, ordered by execution priority)

---

### Group A: Orchestration & Analysis Test Coverage

**Why now:** `cli/orchestration.py` and `cli/analysis.py` are the main integration
paths but lack direct tests. Existing coverage comes indirectly from CLI-level
tests (`test_click_cli.py`, `test_cli.py`) which exercise defaults and branching
but don't test the orchestration/analysis functions in isolation. Missing: dual
mutation mode paths, source tracker construction, error branches, statistics
generation.

| # | Task | File | What's missing |
|---|------|------|----------------|
| A1 | Add direct tests for `run_single_simulation_iteration()` | `cli/orchestration.py` | Dual mode, source tracker branches, error paths |
| A2 | Add tests for `run_orf_prediction()` both modes | `cli/analysis.py` | Dual mutation ORF output, toxic detection integration |
| A3 | Add tests for `run_read_simulation()` both modes | `cli/analysis.py` | Normal + mutated variant dispatch, error wrapping |
| A4 | Add tests for `write_simulation_statistics()` | `cli/analysis.py` | Provenance collection, dual mode stats |
| A5 | Add tests for `_run_toxic_detection()` helper | `cli/analysis.py` | Config extraction, parameter defaults |
| A6 | Add tests for lazy `__getattr__` import | `read_simulator/__init__.py` | Success + AttributeError paths |

**Estimated effort:** 2 sessions. A1-A3 first.

---

### Group B: Doc/Code Consistency Audit

**Why now:** Documentation references nonexistent flags and impossible line numbers
in a post-refactor codebase. This is more urgent than writing new guides.

| # | Task | Evidence | Impact |
|---|------|----------|--------|
| B1 | Remove `--output-stats` references from simulation guide | `docs/guides/simulation.md:82,535,766` — flag does not exist | HIGH |
| B2 | Fix stale `click_main.py` line references in docs | `docs/guides/snapshot-validation.md:671` refs line 1235 in a 79-line file; `docs/guides/toxic-protein-detection.md:564` refs line 790 | HIGH |
| B3 | Audit all `docs/guides/*.md` for references to removed/renamed CLI options | Full sweep needed | HIGH |
| B4 | Enable mkdocs-click CLI reference | Already in deps (`pyproject.toml:35`), page exists (`docs/reference/cli.md.disabled`), just needs `mkdocs.yml:168` uncommented | LOW effort |
| B5 | Add `Documentation` URL to pyproject.toml | Missing from `[project.urls]` | LOW |

**Estimated effort:** 1 session. B1-B3 are the critical items.

---

### Group C: Targeted CLI UX Fixes

**Why now:** Most options already show defaults where it matters. The blanket
"38 options" count was inflated — most missing cases are `default=None` or
`is_flag=True, default=False` where `[default: None]` adds noise, not clarity.

| # | Task | Scope | Impact |
|---|------|-------|--------|
| C1 | Add `show_default=True` to `--log-level` | `cli/click_main.py:30` — the obvious concrete miss | MEDIUM |
| C2 | Audit remaining options: only add `show_default` where default is a meaningful value (not None/False) | `simulate.py`, `reads.py`, `analyze.py` | MEDIUM |
| C3 | Clarify `--mutation-name` help (explain "normal,dupC" dual format) | `simulate.py:77` | MEDIUM |
| C4 | Add config requirement note to `reads` and `analyze` group help | `reads.py`, `analyze.py` group docstrings | LOW |
| C5 | Add constraint info to `--pass-num` (>=2), explain RQ in `--min-rq` | `reads.py` | LOW |

**Estimated effort:** < 1 session.

---

### Group D: Orchestration Cleanup

**Why now:** `orchestration.py:132-234` has 100 lines of source tracker data
reshaping. Depends on test coverage from Group A.

| # | Task | Scope | Impact |
|---|------|-------|--------|
| D1 | Add `ReadSourceTracker.from_simulation_results()` factory | `source_tracking.py` | MEDIUM |
| D2 | Move source tracker construction out of orchestration | `orchestration.py` | MEDIUM |

**Estimated effort:** 1 session. Blocked by A1 (need tests before refactoring).

---

### Group E: Incremental Config Typing

**Why now:** `ConfigDict = dict[str, Any]` is the largest type safety gap, but
the config is a validated-then-mutated dict — CLI commands patch nested sections
in place (`reads.py:193`), and `load_config_raw()` intentionally skips validation.
A monolithic dataclass + `from_dict()` bridge fights this pattern.

**Approach:** Section-level `TypedDict`s that are still plain dicts at runtime.
Start with the most-accessed sections, apply incrementally.

| # | Task | Scope | Impact |
|---|------|-------|--------|
| E1 | Define `LengthModelConfig(TypedDict)` — tightest section, used in `distribution.py` | `type_defs.py` | MEDIUM |
| E2 | Define `ReadSimulationConfig(TypedDict)` — covers `config["read_simulation"]` | `type_defs.py` | MEDIUM |
| E3 | Define `NanosimConfig(TypedDict)` and `PacbioConfig(TypedDict)` | `type_defs.py` | MEDIUM |
| E4 | Define `MutationConfig(TypedDict)` — covers `config["mutations"][name]` | `type_defs.py` | MEDIUM |
| E5 | Define `ConstantsConfig(TypedDict)` — covers `config["constants"][assembly]` | `type_defs.py` | MEDIUM |
| E6 | Migrate callers incrementally: annotate function params with section TypedDicts | `simulate.py`, `mutate.py`, etc. | HIGH |

**Estimated effort:** 2 sessions. E1-E3 first (most accessed), then E4-E6.

**Key benefit over dataclass:** No `from_dict()` bridge needed. Existing code
that does `config["read_simulation"]["coverage"] = 30` keeps working.

---

### Group F: Remaining Test & Cleanup

| # | Task | Scope | Impact |
|---|------|-------|--------|
| F1 | Error path tests for pipeline advanced features | `read_simulator/pipeline.py` (57%) | MEDIUM |
| F2 | Error path tests for subprocess runner | `read_simulator/utils/common_utils.py` (59%) | MEDIUM |
| F3 | Consider enum for simulator types | `read_simulation.py` | LOW |
| F4 | Assess samtools_wrapper.py split (1018 LOC) | `samtools_wrapper.py` | LOW |
| F5 | Remove `cast(Any, ...)` in `determine_simulation_mode()` | `cli/config.py:198` | LOW |

---

## Execution Order

```
Phase 1: A1-A3 (orchestration/analysis tests) — protect code before changes
Phase 2: B1-B4 (doc/code consistency audit)   — fix actively misleading docs
Phase 3: D1-D2 (orchestration cleanup)        — now safe with tests from Phase 1
Phase 4: C1-C3 (targeted CLI UX)              — quick, low risk
Phase 5: E1-E6 (incremental config typing)    — biggest type safety win
Phase 6: A4-A6 + F (remaining tests, cleanup) — mop-up
```
