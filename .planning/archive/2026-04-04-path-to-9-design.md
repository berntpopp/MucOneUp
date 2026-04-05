# Path to ~9/10: Typed Configs, Lazy Dispatch, Orchestration Thinning

## Goal

Address the 3 remaining findings from codebase-review-report-v3.md (findings #3, #4, #6) to bring the codebase from ~8/10 to ~9/10. All changes are backward-compatible refactors with no user-facing behavior changes.

## Scope

Three **independent** improvements. Each produces a separately shippable commit (or PR) with its own verification gate. The implementation plan must keep them isolated — no commit should mix changes from different sections.

1. **Lazy backend dispatch** (smallest) — Import only the selected simulator backend
2. **Typed config views** (medium) — Harden existing TypedDicts and apply them at pipeline boundaries
3. **Orchestration thinning** (medium) — Extract source tracker setup from orchestration

---

## 1. Lazy Backend Dispatch

**Problem:** `_get_simulator_map()` in `read_simulation.py:65` imports all 4 pipeline modules eagerly. Only one is used per invocation. The amplicon pipeline pulls in Biopython, which may not be needed.

**Change:** Replace `_get_simulator_map()` with `_get_simulator(simulator_type: str)` that imports only the selected backend via if/elif dispatch. Same return type (`Callable[..., str]`), same caller interface in `simulate_reads()`.

**Files:**
- Modify: `muc_one_up/read_simulation.py` — replace `_get_simulator_map()` with `_get_simulator()`

**Testing:** Extend `tests/test_lazy_imports.py` with two assertions:
1. After `import muc_one_up.read_simulation`, `"muc_one_up.read_simulator.amplicon_pipeline"` must NOT be in `sys.modules`.
2. After calling `simulate_reads()` with `simulator="amplicon"` (mocked to avoid real tool execution), the amplicon module IS in `sys.modules`.

**Verification:** `python -m pytest tests/test_lazy_imports.py -v` passes. All existing tests pass.

---

## 2. Typed Config Views (TypedDicts)

**Problem:** Pipeline functions accept `config: dict[str, Any]` and access nested keys without type checking. Key mismatches are invisible to mypy (the ONT coverage bug was one such case).

**Change:** `muc_one_up/type_defs.py` already defines `ReadSimulationConfig`, `NanosimConfig`, and `PacbioConfig`. This task hardens and completes them:

### Fixes to existing TypedDicts

- **`NanosimConfig`**: Remove stale `min_len: int` key (line 186). The correct key is `min_read_length`, which is already present. Add missing keys: `other_options: str`, `correction_factor: float`, `enable_split_simulation: bool`, `enable_coverage_correction: bool`.
- **`PacbioConfig`**: Already correct. No changes needed.
- **`ReadSimulationConfig`**: Already has the core keys. Add `seed: int` if missing.

### New TypedDict

- **`AmpliconConfig`**: New TypedDict for `amplicon_params` section:
  ```python
  class AmpliconConfig(TypedDict, total=False):
      forward_primer: str
      reverse_primer: str
      expected_product_range: tuple[int, int] | None
      pcr_bias: dict[str, Any]  # nested config, leave inner structure untyped
  ```

### Annotation at boundaries

Pipeline functions annotate the config sections they extract. Since the root config is `dict[str, Any]`, extracting a section requires a cast:

```python
# In ont_pipeline.py
from ..type_defs import NanosimConfig
ns_params = cast(NanosimConfig, config.get("nanosim_params", {}))
coverage = ns_params.get("coverage")  # mypy knows: float | None
```

CLI commands annotate sections they construct:

```python
# In reads.py ont()
from ...type_defs import NanosimConfig
ns = cast(NanosimConfig, config.setdefault("nanosim_params", {}))
ns["min_read_length"] = min_read_length  # mypy checks key validity
```

**Files:**
- Modify: `muc_one_up/type_defs.py` — fix `NanosimConfig`, add `AmpliconConfig`
- Modify: `muc_one_up/read_simulator/ont_pipeline.py` — cast + annotate `ns_params`
- Modify: `muc_one_up/read_simulator/pacbio_pipeline.py` — cast + annotate `pb_params`
- Modify: `muc_one_up/read_simulator/amplicon_pipeline.py` — cast + annotate `amplicon_params`
- Modify: `muc_one_up/read_simulator/pipeline.py` — cast + annotate `rs_config` (Illumina)
- Modify: `muc_one_up/cli/commands/reads.py` — cast + annotate config sections in ont(), pacbio(), amplicon()
- Modify: `muc_one_up/read_simulator/utils/metadata_writer.py` — cast + annotate config section reads

**Verification:** `python -m mypy muc_one_up/` passes with 0 errors. No dedicated type-test harness — mypy is the verification gate. All existing tests pass unchanged (TypedDicts are structurally compatible with plain dicts at runtime).

---

## 3. Orchestration Thinning

**Problem:** `run_single_simulation_iteration()` in `orchestration.py` builds `ReadSourceTracker` objects, resolves assembly defaults, and writes coordinate map files (~65 lines of domain logic mixed with orchestration, lines 132-195).

**Change:** Extract that entire block into a helper function. The helper owns all three responsibilities (tracker creation, coordinate map file writing, and logging) because they form one cohesive unit.

```python
# muc_one_up/cli/source_tracking_setup.py

from __future__ import annotations

from muc_one_up.read_simulator.source_tracking import ReadSourceTracker
from muc_one_up.type_defs import HaplotypeResult

def build_source_trackers(
    *,
    config: dict[str, Any],
    results: list[HaplotypeResult],
    mutated_results: list[HaplotypeResult] | None,
    mutation_positions: dict[int, list[tuple[int, str]]] | None,
    mutation_name: str | None,
    applied_snp_info_normal: dict | None,
    applied_snp_info_mut: dict | None,
    reference_assembly: str | None,
    out_dir: str,
    out_base: str,
    sim_index: int,
    dual_mutation_mode: bool,
) -> tuple[ReadSourceTracker | None, ReadSourceTracker | None]:
    """Build read source trackers and write coordinate map files.

    Resolves the reference assembly (defaults to hg38), constructs
    ReadSourceTracker objects via from_simulation_results(), writes
    coordinate map TSV files, and logs paths.

    Returns:
        (normal_tracker, mutated_tracker). In non-dual mode,
        mutated_tracker is None.
    """
    ...
```

The orchestration function becomes:

```python
if getattr(args, "track_read_source", False) and getattr(args, "simulate_reads", None):
    source_tracker, source_tracker_mut = build_source_trackers(
        config=config,
        results=results,
        mutated_results=mutated_results,
        mutation_positions=mutation_positions,
        mutation_name=getattr(args, "mutation_name", None),
        applied_snp_info_normal=applied_snp_info_normal,
        applied_snp_info_mut=applied_snp_info_mut,
        reference_assembly=getattr(args, "reference_assembly", None),
        out_dir=out_dir,
        out_base=out_base,
        sim_index=sim_index,
        dual_mutation_mode=dual_mutation_mode,
    )
```

**Side effects owned by the helper:**
- Constructs `ReadSourceTracker` objects (1 or 2 depending on dual mode)
- Resolves `reference_assembly` (falls back to `config.get("reference_assembly", "hg38")`)
- Writes coordinate map TSV files to `out_dir/out_base.NNN[.normal|.mut].repeat_coordinates.tsv`
- Logs file paths via `logging.info`

**Files:**
- Create: `muc_one_up/cli/source_tracking_setup.py` — extracted helper
- Modify: `muc_one_up/cli/orchestration.py` — replace inline block (~65 lines) with helper call (~15 lines)

**Testing:** Existing orchestration tests pass unchanged. Add unit test for `build_source_trackers` verifying:
- Returns `(tracker, None)` in non-dual mode
- Returns `(tracker, mut_tracker)` in dual mode
- Writes coordinate map files to expected paths

**Verification:** `python -m pytest --tb=short -q` passes. `ruff check` passes.

---

## Out of Scope

- Replacing root `ConfigDict` with a typed object (too large, touches every module)
- Refactoring `load_config_raw()` or config validation
- Domain model changes beyond type annotations
- Any user-facing CLI changes

## Success Criteria

Per-task (each verified independently):

1. **Lazy dispatch:** `_get_simulator_map()` no longer exists; importing `read_simulation` does not load `amplicon_pipeline` into `sys.modules`.
2. **Typed configs:** `NanosimConfig.min_len` removed; `AmpliconConfig` added; all pipeline entry points use cast + TypedDict annotations; `python -m mypy muc_one_up/` reports 0 errors.
3. **Orchestration:** `run_single_simulation_iteration()` is ~50 lines shorter; `build_source_trackers()` exists with full test coverage for both modes.

Global: all existing tests pass, ruff clean, mypy green.
