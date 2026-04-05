# Codebase Review Report v2

Date: 2026-04-03
Repository: MucOneUp
Version: 0.36.0
Previous review: `.planning/codebase-review-report.md` (2026-04-01, score 5/10)

## Verification Commands Run

`ruff check muc_one_up/ tests/`

```text
All checks passed!
```

`mypy muc_one_up/`

```text
muc_one_up/read_simulator/utils/common_utils.py:184:16: error: Incompatible types in assignment ...
muc_one_up/read_simulator/utils/common_utils.py:229:9: error: "CompletedProcess[str]" has no attribute "wait"
muc_one_up/read_simulator/utils/common_utils.py:233:34: error: "CompletedProcess[str]" has no attribute "pid"
muc_one_up/read_simulator/utils/samtools.py:245:16: error: Item "None" of "str | None" has no attribute "strip"
muc_one_up/read_simulator/wrappers/samtools_wrapper.py:804:22: error: Incompatible types in assignment ...
Found 8 errors in 3 files (checked 75 source files)
```

`pytest --tb=short -q`

```text
collected 650 items / 34 errors
E   ModuleNotFoundError: No module named 'Bio'
E   ModuleNotFoundError: No module named 'orfipy_core'
!!!!!!!!!!!!!!!!!!! Interrupted: 34 errors during collection !!!!!!!!!!!!!!!!!!!
```

`find muc_one_up/ -name "*.py" | xargs wc -l | sort -n`

```text
     79 muc_one_up/cli/click_main.py
    264 muc_one_up/cli/orchestration.py
    397 muc_one_up/cli/analysis.py
    484 muc_one_up/cli/commands/analyze.py
    579 muc_one_up/cli/commands/reads.py
    711 muc_one_up/read_simulator/pipeline.py
   1018 muc_one_up/read_simulator/wrappers/samtools_wrapper.py
  17724 total
```

## Overall Rating

| Dimension | Score | Rationale |
| --- | ---: | --- |
| Architecture | 6/10 | The CLI root is now thin and command registration is cleaner, but import direction is still leaky: `muc_one_up/cli/click_main.py:63`, `muc_one_up/cli/analysis.py:13`, `muc_one_up/read_simulation.py:63`, `muc_one_up/read_simulator/__init__.py:10`. |
| Modularity | 6/10 | The refactor split the former god module, but large multi-responsibility modules remain, especially `muc_one_up/cli/commands/reads.py:71`, `muc_one_up/cli/orchestration.py:24`, `muc_one_up/read_simulator/pipeline.py:104`, `muc_one_up/read_simulator/wrappers/samtools_wrapper.py:657`. |
| DRY | 6/10 | ORF logic is better centralized, but the three read commands still duplicate batch iteration, output naming, and source-tracker reconstruction in `muc_one_up/cli/commands/reads.py:137`, `muc_one_up/cli/commands/reads.py:292`, `muc_one_up/cli/commands/reads.py:535`. |
| KISS | 6/10 | Several areas are now simpler, but there is still avoidable orchestration/config shaping in `muc_one_up/cli/orchestration.py:132` and compatibility layering in `muc_one_up/read_simulation.py:93` that adds surface area without fully hiding complexity. |
| SOLID | 6/10 | Command modules and `AssemblyContext` improve SRP/OCP, but the orchestration layer still constructs domain artifacts directly in `muc_one_up/cli/orchestration.py:140`, and config access remains infrastructure-shaped in `muc_one_up/simulate.py:90` and `muc_one_up/mutate.py:59`. |
| Type safety | 4/10 | The codebase has more annotations than before, but the typed boundary is not holding: `muc_one_up/read_simulator/utils/common_utils.py:22`, `muc_one_up/read_simulator/utils/common_utils.py:184`, `muc_one_up/read_simulator/utils/samtools.py:240`, `muc_one_up/read_simulator/wrappers/samtools_wrapper.py:804`, `muc_one_up/type_defs.py:114`. |
| Testability | 5/10 | There is broad test coverage by count, but import-time coupling prevents collection in a non-fully-provisioned environment: `muc_one_up/cli/__init__.py:10`, `muc_one_up/translate.py:52`, `muc_one_up/read_simulator/utils/reference_utils.py:13`. |
| Operational robustness | 6/10 | Subprocess handling is more centralized and temp-path handling appears improved, but the runner contract is internally inconsistent and still mixes shell and non-shell paths in `muc_one_up/read_simulator/utils/common_utils.py:57` and `muc_one_up/read_simulator/wrappers/samtools_wrapper.py:740`. |
| CLI UX | 6/10 | Help text and command grouping are materially better, but there is documentation drift and some contract softness: `muc_one_up/cli/click_main.py:47`, `muc_one_up/cli/commands/simulate.py:133`, `muc_one_up/config.py:574`, `muc_one_up/cli/commands/analyze.py:83`. |
| Tooling/release hygiene | 5/10 | `ruff` is clean and packaging metadata is reasonable, but the project currently misses a clean mypy run and the requested test command does not collect successfully in this environment: `pyproject.toml:31`, `pyproject.toml:120`, `muc_one_up/read_simulator/utils/common_utils.py:184`, `muc_one_up/translate.py:52`. |

## Executive Summary

The refactor materially improved the shape of the codebase: the old CLI god module is gone, the command tree is clearer, and the read-simulation stack has more explicit infrastructure such as `AssemblyContext`, `OutputConfig`, and centralized process helpers. This is no longer a 5/10 codebase.

The main weaknesses shifted rather than disappeared. The most important regression is import-time coupling: public CLI and read-simulation imports now eagerly pull in Biopython and `orfipy_core`, which breaks test collection and basic module usability in a partially provisioned environment. Type safety also lags the refactor, with mypy failures concentrated in the new subprocess abstraction.

Compared with the previous review, most structural issues are genuinely improved, but the remaining debt is now concentrated in boundary hygiene: dependency loading, typed subprocess contracts, and CLI/config seams.

## Comparison with Previous Review

| Old Finding | Current Status | Notes |
| --- | --- | --- |
| 1. Import-time CLI coupling breaks basic operability | REGRESSED | `rfc8785` is no longer the problem, but eager imports still break CLI/test collection through `Bio` and `orfipy_core` (`muc_one_up/cli/__init__.py:10`, `muc_one_up/cli/analysis.py:13`, `muc_one_up/translate.py:52`, `muc_one_up/read_simulator/utils/reference_utils.py:13`). |
| 2. `click_main.py` is a god module | RESOLVED | `muc_one_up/cli/click_main.py:1` is now a small registration layer; command logic moved into dedicated modules. |
| 3. CLI contracts do not consistently match backend behavior | PARTIALLY RESOLVED | Output handling looks cleaner, but utility commands still bypass validated config and mutate raw dicts inline (`muc_one_up/config.py:574`, `muc_one_up/cli/commands/reads.py:103`, `muc_one_up/cli/commands/reads.py:482`). |
| 4. ORF/toxic-protein workflows are duplicated and drifting | RESOLVED | Shared ORF logic is centralized through `muc_one_up/cli/analysis.py:22` and standalone helpers are invoked from commands. |
| 5. CLI orchestration is carrying domain logic | STILL OPEN | `muc_one_up/cli/orchestration.py:132` still constructs source-tracking state and reshapes mutation/SNP metadata directly. |
| 6. The core domain is coupled to raw config dictionaries | STILL OPEN | `ConfigDict = dict[str, Any]` remains the dominant contract (`muc_one_up/type_defs.py:114`, `muc_one_up/simulate.py:90`, `muc_one_up/mutate.py:59`). |
| 7. Domain state is encoded implicitly in strings and tuples | PARTIALLY RESOLVED | `RepeatUnit` and `MutationTarget` improved matters, but compatibility handling remains in `muc_one_up/type_defs.py:41` and `muc_one_up/validation.py:146`. |
| 8. Mutation handling has weak source-of-truth boundaries | RESOLVED | The typed chain model and rebuild helpers are a real improvement; no comparable drift is obvious in current review. |
| 9. Read simulation resolves assembly inconsistently | RESOLVED | `AssemblyContext` is now constructed centrally in `muc_one_up/read_simulator/pipeline.py:139` and `muc_one_up/read_simulator/ont_pipeline.py:122`. |
| 10. External tool boundaries are inconsistent and duplicated | PARTIALLY RESOLVED | Centralization improved the common path, but `samtools_wrapper` still keeps a custom shell pipeline path beside `run_pipeline` (`muc_one_up/read_simulator/wrappers/samtools_wrapper.py:740`, `muc_one_up/read_simulator/wrappers/samtools_wrapper.py:804`). |
| 11. Temporary-file lifetime management is fragile | RESOLVED | I did not find the pre-refactor stale-temp-path pattern in the currently reviewed split-flow code. |
| 12. Randomness is global and not injectable | PARTIALLY RESOLVED | The core simulator improved, but global randomness is still used in adjacent workflows (`muc_one_up/cli/mutations.py:85`, `muc_one_up/snp_integrator.py:249`, `muc_one_up/read_simulator/fragment_simulation.py:330`). |
| 13. Error handling is inconsistent and too broad | PARTIALLY RESOLVED | Central CLI mapping exists, but several broad wrappers remain (`muc_one_up/cli/analysis.py:53`, `muc_one_up/cli/analysis.py:199`, `muc_one_up/read_simulator/utils/common_utils.py:194`). |

## New Findings

1. **High**: Import-time dependency coupling has regressed and now breaks CLI/test collection in partially provisioned environments.
   Evidence: `muc_one_up/cli/__init__.py:10`, `muc_one_up/cli/click_main.py:64`, `muc_one_up/cli/analysis.py:13`, `muc_one_up/translate.py:52`, `muc_one_up/read_simulation.py:63`, `muc_one_up/read_simulator/__init__.py:10`, `muc_one_up/read_simulator/utils/reference_utils.py:13`, `muc_one_up/analysis/snapshot_validator.py:23`.
   Impact: importing `muc_one_up.cli`, `muc_one_up.read_simulator`, or tests that touch those packages now eagerly requires `Bio` and `orfipy_core`, which is exactly why `pytest --tb=short -q` aborts during collection in this environment. This is a real architecture/testability regression even though the specific missing package changed from the previous review.

2. **High**: The new subprocess abstraction is not type-sound, and the breakage leaks into wrapper code.
   Evidence: `muc_one_up/read_simulator/utils/common_utils.py:22`, `muc_one_up/read_simulator/utils/common_utils.py:105`, `muc_one_up/read_simulator/utils/common_utils.py:145`, `muc_one_up/read_simulator/utils/common_utils.py:184`, `muc_one_up/read_simulator/utils/common_utils.py:229`, `muc_one_up/read_simulator/utils/samtools.py:240`, `muc_one_up/read_simulator/utils/samtools.py:303`, `muc_one_up/read_simulator/wrappers/samtools_wrapper.py:804`.
   Impact: `RunResult` is the intended unified contract, but the implementation still mixes `CompletedProcess`, `Popen`, and nullable text outputs in ways mypy can no longer reconcile. This lowers trust in the runner at exactly the layer meant to standardize tool execution.

3. **Medium**: CLI utility commands rely on raw, partially validated config dictionaries and mutate them in-place, weakening the new command/backend boundary.
   Evidence: `muc_one_up/config.py:574`, `muc_one_up/cli/commands/reads.py:103`, `muc_one_up/cli/commands/reads.py:256`, `muc_one_up/cli/commands/reads.py:482`, `muc_one_up/cli/commands/analyze.py:89`, `muc_one_up/type_defs.py:114`.
   Impact: the refactor improved command layout, but not the contract. `load_config_raw()` intentionally skips schema validation, then command handlers inject simulator-specific keys into nested dicts. That keeps the CLI flexible, but it makes failures later and less precise, weakens type safety, and reintroduces hidden coupling between command code and config internals.

4. **Medium**: CLI help text has drifted from the actual command surface.
   Evidence: `muc_one_up/cli/click_main.py:47`, `muc_one_up/cli/commands/simulate.py:131`.
   Impact: `simulate` advertises a nonexistent `pipeline` command even though the root help lists only `simulate`, `reads`, and `analyze`. This is small compared with the architecture issues, but it is exactly the sort of UX/documentation drift that causes user confusion and stale docs.

5. **Low-Medium**: The command split introduced a new copy-paste hotspot in `reads.py`.
   Evidence: `muc_one_up/cli/commands/reads.py:103`, `muc_one_up/cli/commands/reads.py:137`, `muc_one_up/cli/commands/reads.py:256`, `muc_one_up/cli/commands/reads.py:292`, `muc_one_up/cli/commands/reads.py:482`, `muc_one_up/cli/commands/reads.py:535`.
   Impact: Illumina, ONT, and PacBio commands now repeat nearly identical batch-processing, output-base derivation, warnings, and source-tracker reconstruction. This is not a correctness bug today, but it raises the cost of keeping the three command variants behaviorally aligned.

## Strengths

- The CLI root is dramatically healthier than before. `muc_one_up/cli/click_main.py:1` is now doing command registration rather than business logic.
- The read-simulation refactor created useful seams. `muc_one_up/read_simulator/assembly_context.py:18`, `muc_one_up/read_simulator/output_config.py:10`, and `muc_one_up/read_simulator/utils/common_utils.py:57` are the right kind of abstractions for this project.
- Typed domain objects are a real improvement. `muc_one_up/type_defs.py:18` and `muc_one_up/type_defs.py:67` are materially better than the pre-refactor tuple/string conventions.
- The project still shows strong test intent and release discipline in configuration. `pyproject.toml:31`, `pyproject.toml:120`, and `pyproject.toml:154` describe a serious toolchain even if the current status is not fully green.
- The split command UX is more discoverable than the old monolith. `muc_one_up/cli/commands/reads.py:18` and `muc_one_up/cli/commands/analyze.py:21` are easier to reason about than the pre-refactor arrangement.

## Recommended Next Steps

1. Break the eager import chain. Lazy-load `Bio`, `orfipy_core`, and read-simulation pipelines at the command/function edge rather than package import time.
2. Fix the `RunResult` contract and make `run_command()`/`run_pipeline()` return shapes unambiguous to mypy and callers.
3. Replace `load_config_raw()` plus ad hoc dict mutation in command handlers with small typed config/view objects for `reads` and `analyze`.
4. Extract the repeated batch/output/source-tracker loop from `muc_one_up/cli/commands/reads.py` into a shared helper before the three subcommands drift.
5. Clean the remaining CLI help drift, starting with the nonexistent `pipeline` reference in `simulate`.
6. Re-run the declared quality gate from a fully provisioned environment and require `ruff`, `mypy`, and `pytest` all to pass together before the next release tag.

## Conclusion

The refactor was worthwhile. Relative to the 2026-04-01 review, the codebase is meaningfully better structured and easier to navigate, and the old CLI monolith problem is genuinely fixed.

The current repository still does not meet the bar of a well-maintained Python CLI project because the new boundaries are not fully hardened: import-time dependency loading regressed, mypy is red, and the requested pytest run does not collect in this environment. My overall assessment moves from **5/10** to **6/10**: improved, but not yet cleanly stable at the package, typing, and tooling boundaries.
