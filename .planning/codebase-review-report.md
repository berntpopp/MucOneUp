# Codebase Review Report

Date: 2026-04-01
Repository: `MucOneUp`

## Scope

This review focused on:

- CLI architecture and command implementation
- Core simulation and mutation domain logic
- Read-simulation pipelines and external-tool wrapper boundaries
- Cross-cutting concerns such as configuration, typing, testability, and robustness

The review criteria were based on practical code-review heuristics around DRY, KISS, SOLID, modularization, correctness risk, and maintainability, informed by:

- Refactoring Guru code smell guidance
- Click command/group structure guidance
- Python Packaging guidance for CLI entry points

## Overall Rating

- Architecture: `4.5/10`
- Modularity: `5/10`
- DRY: `4/10`
- KISS: `5/10`
- SOLID: `4/10`
- Type safety: `4/10`
- Testability: `3.5/10`
- Operational robustness: `4/10`
- CLI UX: `6/10`
- Tooling/release hygiene: `6/10`
- Overall: `5/10`

## Executive Summary

This is a serious codebase with substantial functionality, meaningful test coverage, and clear domain intent. It is not a rewrite candidate.

The main issue is architectural drift. The CLI layer has become a god module, domain logic relies on weakly typed shared dictionaries and string conventions, and the read-simulation stack has inconsistent configuration and tool-execution boundaries. The system still looks recoverable, but further feature work without refactoring will make it steadily harder to test and reason about.

## Key Findings

### 1. Import-time CLI coupling breaks basic operability

Severity: Critical

Importing the CLI package eagerly imports the full Click application and transitive dependencies. That currently hard-fails if `rfc8785` is unavailable.

- `muc_one_up/cli/__init__.py:10`
- `muc_one_up/cli/click_main.py:30`
- `muc_one_up/config_fingerprint.py:37`

Impact:

- CLI import can fail before command execution
- test collection fails in minimal environments
- dependency boundaries are too broad for a package entrypoint

### 2. `click_main.py` is a god module

Severity: High

The CLI entrypoint owns too many responsibilities: Click wiring, JSON loading, backend preparation, subprocess execution, FASTA parsing, output naming, and analysis flow.

- `muc_one_up/cli/click_main.py:472`
- `muc_one_up/cli/click_main.py:625`
- `muc_one_up/cli/click_main.py:807`
- `muc_one_up/cli/click_main.py:1019`
- `muc_one_up/cli/click_main.py:1173`
- `muc_one_up/cli/click_main.py:1299`
- `muc_one_up/cli/click_main.py:1455`

Additional signal:

- `muc_one_up/cli/click_main.py` is 1586 lines

Impact:

- poor SRP
- high cost of change
- expensive command-level testing
- growing duplication between commands

### 3. CLI contracts do not consistently match backend behavior

Severity: High

The read commands expose `--out-dir` and `--out-base`, but those values are not reliably propagated into the underlying simulator outputs. They affect logs and some side files more than the actual simulated read outputs.

- `muc_one_up/cli/click_main.py:434`
- `muc_one_up/cli/click_main.py:539`
- `muc_one_up/cli/click_main.py:587`
- `muc_one_up/cli/click_main.py:691`
- `muc_one_up/cli/click_main.py:739`
- `muc_one_up/cli/click_main.py:931`
- `muc_one_up/read_simulation.py:95`

Impact:

- misleading CLI UX
- hidden surprises for users
- fragile maintenance because the interface promises more control than the backend honors

### 4. ORF/toxic-protein workflows are duplicated and drifting

Severity: High

ORF analysis behavior exists in more than one place and has already diverged. Parts of the code read constants as if they were flat keys, while other parts use assembly-keyed constants.

- `muc_one_up/cli/analysis.py:22`
- `muc_one_up/cli/analysis.py:60`
- `muc_one_up/cli/analysis.py:125`
- `muc_one_up/cli/click_main.py:1045`
- `muc_one_up/cli/click_main.py:1138`
- `muc_one_up/cli/orchestration.py:141`

Impact:

- DRY violation
- higher correctness risk
- future changes likely to diverge further

### 5. CLI orchestration is carrying domain logic

Severity: Medium-High

`run_single_simulation_iteration()` is nominally an orchestrator, but it also reshapes data, derives mutation metadata, constructs source trackers, and makes output decisions.

- `muc_one_up/cli/orchestration.py:24`
- `muc_one_up/cli/orchestration.py:125`

Impact:

- weak cohesion
- hard-to-mock transaction script
- reduced modularity and testability

### 6. The core domain is coupled to raw config dictionaries

Severity: High

Simulation and mutation code depend directly on wide `dict[str, Any]` structures rather than focused typed domain objects.

- `muc_one_up/type_defs.py:21`
- `muc_one_up/config.py:118`
- `muc_one_up/config.py:476`
- `muc_one_up/simulate.py:102`
- `muc_one_up/simulate.py:253`
- `muc_one_up/mutate.py:224`

Impact:

- violates DIP in practice
- spreads infrastructure-shaped knowledge through core logic
- makes unit tests broader than necessary

### 7. Domain state is encoded implicitly in strings and tuples

Severity: High

Mutation state is represented with a `"m"` string suffix, mutation targets are 1-based positional tuples, and haplotypes are loose tuples/lists.

- `muc_one_up/type_defs.py:15`
- `muc_one_up/type_defs.py:29`
- `muc_one_up/simulate.py:113`
- `muc_one_up/simulate.py:161`
- `muc_one_up/mutate.py:157`
- `muc_one_up/mutate.py:196`
- `muc_one_up/validation.py:147`

Impact:

- primitive obsession
- hidden invariants
- higher regression risk when behavior evolves

### 8. Mutation handling has weak source-of-truth boundaries

Severity: High

Mutation logic maintains both `seq` and `chain` representations and updates them through mixed strategies. Assembly rules are not centralized.

- `muc_one_up/mutate.py:149`
- `muc_one_up/mutate.py:185`
- `muc_one_up/mutate.py:211`
- `muc_one_up/mutate.py:233`
- `muc_one_up/mutate.py:344`
- `muc_one_up/simulate.py:85`
- `muc_one_up/simulate.py:289`

Impact:

- duplicate representations can drift
- correctness depends on conventions
- difficult to refactor safely

### 9. Read simulation resolves assembly inconsistently

Severity: High

The Illumina pipeline reads the active assembly from more than one config location during a single run.

- `muc_one_up/read_simulator/pipeline.py:189`
- `muc_one_up/read_simulator/pipeline.py:291`

Impact:

- potential mismatch between sample BAM assembly and alignment reference
- correctness issue, not just style

### 10. External tool boundaries are inconsistent and duplicated

Severity: Medium-High

Some wrappers use shared helpers, others bypass them with ad hoc subprocess handling. Long-read alignment also appears duplicated across wrappers.

Representative refs:

- `muc_one_up/read_simulator/utils/common_utils.py:20`
- `muc_one_up/read_simulator/wrappers/bwa_wrapper.py:55`
- `muc_one_up/read_simulator/wrappers/samtools_wrapper.py:61`
- `muc_one_up/read_simulator/wrappers/nanosim_wrapper.py:144`
- `muc_one_up/read_simulator/wrappers/minimap2_wrapper.py:59`
- `muc_one_up/read_simulator/ont_pipeline.py:297`

Impact:

- inconsistent timeout and cleanup behavior
- weaker observability
- duplicated execution semantics

### 11. Temporary-file lifetime management is fragile

Severity: High

The ONT split-simulation path returns FASTQ paths from inside a temporary directory and later code consumes them after that directory scope may already have ended.

- `muc_one_up/read_simulator/utils/diploid_handler.py:221`
- `muc_one_up/read_simulator/utils/diploid_handler.py:323`
- `muc_one_up/read_simulator/ont_pipeline.py:329`

Impact:

- fragile behavior by construction
- intermittent failures likely

### 12. Randomness is global and not injectable

Severity: Medium

The code uses module-level `random.*` and seeds the global RNG.

- `muc_one_up/simulate.py:80`
- `muc_one_up/simulate.py:211`
- `muc_one_up/simulate.py:271`
- `muc_one_up/distribution.py:41`
- `muc_one_up/probabilities.py:44`
- `muc_one_up/mutate.py:175`

Impact:

- hidden coupling between tests and workflows
- weaker reproducibility guarantees
- harder composition of simulation steps

### 13. Error handling is inconsistent and too broad

Severity: Medium

Commands mix continue-on-error and fail-fast behavior, and many paths collapse all failures into generic exit codes after broad exception handling.

- `muc_one_up/cli/click_main.py:578`
- `muc_one_up/cli/click_main.py:730`
- `muc_one_up/cli/click_main.py:970`
- `muc_one_up/cli/click_main.py:1104`
- `muc_one_up/cli/click_main.py:1151`
- `muc_one_up/cli/analysis.py:105`
- `muc_one_up/cli/analysis.py:152`

Impact:

- uneven CLI behavior
- imprecise failure semantics
- weaker tests for error paths

## Strengths

- The repository has clear domain decomposition: CLI, simulation, analysis, read simulation, wrappers, tests.
- `pyproject.toml` is generally well structured and the console entry point is correctly declared.
- Linting discipline is present and `ruff` passes.
- There is substantial documentation and test coverage, which lowers refactoring risk.
- The code documents biological rules and workflow intent clearly.
- Several reusable seams already exist and can be consolidated rather than invented from scratch.

## Verification Notes

Commands run during review:

- `ruff check muc_one_up tests` -> passed
- `mypy muc_one_up` -> failed due missing import/stub handling around `rfc8785`
- `pytest -q tests/cli/test_help_version_access.py tests/test_cli.py tests/test_simulate.py` -> failed during collection because CLI import transitively requires `rfc8785`

## Recommended Change Order

### Phase 1: Stabilize basic operability

1. Decouple CLI imports from optional-heavy transitive modules.
2. Fix `rfc8785` dependency handling so import, mypy, and basic CLI tests work in minimal environments.
3. Add smoke tests for `muc_one_up.cli` import and `muconeup --help`.

### Phase 2: Fix the CLI boundary

1. Split `muc_one_up/cli/click_main.py` into a thin root app plus per-command modules.
2. Replace the fake argparse namespace with typed command option objects.
3. Make CLI options map directly and consistently to backend behavior.
4. Centralize error mapping at the CLI boundary.

### Phase 3: Repair domain modeling

1. Introduce typed domain models for simulation config, haplotypes, and mutation targets.
2. Replace mutation string suffixes and positional tuples with explicit data structures.
3. Centralize sequence/chain assembly so there is one source of truth.
4. Inject RNG objects rather than using global `random`.

### Phase 4: Consolidate read-simulation infrastructure

1. Unify assembly resolution.
2. Centralize subprocess/tool execution behind one runner.
3. Remove duplicated minimap2 and wrapper execution paths.
4. Fix temporary-directory lifetime issues in ONT handling.

## Recommended Refactoring Targets

- First target: `muc_one_up/cli/click_main.py`
- Second target: `muc_one_up/config_fingerprint.py` import boundary
- Third target: `muc_one_up/cli/orchestration.py`
- Fourth target: `muc_one_up/simulate.py` and `muc_one_up/mutate.py`
- Fifth target: `muc_one_up/read_simulator/pipeline.py`

## Conclusion

The repository has strong domain value and enough structure to improve incrementally. The main technical debt is not lack of functionality, but weak boundaries: the CLI does too much, the domain model is too implicit, and the read-simulation stack has inconsistent contracts.

The highest-value next move is to refactor the CLI boundary first. That will reduce coupling, improve testability, and make subsequent domain and pipeline cleanup much easier to perform safely.
