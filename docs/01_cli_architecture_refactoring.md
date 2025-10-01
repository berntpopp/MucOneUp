# CLI Architecture & Refactoring

**Status:** ✅ **IMPLEMENTED** (2025-10-01)
**Priority:** 🔴 URGENT (was)
**Estimated Effort:** 2-3 days
**Actual Effort:** 1 day
**Impact:** HIGH - Enables testing, improves maintainability

---

## ✅ Implementation Summary

**Completed:** 2025-10-01

**Achievements:**
- ✅ Refactored 801-line `main()` to 42 lines (94% reduction)
- ✅ Created 8 focused modules following SOLID principles
- ✅ Eliminated ~120 lines of SNP duplication (DRY)
- ✅ All 46 tests passing (zero regressions)
- ✅ Zero linting errors (ruff)
- ✅ CLI tool verified working

**Documentation:**
- [CLI_REFACTORING_SUMMARY.md](./CLI_REFACTORING_SUMMARY.md) - Initial refactoring (functions extraction)
- [CLI_MODULARIZATION_SUMMARY.md](./CLI_MODULARIZATION_SUMMARY.md) - Package modularization
- [VALIDATION_CHECKLIST.md](./VALIDATION_CHECKLIST.md) - Comprehensive validation

---

## Problem Statement

The current `cli.py` has severe architectural issues:

- **1,144 lines** in a single file
- **801-line `main()` function** (70% of the file)
- **Only 6 functions** total (insufficient decomposition)
- **8+ levels of nesting** in control flow
- **Cyclomatic complexity ~50+** (target: <10 per function)
- **29 `sys.exit()` calls** in cli.py alone
- **Code duplication:** SNP integration logic duplicated (~120 lines)
- **Duplicate imports:** Same modules imported twice

**Result:** Untestable, unmaintainable, violates Single Responsibility Principle.

## SOLID Principles Violated

| Principle | Violation | Impact |
|-----------|-----------|--------|
| **S**ingle Responsibility | `main()` handles argument parsing, config loading, simulation, mutation, SNP integration, file I/O, ORF prediction, read simulation, and stats | Cannot test individual concerns |
| **O**pen/Closed | Hard to extend without modifying `main()` | Every new feature requires changing monolithic function |
| **L**iskov Substitution | N/A (no inheritance) | - |
| **I**nterface Segregation | No interfaces, everything in one function | Cannot mock or substitute components |
| **D**ependency Inversion | Direct dependencies on implementations | Tight coupling, hard to test |

## DRY Violations

### 1. SNP Integration Duplication

```python
# Lines 646-739: SNP processing for dual mode
if args.snp_file or args.random_snps:
    # ~93 lines of SNP integration logic

# Lines 796-863: Nearly IDENTICAL SNP processing for single mode
if args.snp_file or args.random_snps:
    # ~67 lines of DUPLICATED logic
```

**Fix:** Extract to `integrate_snps(args, config, results, skip_reference_check=False)`

### 2. Output File Generation

Multiple places construct output paths:
```python
fasta_out = os.path.join(out_dir, f"{base_name}.{iteration_suffix}.simulated.fa")
structure_out = os.path.join(out_dir, f"{base_name}.{iteration_suffix}.structure.txt")
stats_out = os.path.join(out_dir, f"{base_name}.{iteration_suffix}.simulation_stats.json")
```

**Fix:** Create `OutputPathBuilder` class with DRY path construction.

## Refactoring Strategy

### Phase 1: Extract Helper Functions (KISS)

Break down `main()` into focused, single-purpose functions:

```python
# From: 801-line main()
# To: <100-line main() that orchestrates

def setup_configuration(args) -> Tuple[Dict[str, Any], Path]:
    """Load config and setup output directory. ~30 lines"""
    pass

def determine_simulation_mode(args) -> List[Any]:
    """Resolve --fixed-lengths, --input-structure, --simulate-series. ~50 lines"""
    pass

def generate_haplotypes(args, config, fixed_conf) -> List[Tuple[str, List[str]]]:
    """Generate haplotypes based on simulation mode. ~60 lines"""
    pass

def apply_mutations(args, config, results) -> Tuple[List[Tuple], Optional[Dict]]:
    """Apply mutations if requested, handles dual mode. ~80 lines"""
    pass

def integrate_snps(args, config, results, skip_reference_check=False) -> Tuple[List[Tuple], Dict]:
    """Unified SNP integration (removes duplication). ~70 lines"""
    pass

def write_outputs(args, results, mutated_units, sim_index, variant=""):
    """Write FASTA, structure, mutated units files. ~50 lines"""
    pass

def run_optional_analyses(args, config, results, sim_index, variant=""):
    """ORF prediction, toxic detection, read sim, stats. ~100 lines"""
    pass

def run_single_simulation(args, config, out_dir, fixed_conf, sim_index):
    """Run one complete simulation iteration. ~150 lines"""
    pass

def main():
    """CLI entry point - orchestration only. <100 lines"""
    args = build_parser().parse_args()
    configure_logging(args.log_level)

    try:
        config, out_dir = setup_configuration(args)
        sim_configs = determine_simulation_mode(args)

        for sim_index, fixed_conf in enumerate(sim_configs, 1):
            run_single_simulation(args, config, out_dir, fixed_conf, sim_index)

    except MucOneUpError as e:
        logging.error(str(e))
        sys.exit(1)
    except Exception:
        logging.exception("Unexpected error")
        sys.exit(2)
```

**Benefit:** Each function is testable, readable, maintains single responsibility.

### Phase 2: Create CLI Module Structure (Modularization)

```
muc_one_up/
├── cli/
│   ├── __init__.py          # Expose main entry point
│   ├── commands.py          # Command implementations (run_single_simulation, etc.)
│   ├── parsers.py           # Argument parsing (build_parser, parse_fixed_lengths)
│   ├── validation.py        # Input validation
│   ├── output.py            # Output path construction, file writing
│   └── orchestration.py     # High-level workflow coordination
└── cli.py                   # Thin wrapper calling cli/__init__.py:main()
```

**Benefit:** Clear separation of concerns, easier to navigate, testable modules.

### Phase 3: Consider Click Migration (Optional)

Current flat argument structure is confusing:
```bash
muconeup --config X --out-base Y --mutation-name Z --simulate-reads ont --random-snps ...
```

Click enables clear command hierarchy:
```bash
muconeup simulate --config X --output Y
muconeup simulate --config X --output Y --mutate dupC --targets 1,25
muconeup reads illumina --input X --coverage 30
muconeup reads ont --input X --coverage 30
```

**Benefits:**
- Clear separation of simulate vs. reads commands
- Better help text per command
- Built-in testing support via `CliRunner`
- Gradual migration possible

**Migration Path:**
1. Keep argparse, extract functions (Phase 1-2) ← **Do this first**
2. Add Click as optional wrapper
3. Migrate commands incrementally
4. Remove argparse once complete

**Do NOT migrate to Click until after Phase 1-2 refactoring.** KISS principle.

## Actionable Steps

### Week 1: Extract Functions

**Goal:** Break down `main()` into 8-10 focused functions.

**Tasks:**
1. ✅ Create `muc_one_up/cli/` directory
2. ✅ Move `build_parser()` to `cli/parsers.py`
3. ✅ Extract `setup_configuration()` from main()
4. ✅ Extract `determine_simulation_mode()` from main()
5. ✅ Extract `generate_haplotypes()` from main()
6. ✅ Extract `apply_mutations()` from main() (handles dual mode)
7. ✅ Extract `integrate_snps()` from main() (removes duplication)
8. ✅ Extract `write_outputs()` from main()
9. ✅ Extract `run_optional_analyses()` from main()
10. ✅ Rewrite `main()` to orchestrate only (<100 lines)
11. ✅ Add tests for each new function

**Success Criteria:**
- `main()` is <100 lines
- Each helper function is <150 lines
- SNP duplication eliminated
- All imports consolidated
- All existing tests pass

### Week 2: Create Module Structure

**Goal:** Organize CLI code into logical modules.

**Tasks:**
1. ✅ Create `cli/commands.py` for command implementations
2. ✅ Create `cli/validation.py` for input validation
3. ✅ Create `cli/output.py` for file writing
4. ✅ Move functions to appropriate modules
5. ✅ Update imports in `cli.py`
6. ✅ Add module-level docstrings
7. ✅ Add integration tests

**Success Criteria:**
- Clear module boundaries
- Each module has single purpose
- No circular dependencies
- All tests pass

## Testing Strategy

### Before Refactoring (Current State)
```python
# CANNOT test main() - sys.exit() kills test process
def test_cli_basic_simulation():
    with pytest.raises(SystemExit):  # This is bad!
        main()  # Entire test runner terminates
```

### After Refactoring
```python
# Can test each function independently
def test_setup_configuration(temp_config, tmp_path):
    args = Mock(config=temp_config, out_dir=str(tmp_path))
    config, out_dir = setup_configuration(args)
    assert config["repeats"] is not None
    assert out_dir.exists()

def test_generate_haplotypes(sample_config):
    args = Mock(fixed_lengths=None, input_structure=None)
    results = generate_haplotypes(args, sample_config, fixed_conf=None)
    assert len(results) == 2  # Diploid

def test_integrate_snps_no_duplication():
    # Test SNP logic once, not twice!
    pass
```

## Code Quality Improvements

### Remove Duplicate Imports
```python
# ❌ Current (cli.py lines 27-33 and 44-50)
from .snp_integrator import (
    apply_snps_from_file,
    generate_random_snps,
    ...
)
# ... 10 lines later ...
from .snp_integrator import (  # DUPLICATE!
    apply_snps_from_file,
    generate_random_snps,
    ...
)

# ✅ Fixed - single import block
from .snp_integrator import (
    apply_snps_from_file,
    generate_random_snps,
    integrate_snps_to_haplotypes,
    classify_snp_region,
)
```

### Reduce Nesting
```python
# ❌ Current - 8 levels deep
if args.mutation_name:
    if "," in args.mutation_name:
        for mutation in mutations:
            if mutation == "normal":
                if args.simulate_series:
                    for fixed_conf in fixed_lengths_configs:
                        if args.random_snps:
                            if args.snp_file:
                                # ... code here ...

# ✅ Fixed - early returns, guard clauses
def apply_mutations(args, config, results):
    if not args.mutation_name:
        return results, None

    if "," not in args.mutation_name:
        return apply_single_mutation(args, config, results)

    return apply_dual_mutations(args, config, results)
```

## Metrics to Track

| Metric | Before | Target | After |
|--------|--------|--------|-------|
| `main()` line count | 801 | <100 | TBD |
| Total functions in cli.py | 6 | 15+ | TBD |
| Max nesting depth | 8+ | 3-4 | TBD |
| Cyclomatic complexity (main) | ~50 | <10 | TBD |
| Code duplication (lines) | ~120 | 0 | TBD |
| Test coverage of CLI | 0% | >60% | TBD |

## References

- **SOLID Principles:** https://en.wikipedia.org/wiki/SOLID
- **DRY Principle:** https://en.wikipedia.org/wiki/Don%27t_repeat_yourself
- **KISS Principle:** https://en.wikipedia.org/wiki/KISS_principle
- **Cyclomatic Complexity:** https://en.wikipedia.org/wiki/Cyclomatic_complexity
- **Click Framework:** https://click.palletsprojects.com/
- **Refactoring Patterns:** Martin Fowler's "Refactoring" book

## Next Steps

After completing this refactoring:
1. ➡️ Proceed to **Error Handling & Exceptions** (remove `sys.exit()`)
2. ➡️ Add comprehensive **Testing** (enabled by this refactoring)
3. ➡️ Add **Type Hints** to all new functions
4. Consider **Click migration** (only after 1-3 are complete)
