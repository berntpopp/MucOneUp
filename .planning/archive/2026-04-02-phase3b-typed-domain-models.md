# Phase 3B: Typed Domain Models — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace implicit string conventions (`"m"` suffix) and raw tuples (`tuple[str, list[str]]`) with explicit typed dataclasses (`RepeatUnit`, `HaplotypeResult`, `MutationTarget`), eliminating fragile string manipulation throughout the domain layer.

**Architecture:** Introduce three dataclasses in `muc_one_up/type_defs.py`. `RepeatUnit(symbol, mutated)` replaces the `str + "m"` convention. `HaplotypeResult(sequence, chain)` replaces `tuple[str, list[str]]`. `MutationTarget(haplotype_index, repeat_index)` replaces `tuple[int, int]`. Each type provides conversion helpers for backward-compatible serialization (structure files still use `"1-2-Xm-9"` format). Migration proceeds bottom-up: type_defs → domain modules → CLI layer → peripheral modules.

**Tech Stack:** Python 3.10+, dataclasses, pytest

**Spec:** `.planning/specs/2026-04-01-codebase-refactoring-design.md` (sections 3.1, 3.2, 3.5)

**Deferred to Phase 3C:** `SimulationConfig` dataclass (typed wrapper for `ConfigDict`) — too large a change to combine with the type migration; every domain function signature would change simultaneously.

---

## File Structure

### Modified Files
| File | Changes |
|------|---------|
| `muc_one_up/type_defs.py` | Add `RepeatUnit`, `HaplotypeResult`, `MutationTarget` dataclasses; update type aliases |
| `muc_one_up/assembly.py` | Accept `list[RepeatUnit]` chains; remove `.rstrip("m")` |
| `muc_one_up/simulate.py` | Build chains as `list[RepeatUnit]`; return `HaplotypeResult` |
| `muc_one_up/mutate.py` | Use `RepeatUnit.mutated` flag; accept `list[MutationTarget]`; return `HaplotypeResult` |
| `muc_one_up/distribution.py` | No changes (no chain/tuple usage) |
| `muc_one_up/probabilities.py` | No changes (no chain/tuple usage) |
| `muc_one_up/io.py` | Parse structure files into `list[RepeatUnit]` |
| `muc_one_up/simulation_statistics.py` | Use `RepeatUnit` attributes instead of `.rstrip("m")` / `.endswith("m")` |
| `muc_one_up/cli/haplotypes.py` | Update return type to `list[HaplotypeResult]` |
| `muc_one_up/cli/mutations.py` | Use `MutationTarget`; update tuple unpacking |
| `muc_one_up/cli/outputs.py` | Use `HaplotypeResult` attributes; serialize chains via `RepeatUnit.__str__` |
| `muc_one_up/cli/orchestration.py` | Update tuple unpacking to use `HaplotypeResult` |
| `muc_one_up/cli/snps.py` | Update type signatures |
| `muc_one_up/translate.py` | Update type signatures |
| `muc_one_up/snp_integrator.py` | Update type signatures |
| `muc_one_up/validation.py` | Use `MutationTarget` type |
| `muc_one_up/read_simulator/source_tracking.py` | Use `RepeatUnit` attributes |

### New Test Files
| File | Responsibility |
|------|---------------|
| `tests/test_domain_models.py` | Tests for `RepeatUnit`, `HaplotypeResult`, `MutationTarget` |

---

## Task 1: Create RepeatUnit Dataclass

**Files:**
- Modify: `muc_one_up/type_defs.py`
- Create: `tests/test_domain_models.py`

- [ ] **Step 1: Write failing tests for RepeatUnit**

```python
# tests/test_domain_models.py
"""Tests for typed domain models."""

from __future__ import annotations

import pytest

from muc_one_up.type_defs import RepeatUnit


class TestRepeatUnit:
    """Tests for RepeatUnit dataclass."""

    def test_create_normal_repeat(self):
        ru = RepeatUnit(symbol="1", mutated=False)
        assert ru.symbol == "1"
        assert ru.mutated is False

    def test_create_mutated_repeat(self):
        ru = RepeatUnit(symbol="X", mutated=True)
        assert ru.symbol == "X"
        assert ru.mutated is True

    def test_str_normal(self):
        """Normal repeat serializes as just the symbol."""
        ru = RepeatUnit(symbol="1", mutated=False)
        assert str(ru) == "1"

    def test_str_mutated(self):
        """Mutated repeat serializes with 'm' suffix."""
        ru = RepeatUnit(symbol="X", mutated=True)
        assert str(ru) == "Xm"

    def test_from_str_normal(self):
        ru = RepeatUnit.from_str("1")
        assert ru.symbol == "1"
        assert ru.mutated is False

    def test_from_str_mutated(self):
        ru = RepeatUnit.from_str("Xm")
        assert ru.symbol == "X"
        assert ru.mutated is True

    def test_from_str_multi_char_symbol(self):
        ru = RepeatUnit.from_str("6p")
        assert ru.symbol == "6p"
        assert ru.mutated is False

    def test_from_str_multi_char_mutated(self):
        ru = RepeatUnit.from_str("6pm")
        assert ru.symbol == "6p"
        assert ru.mutated is True

    def test_equality(self):
        assert RepeatUnit("1", False) == RepeatUnit("1", False)
        assert RepeatUnit("1", True) != RepeatUnit("1", False)

    def test_chain_to_str_list(self):
        """Helper to serialize a chain for structure file output."""
        chain = [RepeatUnit("1", False), RepeatUnit("X", True), RepeatUnit("9", False)]
        assert [str(ru) for ru in chain] == ["1", "Xm", "9"]

    def test_chain_from_str_list(self):
        """Helper to parse a chain from structure file input."""
        raw = ["1", "Xm", "9"]
        chain = [RepeatUnit.from_str(s) for s in raw]
        assert chain[0] == RepeatUnit("1", False)
        assert chain[1] == RepeatUnit("X", True)
        assert chain[2] == RepeatUnit("9", False)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_domain_models.py -v`
Expected: FAIL with ImportError

- [ ] **Step 3: Implement RepeatUnit**

Add to `muc_one_up/type_defs.py` (after the imports, before existing aliases):

```python
from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class RepeatUnit:
    """A single repeat unit in a VNTR chain.

    Replaces the convention of using plain strings with 'm' suffix
    for mutated repeats. Provides typed access to the symbol and
    mutation status.

    Attributes:
        symbol: The repeat type identifier (e.g., "1", "6p", "X").
        mutated: Whether this repeat has been mutated.
    """

    symbol: str
    mutated: bool = False

    def __str__(self) -> str:
        """Serialize to legacy string format (e.g., 'Xm' if mutated)."""
        return f"{self.symbol}m" if self.mutated else self.symbol

    @classmethod
    def from_str(cls, s: str) -> RepeatUnit:
        """Parse from legacy string format (e.g., 'Xm' -> RepeatUnit('X', True))."""
        if s.endswith("m"):
            return cls(symbol=s[:-1], mutated=True)
        return cls(symbol=s, mutated=False)
```

Update the `RepeatChain` alias:

```python
# Type aliases for haplotype representation — legacy (will migrate consumers)
RepeatChain = list[str]
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_domain_models.py -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/type_defs.py tests/test_domain_models.py
git commit -m "feat: add RepeatUnit dataclass replacing 'm' suffix convention"
```

---

## Task 2: Create HaplotypeResult and MutationTarget Dataclasses

**Files:**
- Modify: `muc_one_up/type_defs.py`
- Modify: `tests/test_domain_models.py`

- [ ] **Step 1: Write failing tests**

Append to `tests/test_domain_models.py`:

```python
from muc_one_up.type_defs import HaplotypeResult, MutationTarget


class TestHaplotypeResult:
    """Tests for HaplotypeResult dataclass."""

    def test_create(self):
        chain = [RepeatUnit("1", False), RepeatUnit("9", False)]
        hr = HaplotypeResult(sequence="ACGT", chain=chain)
        assert hr.sequence == "ACGT"
        assert len(hr.chain) == 2

    def test_chain_symbols(self):
        """chain_strs() returns legacy string list for backward compat."""
        chain = [RepeatUnit("1", False), RepeatUnit("X", True), RepeatUnit("9", False)]
        hr = HaplotypeResult(sequence="ACGT", chain=chain)
        assert hr.chain_strs() == ["1", "Xm", "9"]

    def test_as_tuple(self):
        """as_tuple() returns legacy (sequence, chain_strs) format."""
        chain = [RepeatUnit("1", False), RepeatUnit("9", False)]
        hr = HaplotypeResult(sequence="ACGT", chain=chain)
        seq, chain_list = hr.as_tuple()
        assert seq == "ACGT"
        assert chain_list == ["1", "9"]

    def test_from_tuple(self):
        """from_tuple() parses legacy (sequence, chain_strs) format."""
        hr = HaplotypeResult.from_tuple(("ACGT", ["1", "Xm", "9"]))
        assert hr.sequence == "ACGT"
        assert hr.chain[1].symbol == "X"
        assert hr.chain[1].mutated is True


class TestMutationTarget:
    """Tests for MutationTarget dataclass."""

    def test_create(self):
        mt = MutationTarget(haplotype_index=1, repeat_index=25)
        assert mt.haplotype_index == 1
        assert mt.repeat_index == 25

    def test_zero_haplotype_raises(self):
        with pytest.raises(ValueError, match="1-based"):
            MutationTarget(haplotype_index=0, repeat_index=1)

    def test_zero_repeat_raises(self):
        with pytest.raises(ValueError, match="1-based"):
            MutationTarget(haplotype_index=1, repeat_index=0)

    def test_negative_raises(self):
        with pytest.raises(ValueError, match="1-based"):
            MutationTarget(haplotype_index=-1, repeat_index=1)

    def test_as_tuple(self):
        mt = MutationTarget(haplotype_index=1, repeat_index=25)
        assert mt.as_tuple() == (1, 25)

    def test_from_tuple(self):
        mt = MutationTarget.from_tuple((2, 30))
        assert mt.haplotype_index == 2
        assert mt.repeat_index == 30
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_domain_models.py -v`
Expected: FAIL with ImportError

- [ ] **Step 3: Implement HaplotypeResult and MutationTarget**

Add to `muc_one_up/type_defs.py` after `RepeatUnit`:

```python
@dataclass(frozen=True, slots=True)
class HaplotypeResult:
    """Result of simulating a single haplotype.

    Replaces the raw tuple[str, list[str]] convention with named fields.

    Attributes:
        sequence: The assembled DNA sequence.
        chain: The repeat chain as typed RepeatUnit objects.
    """

    sequence: str
    chain: list[RepeatUnit]

    def chain_strs(self) -> list[str]:
        """Return chain as legacy string list (e.g., ['1', 'Xm', '9'])."""
        return [str(ru) for ru in self.chain]

    def as_tuple(self) -> tuple[str, list[str]]:
        """Convert to legacy (sequence, chain_strs) tuple."""
        return (self.sequence, self.chain_strs())

    @classmethod
    def from_tuple(cls, t: tuple[str, list[str]]) -> HaplotypeResult:
        """Parse from legacy (sequence, chain_strs) tuple."""
        return cls(
            sequence=t[0],
            chain=[RepeatUnit.from_str(s) for s in t[1]],
        )


@dataclass(frozen=True, slots=True)
class MutationTarget:
    """A target position for mutation application.

    Replaces raw tuple[int, int] with named, validated fields.
    Uses 1-based indexing matching biological conventions.

    Attributes:
        haplotype_index: 1-based haplotype number.
        repeat_index: 1-based repeat position within the chain.
    """

    haplotype_index: int
    repeat_index: int

    def __post_init__(self) -> None:
        if self.haplotype_index < 1 or self.repeat_index < 1:
            raise ValueError(
                f"MutationTarget uses 1-based indexing, got "
                f"haplotype_index={self.haplotype_index}, "
                f"repeat_index={self.repeat_index}"
            )

    def as_tuple(self) -> tuple[int, int]:
        """Convert to legacy (haplotype_index, repeat_index) tuple."""
        return (self.haplotype_index, self.repeat_index)

    @classmethod
    def from_tuple(cls, t: tuple[int, int]) -> MutationTarget:
        """Parse from legacy (haplotype_index, repeat_index) tuple."""
        return cls(haplotype_index=t[0], repeat_index=t[1])
```

Update existing aliases — keep old aliases as deprecated references:

```python
# Legacy type aliases (kept for gradual migration)
RepeatChain = list[str]
Haplotype = tuple[str, RepeatChain]
HaplotypeList = list[Haplotype]
MutationTargets = list[tuple[int, int]]
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_domain_models.py -v`
Expected: All PASS

- [ ] **Step 5: Run full test suite to verify no regressions**

Run: `pytest --tb=short -q`
Expected: All PASS (existing code still uses legacy types)

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/type_defs.py tests/test_domain_models.py
git commit -m "feat: add HaplotypeResult and MutationTarget dataclasses"
```

---

## Task 3: Migrate assembly.py to Use RepeatUnit

**Files:**
- Modify: `muc_one_up/assembly.py`
- Modify: `tests/test_assembly.py`

- [ ] **Step 1: Update assemble_sequence to accept list[RepeatUnit]**

In `muc_one_up/assembly.py`:

1. Update imports:
```python
from .type_defs import ConfigDict, DNASequence, RepeatUnit
```

2. Change function signature and body — accept `list[RepeatUnit]` instead of `list[str]`:

```python
def assemble_sequence(chain: list[RepeatUnit], config: ConfigDict) -> DNASequence:
    """Assemble a DNA sequence from a repeat chain and flanking constants.

    Concatenates: left_constant + repeat_units + right_constant.
    The right constant is only appended if the chain ends with repeat "9"
    (the canonical terminal repeat in MUC1 VNTR).

    Args:
        chain: List of RepeatUnit objects.
        config: Configuration dict with 'repeats', 'constants', and
                optionally 'reference_assembly' keys.

    Returns:
        Assembled DNA sequence string.

    Raises:
        KeyError: If a repeat symbol is not in config,
            or if required constants for the reference assembly are missing.
    """
    repeats_dict = config["repeats"]
    ref_assembly = config.get("reference_assembly", "hg38")
    constants = config["constants"][ref_assembly]
    left_const = constants["left"]
    right_const = constants["right"]

    # Build repeat region
    parts: list[str] = [left_const]
    for unit in chain:
        if unit.symbol not in repeats_dict:
            raise KeyError(
                f"Repeat symbol '{unit.symbol}' not found in config repeats"
            )
        parts.append(repeats_dict[unit.symbol])

    # Right constant only if chain ends with canonical terminal repeat "9"
    if chain and chain[-1].symbol == "9":
        parts.append(right_const)
    elif chain:
        logger.debug("Chain does not end with '9'; right constant omitted.")

    return "".join(parts)
```

- [ ] **Step 2: Update tests/test_assembly.py to use RepeatUnit**

Replace string chains with `RepeatUnit` objects throughout. Example for `test_simple_chain`:

```python
from muc_one_up.type_defs import RepeatUnit

# In each test, replace:
#   chain = ["1", "2", "9"]
# With:
#   chain = [RepeatUnit("1"), RepeatUnit("2"), RepeatUnit("9")]
#
# And for mutated:
#   chain = ["1", "Xm", "9"]
# With:
#   chain = [RepeatUnit("1"), RepeatUnit("X", mutated=True), RepeatUnit("9")]
```

Update all test methods in `TestAssembleSequence`:
- `test_simple_chain`: chain = `[RepeatUnit("1"), RepeatUnit("2"), RepeatUnit("9")]`
- `test_mutation_marker_stripped`: chain = `[RepeatUnit("1"), RepeatUnit("X", True), RepeatUnit("9")]`
- `test_empty_chain`: chain = `[]` (unchanged)
- `test_chain_not_ending_with_9_omits_right_constant`: chain = `[RepeatUnit("1"), RepeatUnit("2")]`
- `test_unknown_symbol_raises`: chain = `[RepeatUnit("UNKNOWN")]`
- `test_hg19_assembly`: chain = `[RepeatUnit("1")]`
- `test_single_repeat_ending_with_9`: chain = `[RepeatUnit("9")]`

- [ ] **Step 3: Update simulate.py and mutate.py thin wrappers**

In `muc_one_up/simulate.py`, update `assemble_haplotype_from_chain` to convert legacy `list[str]` to `list[RepeatUnit]` before delegating:

```python
def assemble_haplotype_from_chain(chain: RepeatChain, config: ConfigDict) -> DNASequence:
    """Assemble complete haplotype sequence from repeat chain.

    Delegates to assembly.assemble_sequence().
    Converts legacy string chain to RepeatUnit list.
    """
    from .type_defs import RepeatUnit

    typed_chain = [RepeatUnit.from_str(s) for s in chain]
    return assemble_sequence(typed_chain, config)
```

In `muc_one_up/mutate.py`, update `rebuild_haplotype_sequence` similarly:

```python
def rebuild_haplotype_sequence(chain: RepeatChain, config: ConfigDict) -> DNASequence:
    """Rebuild haplotype sequence after mutation.

    Delegates to assembly.assemble_sequence().
    Converts legacy string chain to RepeatUnit list.
    """
    from .type_defs import RepeatUnit

    typed_chain = [RepeatUnit.from_str(s) for s in chain]
    return assemble_sequence(typed_chain, config)
```

- [ ] **Step 4: Run tests**

Run: `pytest tests/test_assembly.py tests/test_simulate.py tests/test_mutate.py --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/assembly.py muc_one_up/simulate.py muc_one_up/mutate.py tests/test_assembly.py
git commit -m "refactor: migrate assembly.py to use RepeatUnit dataclass"
```

---

## Task 4: Migrate simulate.py Chain Building to RepeatUnit

**Files:**
- Modify: `muc_one_up/simulate.py`
- Modify: `tests/test_simulate.py`

- [ ] **Step 1: Update simulate_single_haplotype to build RepeatUnit chains**

In `muc_one_up/simulate.py`:

1. Add import: `from .type_defs import RepeatUnit`

2. In `simulate_single_haplotype()`, change chain building from `list[str]` to `list[RepeatUnit]`. Replace:
   - `repeat_chain = []` stays
   - `repeat_chain.append(current_symbol)` → `repeat_chain.append(RepeatUnit(current_symbol))`
   - `repeat_chain.append(forced_6)` → `repeat_chain.append(RepeatUnit(forced_6))`
   - Same for "7", "8", "9" appends
   - Use `assemble_sequence(repeat_chain, config)` directly instead of inline assembly
   - Return `HaplotypeResult(sequence=assembled_seq, chain=repeat_chain)`

3. Update `simulate_diploid()` to return `list[HaplotypeResult]`:
   - Change return type annotation
   - Replace `haplotypes.append((seq, chain))` with `haplotypes.append(result)`

4. Update `simulate_from_chains()`:
   - Convert incoming `list[str]` chains to `list[RepeatUnit]`
   - Mark mutations with `RepeatUnit(symbol, mutated=True)` instead of `+ "m"`
   - Return `list[HaplotypeResult]`

- [ ] **Step 2: Update tests/test_simulate.py**

Update tests that unpack `(seq, chain)` tuples to use `HaplotypeResult` attributes:
- Replace `for seq, chain in results:` with `for hr in results:` then use `hr.sequence`, `hr.chain`
- Replace `results[0][1]` (chain access) with `results[0].chain`
- Update chain assertions to compare `RepeatUnit` objects

- [ ] **Step 3: Run tests**

Run: `pytest tests/test_simulate.py tests/test_assembly.py --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 4: Run seed reproducibility test**

Run: `pytest tests/test_simulate.py::TestSimulateDiploidAdvanced::test_seed_ensures_reproducibility -v --no-cov`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/simulate.py tests/test_simulate.py
git commit -m "refactor: migrate simulate.py to build RepeatUnit chains and return HaplotypeResult"
```

---

## Task 5: Migrate mutate.py to RepeatUnit and MutationTarget

**Files:**
- Modify: `muc_one_up/mutate.py`
- Modify: `tests/test_mutate.py`

- [ ] **Step 1: Update apply_mutations to use typed models**

In `muc_one_up/mutate.py`:

1. Update imports: add `RepeatUnit, HaplotypeResult, MutationTarget`

2. Change `apply_mutations()` signature:
   - `results: list[HaplotypeResult]` instead of `HaplotypeList`
   - `targets: list[MutationTarget]` instead of `MutationTargets`
   - Return `tuple[list[HaplotypeResult], dict[int, list[tuple[int, str]]]]`

3. Replace "m" suffix logic:
   - `current_symbol.replace("m", "")` → `chain[repeat_index].symbol`
   - `chain[repeat_index] = chain[repeat_index] + "m"` → create new `RepeatUnit(symbol, mutated=True)` and replace in chain (note: chains need to be `list[RepeatUnit]`, mutable)
   - `current_symbol_clean not in allowed_repeats` → `chain[repeat_index].symbol not in allowed_repeats`

4. Update `apply_changes_to_repeat()`:
   - `chain[i].replace("m", "")` → `chain[i].symbol`
   - `chain[repeat_index].replace("m", "")` → `chain[repeat_index].symbol`

5. Iterate `targets` as `MutationTarget` objects:
   - `for hap_i, rep_i in targets:` → `for target in targets:` then use `target.haplotype_index`, `target.repeat_index`

- [ ] **Step 2: Update tests/test_mutate.py**

Update test fixtures and assertions:
- Construct `HaplotypeResult` objects instead of `(seq, chain)` tuples
- Construct `MutationTarget` objects instead of `(int, int)` tuples
- Assert on `RepeatUnit.mutated` instead of `str.endswith("m")`

- [ ] **Step 3: Run tests**

Run: `pytest tests/test_mutate.py tests/test_assembly.py --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 4: Commit**

```bash
git add muc_one_up/mutate.py tests/test_mutate.py
git commit -m "refactor: migrate mutate.py to RepeatUnit, HaplotypeResult, and MutationTarget"
```

---

## Task 6: Migrate CLI Layer to Typed Models

**Files:**
- Modify: `muc_one_up/cli/haplotypes.py`
- Modify: `muc_one_up/cli/mutations.py`
- Modify: `muc_one_up/cli/outputs.py`
- Modify: `muc_one_up/cli/orchestration.py`
- Modify: `muc_one_up/cli/snps.py`

- [ ] **Step 1: Update cli/haplotypes.py**

Change return type from `list[tuple[str, list[str]]]` to `list[HaplotypeResult]`:

```python
from ..type_defs import HaplotypeResult

def generate_haplotypes(args, config, fixed_conf, predefined_chains) -> list[HaplotypeResult]:
```

No other changes needed — `simulate_diploid` and `simulate_from_chains` already return the new types after Task 4.

- [ ] **Step 2: Update cli/mutations.py**

1. `parse_mutation_targets()`: return `list[MutationTarget]`:
```python
from ..type_defs import MutationTarget

def parse_mutation_targets(mutation_targets_list: list[str | tuple[int, int]]) -> list[MutationTarget]:
    # ... parse as before but return MutationTarget objects
    mutation_positions.append(MutationTarget(hap_i, rep_i))
```

2. `find_random_mutation_target()`: accept `list[HaplotypeResult]`, return `list[MutationTarget]`:
```python
def find_random_mutation_target(
    results: list[HaplotypeResult], config: dict, mutation_name: str
) -> list[MutationTarget]:
```
Replace `for hap_idx, (_seq, chain) in enumerate(results, start=1):` with `for hap_idx, hr in enumerate(results, start=1):` and use `hr.chain`.
Replace `pure_sym = sym.replace("m", "")` with `unit.symbol`.

3. `apply_mutation_pipeline()`: update tuple unpacking to use `HaplotypeResult`:
Replace `[(seq, chain.copy()) for seq, chain in results]` with `[HaplotypeResult(hr.sequence, list(hr.chain)) for hr in results]`.

- [ ] **Step 3: Update cli/outputs.py**

Replace all `(seq, chain)` unpacking with `HaplotypeResult` attributes:
- `[seq for seq, chain in results]` → `[hr.sequence for hr in results]`
- `for i, (_sequence, chain) in enumerate(results, start=1):` → `for i, hr in enumerate(results, start=1):`
- `chain_str = "-".join(chain)` → `chain_str = "-".join(str(ru) for ru in hr.chain)`

- [ ] **Step 4: Update cli/orchestration.py**

Replace tuple unpacking:
- `_sequence, chain = result` → use `result.chain`
- Update type comments

- [ ] **Step 5: Update cli/snps.py**

Update type signatures from `list[tuple[str, list[str]]]` to `list[HaplotypeResult]`.
Update unpacking patterns.

- [ ] **Step 6: Run tests**

Run: `pytest tests/ --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 7: Commit**

```bash
git add muc_one_up/cli/haplotypes.py muc_one_up/cli/mutations.py muc_one_up/cli/outputs.py muc_one_up/cli/orchestration.py muc_one_up/cli/snps.py
git commit -m "refactor: migrate CLI layer to HaplotypeResult and MutationTarget"
```

---

## Task 7: Migrate Peripheral Modules

**Files:**
- Modify: `muc_one_up/io.py`
- Modify: `muc_one_up/simulation_statistics.py`
- Modify: `muc_one_up/translate.py`
- Modify: `muc_one_up/snp_integrator.py`
- Modify: `muc_one_up/validation.py`
- Modify: `muc_one_up/read_simulator/source_tracking.py`

- [ ] **Step 1: Update io.py**

In `parse_structure_file()`, convert parsed chains to `list[RepeatUnit]`:
- Replace `chain = chain_str.split("-")` with parsing through `RepeatUnit.from_str`:
```python
chain = [RepeatUnit.from_str(s) for s in chain_str.split("-")]
```
- Remove `pure_symbol = symbol.rstrip("m")` validation — use `unit.symbol` instead.
- Return `list[list[RepeatUnit]]` instead of `list[list[str]]`.

- [ ] **Step 2: Update simulation_statistics.py**

Replace all "m" suffix operations:
- `get_repeat_lengths()`: `pure_symbol = symbol.rstrip("m")` → accept `list[RepeatUnit]`, use `unit.symbol`
- `count_repeat_types()`: same pattern
- `get_mutation_details()`: `symbol.endswith("m")` → `unit.mutated`; `symbol.rstrip("m")` → `unit.symbol`
- `generate_haplotype_stats()`: `for seq, chain in simulation_results:` → `for hr in simulation_results:` using `hr.sequence`, `hr.chain`
- `sum(1 for r in chain if r.endswith("m"))` → `sum(1 for ru in hr.chain if ru.mutated)`

- [ ] **Step 3: Update translate.py**

Update `predict_orfs_in_haplotypes()` signature:
- `results: list[tuple[str, list[str]]]` → `results: list[HaplotypeResult]`
- Update unpacking accordingly

- [ ] **Step 4: Update snp_integrator.py**

Update `get_vntr_boundaries()` signature:
- `simulation_results: list[tuple[str, list[str]]]` → `simulation_results: list[HaplotypeResult]`

- [ ] **Step 5: Update validation.py**

Update `MutationTargets` usage to `list[MutationTarget]`.

- [ ] **Step 6: Update read_simulator/source_tracking.py**

Replace:
- `lookup_key = symbol.rstrip("m")` → `unit.symbol`
- `symbol.endswith("m")` → `unit.mutated`

- [ ] **Step 7: Run tests**

Run: `pytest --tb=short -q`
Expected: All PASS

- [ ] **Step 8: Commit**

```bash
git add muc_one_up/io.py muc_one_up/simulation_statistics.py muc_one_up/translate.py muc_one_up/snp_integrator.py muc_one_up/validation.py muc_one_up/read_simulator/source_tracking.py
git commit -m "refactor: migrate peripheral modules to typed domain models"
```

---

## Task 8: Clean Up type_defs.py and Update RNG Threading Tests

**Files:**
- Modify: `muc_one_up/type_defs.py`
- Modify: `tests/test_rng_threading.py`

- [ ] **Step 1: Remove deprecated aliases from type_defs.py**

Remove or mark as deprecated:
- `Haplotype = tuple[HaplotypeName, RepeatChain]` — replaced by `HaplotypeResult`
- `HaplotypeList = list[Haplotype]` — replaced by `list[HaplotypeResult]`
- `MutationTargets = list[tuple[int, int]]` — replaced by `list[MutationTarget]`

Keep `RepeatChain = list[str]` as it's used by `io.py` for raw parsing.

- [ ] **Step 2: Update test_rng_threading.py**

Update `TestSimulateRNG.test_simulate_diploid_reproducible_with_rng` to use `HaplotypeResult`:
- Replace `r1[0][1] == r2[0][1]` with `r1[0].chain == r2[0].chain`

- [ ] **Step 3: Run tests**

Run: `pytest --tb=short -q`
Expected: All PASS

- [ ] **Step 4: Commit**

```bash
git add muc_one_up/type_defs.py tests/test_rng_threading.py
git commit -m "refactor: clean up deprecated type aliases in type_defs.py"
```

---

## Task 9: Full Verification and CI Readiness

**Files:** None (verification only)

- [ ] **Step 1: Run full test suite**

Run: `pytest --tb=short -q`
Expected: All tests pass.

- [ ] **Step 2: Run ruff linter and formatter**

Run: `ruff check muc_one_up/ tests/ && ruff format --check muc_one_up/ tests/`
Expected: Clean.

- [ ] **Step 3: Run mypy**

Run: `mypy muc_one_up/`
Expected: No new errors.

- [ ] **Step 4: Verify success criteria**

```bash
# No .rstrip("m"), .replace("m", ""), + "m", or .endswith("m") in domain modules
grep -rn 'rstrip("m")\|replace("m"\|+ "m"\|endswith("m")' muc_one_up/simulate.py muc_one_up/mutate.py muc_one_up/assembly.py muc_one_up/simulation_statistics.py muc_one_up/io.py
# Expected: zero hits

# No raw tuple[str, list[str]] in domain function signatures
grep -rn "tuple\[str, list\[str\]\]" muc_one_up/
# Expected: zero hits (all migrated to HaplotypeResult)

# Seed reproducibility preserved
pytest tests/test_simulate.py::TestSimulateDiploidAdvanced::test_seed_ensures_reproducibility -v --no-cov
# Expected: PASS
```

- [ ] **Step 5: Commit any formatting fixes**

```bash
git add -u
git commit -m "style: fix formatting after Phase 3B changes"
```

- [ ] **Step 6: Final test run**

Run: `pytest --tb=short -q`
Expected: All pass. Phase 3B is complete.
