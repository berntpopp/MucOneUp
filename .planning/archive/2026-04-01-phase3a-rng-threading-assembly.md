# Phase 3A: Thread Explicit RNG & Centralize Assembly — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace global `random` module usage with explicit `rng: random.Random` instances in all domain functions, and consolidate the two redundant assembly functions into a single `assemble_sequence()`.

**Architecture:** Add `rng: random.Random` parameter to every domain function that uses randomness. The CLI entry point creates `rng = random.Random(seed)` and threads it through. Merge `assemble_haplotype_from_chain()` (simulate.py:85) and `rebuild_haplotype_sequence()` (mutate.py:211) into a single `assemble_sequence()` in a new `assembly.py` module. Both callers switch to the shared function.

**Tech Stack:** Python 3.10+, random.Random, dataclasses, pytest

**Spec:** `.planning/specs/2026-04-01-codebase-refactoring-design.md` (sections 3.3, 3.4)

---

## File Structure

### New Files
| File | Responsibility |
|------|---------------|
| `muc_one_up/assembly.py` | Single `assemble_sequence()` function for chain → sequence conversion |
| `tests/test_assembly.py` | Tests for centralized assembly |

### Modified Files
| File | Changes |
|------|---------|
| `muc_one_up/simulate.py` | Accept `rng` parameter; call `assemble_sequence()` instead of inline assembly |
| `muc_one_up/mutate.py` | Accept `rng` parameter; call `assemble_sequence()` instead of `rebuild_haplotype_sequence()` |
| `muc_one_up/distribution.py` | Accept `rng` parameter instead of global `random` |
| `muc_one_up/probabilities.py` | Accept `rng` parameter instead of global `random` |
| `muc_one_up/cli/haplotypes.py` | Create `rng = random.Random(seed)` and pass to simulate functions |

---

## Task 1: Create Centralized assembly.py

**Files:**
- Create: `muc_one_up/assembly.py`
- Create: `tests/test_assembly.py`

- [ ] **Step 1: Write failing tests for assemble_sequence**

```python
# tests/test_assembly.py
"""Tests for centralized sequence assembly."""

from __future__ import annotations

import pytest

from muc_one_up.assembly import assemble_sequence


class TestAssembleSequence:
    """Tests for assemble_sequence()."""

    @pytest.fixture()
    def config(self):
        return {
            "repeats": {
                "1": "GCCACCACCTCCAACTCCT",
                "2": "GCCACCACC",
                "X": "GCCACCACCTCCAACTCCTGCCACC",
                "9": "TCTAGGACCTAGCTCCT",
            },
            "constants": {
                "hg38": {
                    "left": "ATGGCCCC",
                    "right": "CAATGGTG",
                }
            },
            "reference_assembly": "hg38",
        }

    def test_simple_chain(self, config):
        chain = ["1", "2", "9"]
        seq = assemble_sequence(chain, config)
        left = config["constants"]["hg38"]["left"]
        right = config["constants"]["hg38"]["right"]
        r1 = config["repeats"]["1"]
        r2 = config["repeats"]["2"]
        r9 = config["repeats"]["9"]
        assert seq == left + r1 + r2 + r9 + right

    def test_mutation_marker_stripped(self, config):
        """Chains with 'm' suffix should look up the base symbol."""
        chain = ["1", "Xm", "9"]
        seq = assemble_sequence(chain, config)
        left = config["constants"]["hg38"]["left"]
        right = config["constants"]["hg38"]["right"]
        r1 = config["repeats"]["1"]
        rX = config["repeats"]["X"]
        r9 = config["repeats"]["9"]
        assert seq == left + r1 + rX + r9 + right

    def test_empty_chain(self, config):
        seq = assemble_sequence([], config)
        left = config["constants"]["hg38"]["left"]
        right = config["constants"]["hg38"]["right"]
        assert seq == left + right

    def test_chain_not_ending_with_9_omits_right_constant(self, config):
        """If chain doesn't end with 9, right constant is omitted."""
        chain = ["1", "2"]
        seq = assemble_sequence(chain, config)
        left = config["constants"]["hg38"]["left"]
        r1 = config["repeats"]["1"]
        r2 = config["repeats"]["2"]
        assert seq == left + r1 + r2

    def test_unknown_symbol_raises(self, config):
        with pytest.raises(KeyError):
            assemble_sequence(["UNKNOWN"], config)

    def test_hg19_assembly(self):
        config = {
            "repeats": {"1": "ACGT"},
            "constants": {
                "hg19": {"left": "LEFT", "right": "RIGHT"},
                "hg38": {"left": "OTHER", "right": "OTHER"},
            },
            "reference_assembly": "hg19",
        }
        # Chain doesn't end with "9", so no right constant
        seq = assemble_sequence(["1"], config)
        assert seq == "LEFT" + "ACGT"

    def test_single_repeat_ending_with_9(self, config):
        chain = ["9"]
        seq = assemble_sequence(chain, config)
        left = config["constants"]["hg38"]["left"]
        right = config["constants"]["hg38"]["right"]
        r9 = config["repeats"]["9"]
        assert seq == left + r9 + right
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_assembly.py -v`
Expected: FAIL with ImportError

- [ ] **Step 3: Implement assemble_sequence**

Create `muc_one_up/assembly.py` by consolidating `assemble_haplotype_from_chain()` (simulate.py:85-124) and `rebuild_haplotype_sequence()` (mutate.py:211-236). Both do the same thing: strip "m" suffix, look up repeat sequences, concatenate with constants.

```python
# muc_one_up/assembly.py
"""Centralized sequence assembly from repeat chains.

Single source of truth for chain → DNA sequence conversion.
Replaces assemble_haplotype_from_chain() in simulate.py and
rebuild_haplotype_sequence() in mutate.py.
"""

from __future__ import annotations

import logging

from .type_defs import ConfigDict, DNASequence, RepeatChain

logger = logging.getLogger(__name__)


def assemble_sequence(chain: RepeatChain, config: ConfigDict) -> DNASequence:
    """Assemble a DNA sequence from a repeat chain and flanking constants.

    Concatenates: left_constant + repeat_units + right_constant.
    The right constant is only appended if the chain ends with repeat "9"
    (the canonical terminal repeat in MUC1 VNTR).

    Mutation markers (trailing 'm' on symbols like 'Xm') are stripped
    before looking up repeat sequences, so mutated and unmutated repeats
    use the same base sequence.

    Args:
        chain: List of repeat symbols (e.g., ["1", "2", "Xm", "7", "8", "9"]).
        config: Configuration dict with 'repeats', 'constants', and
                optionally 'reference_assembly' keys.

    Returns:
        Assembled DNA sequence string.

    Raises:
        KeyError: If a repeat symbol (after stripping 'm') is not in config.
    """
    repeats_dict = config["repeats"]
    ref_assembly = config.get("reference_assembly", "hg38")
    constants = config.get("constants", {}).get(ref_assembly, {})
    left_const = constants.get("left", "")
    right_const = constants.get("right", "")

    # Build repeat region
    parts: list[str] = [left_const]
    for symbol in chain:
        base_symbol = symbol.rstrip("m")
        parts.append(repeats_dict[base_symbol])

    # Right constant only if chain ends with canonical terminal repeat "9"
    if chain and chain[-1].replace("m", "") == "9":
        parts.append(right_const)

    return "".join(parts)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_assembly.py -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/assembly.py tests/test_assembly.py
git commit -m "feat: add centralized assemble_sequence replacing dual assembly functions"
```

---

## Task 2: Migrate simulate.py to Use assemble_sequence

**Files:**
- Modify: `muc_one_up/simulate.py:85-124` (replace `assemble_haplotype_from_chain`)

- [ ] **Step 1: Write test verifying assembly delegation**

Add to `tests/test_assembly.py`:

```python
def test_simulate_uses_centralized_assembly():
    """simulate.py should delegate to assembly.assemble_sequence."""
    import ast
    import inspect

    import muc_one_up.simulate as mod

    source = inspect.getsource(mod)
    assert "from .assembly import assemble_sequence" in source or \
           "from muc_one_up.assembly import assemble_sequence" in source, (
        "simulate.py should import assemble_sequence from assembly module"
    )
```

- [ ] **Step 2: Replace assemble_haplotype_from_chain with a thin wrapper**

In `muc_one_up/simulate.py`:

1. Add import at top: `from .assembly import assemble_sequence`

2. Replace the body of `assemble_haplotype_from_chain()` (lines 85-124) with a delegation to `assemble_sequence`:

```python
def assemble_haplotype_from_chain(chain: RepeatChain, config: ConfigDict) -> DNASequence:
    """Assemble a haplotype sequence from a repeat chain.

    Delegates to assembly.assemble_sequence().
    Kept as a thin wrapper for backward compatibility.
    """
    return assemble_sequence(chain, config)
```

3. Update `simulate_single_haplotype()` (line ~302) which calls `assemble_haplotype_from_chain(chain, config)` — this still works since the wrapper delegates.

- [ ] **Step 3: Run tests**

Run: `pytest tests/test_simulate.py tests/test_assembly.py --tb=short -q --no-cov`
Expected: All PASS (behavior unchanged)

- [ ] **Step 4: Commit**

```bash
git add muc_one_up/simulate.py tests/test_assembly.py
git commit -m "refactor: delegate simulate.py assembly to centralized assemble_sequence"
```

---

## Task 3: Migrate mutate.py to Use assemble_sequence

**Files:**
- Modify: `muc_one_up/mutate.py:211-236` (replace `rebuild_haplotype_sequence`)

- [ ] **Step 1: Write test verifying assembly delegation**

Add to `tests/test_assembly.py`:

```python
def test_mutate_uses_centralized_assembly():
    """mutate.py should delegate to assembly.assemble_sequence."""
    import inspect

    import muc_one_up.mutate as mod

    source = inspect.getsource(mod)
    assert "from .assembly import assemble_sequence" in source or \
           "from muc_one_up.assembly import assemble_sequence" in source, (
        "mutate.py should import assemble_sequence from assembly module"
    )
```

- [ ] **Step 2: Replace rebuild_haplotype_sequence with a thin wrapper**

In `muc_one_up/mutate.py`:

1. Add import at top: `from .assembly import assemble_sequence`

2. Replace the body of `rebuild_haplotype_sequence()` (lines 211-236) with:

```python
def rebuild_haplotype_sequence(chain: RepeatChain, config: ConfigDict) -> DNASequence:
    """Rebuild a haplotype sequence after mutation.

    Delegates to assembly.assemble_sequence().
    Kept as a thin wrapper for backward compatibility.
    """
    return assemble_sequence(chain, config)
```

- [ ] **Step 3: Run tests**

Run: `pytest tests/test_mutate.py tests/test_assembly.py --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 4: Run the seed reproducibility test**

Run: `pytest tests/test_simulate.py::TestSimulateDiploidAdvanced::test_seed_ensures_reproducibility -v --no-cov`
Expected: PASS (byte-identical output preserved)

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/mutate.py tests/test_assembly.py
git commit -m "refactor: delegate mutate.py assembly to centralized assemble_sequence"
```

---

## Task 4: Thread RNG Through distribution.py and probabilities.py

**Files:**
- Modify: `muc_one_up/distribution.py:9-47`
- Modify: `muc_one_up/probabilities.py:9-47`

- [ ] **Step 1: Write failing tests for RNG parameter**

```python
# tests/test_rng_threading.py
"""Tests for explicit RNG threading through domain functions."""

from __future__ import annotations

import random

from muc_one_up.distribution import sample_repeat_count
from muc_one_up.probabilities import pick_next_repeat


class TestDistributionRNG:
    def test_sample_repeat_count_accepts_rng(self):
        """sample_repeat_count should accept an rng parameter."""
        rng = random.Random(42)
        length_model = {
            "distribution": "normal",
            "min_repeats": 20,
            "max_repeats": 60,
            "mean_repeats": 40,
        }
        result = sample_repeat_count(length_model, rng=rng)
        assert 20 <= result <= 60

    def test_sample_repeat_count_reproducible_with_rng(self):
        """Same RNG seed should produce same results."""
        length_model = {
            "distribution": "normal",
            "min_repeats": 20,
            "max_repeats": 60,
            "mean_repeats": 40,
        }
        results_a = [sample_repeat_count(length_model, rng=random.Random(42)) for _ in range(5)]
        results_b = [sample_repeat_count(length_model, rng=random.Random(42)) for _ in range(5)]
        assert results_a == results_b

    def test_sample_repeat_count_defaults_to_global_random(self):
        """Without rng parameter, should still work (backward compat)."""
        length_model = {
            "distribution": "normal",
            "min_repeats": 20,
            "max_repeats": 60,
            "mean_repeats": 40,
        }
        result = sample_repeat_count(length_model)
        assert 20 <= result <= 60


class TestProbabilitiesRNG:
    def test_pick_next_repeat_accepts_rng(self):
        """pick_next_repeat should accept an rng parameter."""
        rng = random.Random(42)
        probs = {"1": {"2": 0.5, "3": 0.5}, "2": {}, "3": {}}
        result = pick_next_repeat(probs, "1", rng=rng)
        assert result in ("2", "3")

    def test_pick_next_repeat_reproducible_with_rng(self):
        """Same RNG seed should produce same sequence."""
        probs = {"1": {"2": 0.5, "3": 0.5}, "2": {}, "3": {}}
        results_a = [pick_next_repeat(probs, "1", rng=random.Random(42)) for _ in range(10)]
        results_b = [pick_next_repeat(probs, "1", rng=random.Random(42)) for _ in range(10)]
        assert results_a == results_b

    def test_pick_next_repeat_defaults_to_global_random(self):
        """Without rng parameter, should still work (backward compat)."""
        probs = {"1": {"2": 1.0}, "2": {}}
        result = pick_next_repeat(probs, "1")
        assert result == "2"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_rng_threading.py -v`
Expected: FAIL (functions don't accept `rng` parameter yet)

- [ ] **Step 3: Add rng parameter to sample_repeat_count**

In `muc_one_up/distribution.py`, modify `sample_repeat_count()`:

```python
import random as _random_module

def sample_repeat_count(
    length_model: LengthModelDict,
    rng: _random_module.Random | None = None,
) -> int:
```

Inside the function, replace `random.gauss(...)` with:

```python
    _rng = rng if rng is not None else _random_module
    # ... then use _rng.gauss(...) instead of random.gauss(...)
```

The `_rng` variable works because both `random.Random` instances and the `random` module itself expose the same `.gauss()`, `.randint()`, etc. API.

- [ ] **Step 4: Add rng parameter to pick_next_repeat**

In `muc_one_up/probabilities.py`, modify `pick_next_repeat()`:

```python
import random as _random_module

def pick_next_repeat(
    probabilities: ProbabilitiesDict,
    current_symbol: str,
    force_end: bool = False,
    rng: _random_module.Random | None = None,
) -> str:
```

Replace `random.choices(...)` with `_rng.choices(...)` where `_rng = rng if rng is not None else _random_module`.

- [ ] **Step 5: Run tests**

Run: `pytest tests/test_rng_threading.py tests/test_distribution.py tests/test_probabilities.py --tb=short -q --no-cov`
Expected: All PASS (new tests pass, existing tests pass via backward-compat default)

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/distribution.py muc_one_up/probabilities.py tests/test_rng_threading.py
git commit -m "feat: add explicit rng parameter to distribution and probabilities functions"
```

---

## Task 5: Thread RNG Through simulate.py

**Files:**
- Modify: `muc_one_up/simulate.py`

- [ ] **Step 1: Write failing test for RNG in simulate**

Add to `tests/test_rng_threading.py`:

```python
from muc_one_up.simulate import simulate_diploid


class TestSimulateRNG:
    def test_simulate_diploid_accepts_rng(self):
        """simulate_diploid should accept an rng parameter."""
        import inspect
        sig = inspect.signature(simulate_diploid)
        assert "rng" in sig.parameters

    def test_simulate_diploid_reproducible_with_rng(self):
        """Same RNG produces identical haplotypes."""
        from tests.conftest import minimal_config
        # Use a real config fixture
        config = {
            "repeats": {
                "1": "GCCACCACCTCCAACTCCT",
                "2": "GCCACCACC",
                "3": "TCCAACTCCT",
                "6": "TCCGACTCCT",
                "6p": "TCCAACTCCT",
                "7": "TCTAGGACCTCCTAGCTCCCAGAACTCCT",
                "8": "TCTAGGACCTCCTAGCTCCTAGCTCCT",
                "9": "TCTAGGACCTAGCTCCT",
            },
            "constants": {
                "hg38": {
                    "left": "ATGGCCCCATCTCTCACCGTCTCGGTCATCTCCTTGATG",
                    "right": "CAATGGTGTCTTGGGTAGCTTCGTCACGGTTTTCCAG",
                }
            },
            "probabilities": {
                "1": {"2": 1.0},
                "2": {"3": 1.0},
                "3": {"6": 0.5, "6p": 0.5},
                "6": {"7": 1.0},
                "6p": {"7": 1.0},
                "7": {"8": 1.0},
                "8": {"9": 1.0},
                "9": {},
            },
            "length_model": {
                "distribution": "normal",
                "min_repeats": 20,
                "max_repeats": 40,
                "mean_repeats": 30,
            },
            "reference_assembly": "hg38",
        }
        rng1 = random.Random(42)
        rng2 = random.Random(42)
        r1 = simulate_diploid(config, seed=None, rng=rng1)
        r2 = simulate_diploid(config, seed=None, rng=rng2)
        assert r1[0][1] == r2[0][1]  # same chains
        assert r1[1][1] == r2[1][1]
```

- [ ] **Step 2: Add rng parameter to simulate functions**

In `muc_one_up/simulate.py`:

1. Add `rng` parameter to `pick_next_symbol_no_end()`:
```python
def pick_next_symbol_no_end(
    probabilities: ProbabilitiesDict,
    current_symbol: str,
    rng: _random_module.Random | None = None,
) -> str | None:
```
Replace `random.choices(...)` with `_rng.choices(...)`.

2. Add `rng` parameter to `simulate_single_haplotype()`:
```python
def simulate_single_haplotype(
    config: ConfigDict,
    target_length: int,
    min_length: int = 10,
    rng: _random_module.Random | None = None,
) -> Haplotype:
```
Replace `random.choice(["6", "6p"])` with `_rng.choice(["6", "6p"])`.
Pass `rng=rng` to `pick_next_symbol_no_end()`.

3. Add `rng` parameter to `simulate_diploid()`:
```python
def simulate_diploid(
    config: ConfigDict,
    num_haplotypes: int = 2,
    fixed_lengths: list[int] | None = None,
    seed: int | None = None,
    rng: _random_module.Random | None = None,
) -> HaplotypeList:
```
If `rng` is None and `seed` is not None, create `rng = random.Random(seed)`.
If `rng` is None and `seed` is None, use global `random` (backward compat).
Pass `rng=rng` to `sample_repeat_count()` and `simulate_single_haplotype()`.

Remove `random.seed(seed)` call — the RNG instance handles seeding.

- [ ] **Step 3: Run tests**

Run: `pytest tests/test_rng_threading.py tests/test_simulate.py --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 4: Run seed reproducibility test**

Run: `pytest tests/test_simulate.py::TestSimulateDiploidAdvanced::test_seed_ensures_reproducibility -v --no-cov`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/simulate.py tests/test_rng_threading.py
git commit -m "feat: thread explicit rng through simulate.py domain functions"
```

---

## Task 6: Thread RNG Through mutate.py

**Files:**
- Modify: `muc_one_up/mutate.py:85-209`

- [ ] **Step 1: Add rng parameter to apply_mutations**

In `muc_one_up/mutate.py`, modify `apply_mutations()`:

```python
import random as _random_module

def apply_mutations(
    config: ConfigDict,
    results: HaplotypeList,
    mutation_name: MutationName,
    targets: MutationTargets,
    rng: _random_module.Random | None = None,
) -> tuple[HaplotypeList, dict[int, list[tuple[int, str]]]]:
```

Replace `random.choice(list(allowed_repeats))` (line 175) with `_rng.choice(list(allowed_repeats))` where `_rng = rng if rng is not None else _random_module`.

- [ ] **Step 2: Run tests**

Run: `pytest tests/test_mutate.py --tb=short -q --no-cov`
Expected: All PASS (backward compat via default None)

- [ ] **Step 3: Commit**

```bash
git add muc_one_up/mutate.py
git commit -m "feat: thread explicit rng through mutate.py apply_mutations"
```

---

## Task 7: Update CLI Entry Point to Create and Pass RNG

**Files:**
- Modify: `muc_one_up/cli/haplotypes.py`

- [ ] **Step 1: Update generate_haplotypes to create and pass RNG**

In `muc_one_up/cli/haplotypes.py`, modify `generate_haplotypes()`:

```python
import random

def generate_haplotypes(args, config, fixed_conf, predefined_chains):
    """Generate haplotypes using typed RNG instance."""
    # Create explicit RNG from seed
    rng = random.Random(args.seed) if args.seed is not None else None

    if fixed_conf == "from_structure":
        results = simulate_from_chains(
            predefined_chains=predefined_chains, config=config
        )
    else:
        results = simulate_diploid(
            config=config,
            num_haplotypes=args.num_haplotypes,
            fixed_lengths=fixed_conf,
            seed=args.seed,
            rng=rng,
        )
    return results
```

Note: `simulate_from_chains` does not use randomness (chains are predefined), so no RNG needed there.

- [ ] **Step 2: Run tests**

Run: `pytest tests/test_click_cli.py tests/test_simulate.py --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 3: Commit**

```bash
git add muc_one_up/cli/haplotypes.py
git commit -m "feat: create explicit RNG at CLI entry point and pass through pipeline"
```

---

## Task 8: Verify Full Suite and CI Readiness

**Files:** None (verification only)

- [ ] **Step 1: Run full test suite**

Run: `pytest --tb=short -q`
Expected: All tests pass.

- [ ] **Step 2: Run ruff linter and formatter**

Run: `ruff check muc_one_up/ tests/ && ruff format --check muc_one_up/ tests/`
Expected: Clean. Fix any issues.

- [ ] **Step 3: Run mypy**

Run: `mypy muc_one_up/`
Expected: No errors.

- [ ] **Step 4: Verify success criteria**

```bash
# No global random usage in core domain files
grep -rn "random\.\(choice\|choices\|seed\|gauss\|randint\)" muc_one_up/simulate.py muc_one_up/mutate.py muc_one_up/distribution.py muc_one_up/probabilities.py | grep -v "_rng\." | grep -v "_random_module"
# Expected: zero hits (all go through _rng)

# Assembly functions delegate to assembly.py
grep -rn "assemble_sequence" muc_one_up/simulate.py muc_one_up/mutate.py muc_one_up/assembly.py
# Expected: imports in simulate.py and mutate.py, definition in assembly.py

# Seed reproducibility preserved
pytest tests/test_simulate.py::TestSimulateDiploidAdvanced::test_seed_ensures_reproducibility -v --no-cov
# Expected: PASS
```

- [ ] **Step 5: Commit any formatting fixes**

```bash
git add -u
git commit -m "style: fix formatting after Phase 3A changes"
```

- [ ] **Step 6: Final test run**

Run: `pytest --tb=short -q`
Expected: All pass. Phase 3A is complete.
