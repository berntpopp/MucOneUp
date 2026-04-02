# Phase 4B: Assembly Context, Temp-File Fix & Alignment Dedup — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate duplicated assembly resolution across pipeline files, fix stale temp-file paths in diploid handler, and consolidate overlapping alignment logic between minimap2 and nanosim wrappers.

**Architecture:** Introduce `AssemblyContext` dataclass holding resolved assembly name, constants, and reference paths — constructed once at pipeline entry. Fix `DiploidSimulationResult` to set temp paths to `None` when `keep_intermediate=False`. Audit minimap2/nanosim for shared alignment logic and extract into a common function.

**Tech Stack:** Python 3.10+, dataclasses, pytest

**Spec:** `.planning/specs/2026-04-01-codebase-refactoring-design.md` (sections 4.1, 4.4, 4.5)

**Prerequisite:** Phase 4A (execution abstraction) should be completed first, but this plan can be executed independently if needed.

---

## File Structure

### New Files
| File | Responsibility |
|------|---------------|
| `muc_one_up/read_simulator/assembly_context.py` | `AssemblyContext` dataclass and factory |

### Modified Files
| File | Changes |
|------|---------|
| `muc_one_up/read_simulator/pipeline.py` | Use `AssemblyContext` instead of re-deriving assembly 3 times |
| `muc_one_up/read_simulator/ont_pipeline.py` | Use `AssemblyContext` instead of local assembly derivation |
| `muc_one_up/read_simulator/utils/diploid_handler.py` | Set stale temp paths to `None` when `keep_intermediate=False` |
| `muc_one_up/read_simulator/ont_pipeline.py` | Handle `None` hap reference fields |

### New Test Files
| File | Responsibility |
|------|---------------|
| `tests/read_simulator/test_assembly_context.py` | Tests for AssemblyContext construction and resolution |

---

## Task 1: Create AssemblyContext Dataclass

**Files:**
- Create: `muc_one_up/read_simulator/assembly_context.py`
- Create: `tests/read_simulator/test_assembly_context.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/read_simulator/test_assembly_context.py
"""Tests for AssemblyContext."""

from __future__ import annotations

import pytest

from muc_one_up.read_simulator.assembly_context import AssemblyContext


class TestAssemblyContext:
    """Tests for AssemblyContext dataclass."""

    def test_create_from_config(self):
        config = {
            "reference_assembly": "hg38",
            "constants": {
                "hg38": {"left": "ATGGCCCC", "right": "CAATGGTG"}
            },
        }
        rs_config = {
            "sample_bam_hg38": "/path/to/sample.bam",
            "vntr_region_hg38": "chr1:155000-156000",
            "human_reference_hg38": "/path/to/hg38.fa",
        }
        ctx = AssemblyContext.from_configs(config, rs_config)
        assert ctx.assembly_name == "hg38"
        assert ctx.left_constant == "ATGGCCCC"
        assert ctx.right_constant == "CAATGGTG"
        assert ctx.sample_bam == "/path/to/sample.bam"
        assert ctx.vntr_region == "chr1:155000-156000"
        assert ctx.human_reference == "/path/to/hg38.fa"

    def test_defaults_to_hg38(self):
        config = {
            "constants": {
                "hg38": {"left": "L", "right": "R"}
            },
        }
        rs_config = {}
        ctx = AssemblyContext.from_configs(config, rs_config)
        assert ctx.assembly_name == "hg38"

    def test_hg19_assembly(self):
        config = {
            "reference_assembly": "hg19",
            "constants": {
                "hg19": {"left": "L19", "right": "R19"}
            },
        }
        rs_config = {
            "sample_bam_hg19": "/bam19",
            "vntr_region_hg19": "region19",
            "human_reference_hg19": "/ref19",
        }
        ctx = AssemblyContext.from_configs(config, rs_config)
        assert ctx.assembly_name == "hg19"
        assert ctx.left_constant == "L19"
        assert ctx.sample_bam == "/bam19"

    def test_missing_assembly_key_returns_none(self):
        """Missing assembly-specific keys should return None, not raise."""
        config = {
            "reference_assembly": "hg38",
            "constants": {
                "hg38": {"left": "L", "right": "R"}
            },
        }
        rs_config = {}  # No assembly-specific keys
        ctx = AssemblyContext.from_configs(config, rs_config)
        assert ctx.sample_bam is None
        assert ctx.vntr_region is None
        assert ctx.human_reference is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/read_simulator/test_assembly_context.py -v --no-cov`
Expected: FAIL with ImportError

- [ ] **Step 3: Implement AssemblyContext**

```python
# muc_one_up/read_simulator/assembly_context.py
"""Centralized assembly resolution for read simulation pipelines.

Resolves assembly name, constants, and reference paths once at pipeline
entry. Eliminates repeated config.get("reference_assembly", "hg38")
calls scattered across pipeline stages.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True, slots=True)
class AssemblyContext:
    """Resolved assembly configuration for a pipeline run.

    Attributes:
        assembly_name: Reference assembly identifier (e.g., "hg38", "hg19").
        left_constant: Left flanking constant sequence.
        right_constant: Right flanking constant sequence.
        sample_bam: Path to sample BAM for the assembly (optional).
        vntr_region: VNTR region string for the assembly (optional).
        human_reference: Path to human reference FASTA (optional).
    """

    assembly_name: str
    left_constant: str
    right_constant: str
    sample_bam: str | None = None
    vntr_region: str | None = None
    human_reference: str | None = None

    @classmethod
    def from_configs(
        cls,
        config: dict[str, Any],
        rs_config: dict[str, Any] | None = None,
    ) -> AssemblyContext:
        """Construct from simulation config and read-simulator config.

        Args:
            config: Main simulation config with 'reference_assembly' and 'constants'.
            rs_config: Read-simulator config with assembly-keyed paths (optional).

        Returns:
            Resolved AssemblyContext.
        """
        rs_config = rs_config or {}
        assembly = config.get("reference_assembly", "hg38")
        constants = config.get("constants", {}).get(assembly, {})

        return cls(
            assembly_name=assembly,
            left_constant=constants.get("left", ""),
            right_constant=constants.get("right", ""),
            sample_bam=rs_config.get(f"sample_bam_{assembly}"),
            vntr_region=rs_config.get(f"vntr_region_{assembly}"),
            human_reference=rs_config.get(f"human_reference_{assembly}"),
        )
```

- [ ] **Step 4: Run tests**

Run: `uv run pytest tests/read_simulator/test_assembly_context.py -v --no-cov`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/read_simulator/assembly_context.py tests/read_simulator/test_assembly_context.py
git commit -m "feat: add AssemblyContext dataclass for centralized assembly resolution"
```

---

## Task 2: Migrate Pipeline Files to AssemblyContext

**Files:**
- Modify: `muc_one_up/read_simulator/pipeline.py`
- Modify: `muc_one_up/read_simulator/ont_pipeline.py`

- [ ] **Step 1: Read both pipeline files**

Identify all locations where `config.get("reference_assembly", "hg38")` is called and assembly-specific config keys are constructed.

- [ ] **Step 2: Update pipeline.py**

1. Import `AssemblyContext` at top
2. Construct `AssemblyContext` once at the top of `simulate_reads_pipeline()`:
```python
from ..assembly_context import AssemblyContext
assembly_ctx = AssemblyContext.from_configs(config, rs_config)
```
3. Replace all `reference_assembly = rs_config.get("reference_assembly", "hg38")` with `assembly_ctx.assembly_name`
4. Replace all `rs_config.get(f"sample_bam_{assembly}")` with `assembly_ctx.sample_bam`
5. Replace all `rs_config.get(f"vntr_region_{assembly}")` with `assembly_ctx.vntr_region`
6. Replace all `rs_config.get(f"human_reference_{assembly}")` with `assembly_ctx.human_reference`

- [ ] **Step 3: Update ont_pipeline.py**

Same pattern — construct `AssemblyContext` once at entry, use throughout.

- [ ] **Step 4: Run tests**

Run: `uv run pytest tests/read_simulator/ --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/read_simulator/pipeline.py muc_one_up/read_simulator/ont_pipeline.py
git commit -m "refactor: migrate pipeline files to use AssemblyContext"
```

---

## Task 3: Fix Stale Temp-File Paths in DiploidSimulationResult

**Files:**
- Modify: `muc_one_up/read_simulator/utils/diploid_handler.py`
- Modify: `muc_one_up/read_simulator/ont_pipeline.py` (handle None fields)

- [ ] **Step 1: Write failing test**

Add to existing diploid handler tests (or create new test):

```python
def test_stale_paths_none_when_not_keeping_intermediate():
    """When keep_intermediate=False, temp paths should be None."""
    # This test verifies the fix for stale temp-file paths
    # After run_split_simulation with keep_intermediate=False,
    # hap1_reference and hap2_reference should be None
    # (they pointed to deleted temp files before this fix)
```

The exact test depends on the existing test infrastructure. Read the current tests first.

- [ ] **Step 2: Fix DiploidSimulationResult construction**

In `diploid_handler.py`, in the `run_split_simulation()` function, when `keep_intermediate=False`:

Change the `DiploidSimulationResult` construction to set reference paths to `None`:

```python
if not keep_intermediate:
    # Temp paths are invalid after context manager exits
    return DiploidSimulationResult(
        merged_fastq=output_fastq,
        hap1_fastq=None,  # was: temp path (now deleted)
        hap2_fastq=None,  # was: temp path (now deleted)
        hap1_reference=None,  # was: temp path (now deleted)
        hap2_reference=None,  # was: temp path (now deleted)
        reads_hap1=reads_hap1,
        reads_hap2=reads_hap2,
    )
```

**Note:** This requires changing `DiploidSimulationResult` fields from `str` to `str | None`. Since it's a `NamedTuple`, the type hints need updating:

```python
class DiploidSimulationResult(NamedTuple):
    merged_fastq: str
    hap1_fastq: str | None
    hap2_fastq: str | None
    hap1_reference: str | None
    hap2_reference: str | None
    reads_hap1: int
    reads_hap2: int
```

- [ ] **Step 3: Update ont_pipeline.py to handle None fields**

Check any code that accesses `result.hap1_reference` or `result.hap2_reference` and add None checks:

```python
if result.hap1_reference is not None:
    # use reference path
```

- [ ] **Step 4: Run tests**

Run: `uv run pytest tests/read_simulator/ tests/test_diploid_handler.py --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/read_simulator/utils/diploid_handler.py muc_one_up/read_simulator/ont_pipeline.py
git commit -m "fix: set stale temp-file paths to None in DiploidSimulationResult"
```

---

## Task 4: Audit and Consolidate Alignment Logic

**Files:**
- Modify: `muc_one_up/read_simulator/wrappers/minimap2_wrapper.py`
- Modify: `muc_one_up/read_simulator/wrappers/nanosim_wrapper.py`

- [ ] **Step 1: Read both files and identify overlap**

Read `minimap2_wrapper.py` and `nanosim_wrapper.py` completely. Identify any shared alignment patterns (e.g., both call minimap2 or both do SAM→BAM conversion).

- [ ] **Step 2: Extract shared logic if significant overlap exists**

If the overlap is minimal (just a SAM→BAM subprocess call), the consolidation may not be worth the complexity. Document the finding and move on.

If significant overlap exists, extract a shared function (e.g., `align_and_convert()`) that both wrappers call.

- [ ] **Step 3: Run tests**

Run: `uv run pytest tests/read_simulator/test_minimap2_wrapper.py tests/read_simulator/test_nanosim_wrapper.py --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 4: Commit (if changes were made)**

```bash
git add muc_one_up/read_simulator/wrappers/
git commit -m "refactor: consolidate shared alignment logic between minimap2 and nanosim wrappers"
```

---

## Task 5: Full Verification

**Files:** None (verification only)

- [ ] **Step 1: Run full test suite**

Run: `uv run pytest --tb=short -q`
Expected: All tests pass.

- [ ] **Step 2: Run ruff and mypy**

Run: `uv run ruff check muc_one_up/ tests/ && uv run ruff format --check muc_one_up/ tests/`
Run: `uv run mypy muc_one_up/`
Expected: Clean (no new errors).

- [ ] **Step 3: Verify assembly resolution centralization**

```bash
# Count assembly derivation locations in pipeline files
grep -rn "reference_assembly" muc_one_up/read_simulator/pipeline.py muc_one_up/read_simulator/ont_pipeline.py
# Expected: only the single AssemblyContext.from_configs() call per file
```

- [ ] **Step 4: Verify no stale temp paths**

```bash
# DiploidSimulationResult fields should be Optional
grep -A5 "class DiploidSimulationResult" muc_one_up/read_simulator/utils/diploid_handler.py
# Expected: str | None for hap*_fastq and hap*_reference fields
```

- [ ] **Step 5: Final test run**

Run: `uv run pytest --tb=short -q`
Expected: All pass. Phase 4B is complete.
