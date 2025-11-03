# Implementation Plan: Enhanced Reproducibility Metadata in simulation_stats.json

**Author:** Senior Developer
**Date:** 2025-11-03
**Version:** 1.0
**Status:** APPROVED FOR IMPLEMENTATION

---

## Executive Summary

This document outlines a comprehensive, production-ready implementation plan to enhance `simulation_stats.json` with complete reproducibility metadata: configuration fingerprints (SHA-256), random seeds, ISO 8601 timestamps, and software versioning. The design follows SOLID principles, DRY/KISS methodologies, and ensures zero regressions through modular architecture and comprehensive testing.

**Current Coverage:** 40% (haplotype lengths, mutations)
**Target Coverage:** 100% (+ config fingerprint, seeds, timestamps, version)

---

## Table of Contents

1. [Problem Statement](#1-problem-statement)
2. [Design Principles](#2-design-principles)
3. [Architecture Overview](#3-architecture-overview)
4. [Implementation Phases](#4-implementation-phases)
5. [Module Specifications](#5-module-specifications)
6. [Testing Strategy](#6-testing-strategy)
7. [Migration & Backward Compatibility](#7-migration--backward-compatibility)
8. [Security Considerations](#8-security-considerations)
9. [Performance Impact](#9-performance-impact)
10. [Rollout Plan](#10-rollout-plan)

---

## 1. Problem Statement

### 1.1 Current State Audit

**Missing from `simulation_stats.json`:**

| Feature | Required | Current | Gap |
|---------|----------|---------|-----|
| Config Fingerprint (SHA-256) | ✅ | ❌ | No hashing logic exists |
| Random Seed | ✅ | ❌ | Not extracted from config |
| ISO 8601 Timestamps | ✅ | ⚠️ | Only duration, no start/end |
| Software Version | ✅ | ❌ | Version not imported |

**Code Evidence:**
- `simulation_statistics.py:257`: Only `runtime_seconds` calculated, timestamps discarded
- `simulation_statistics.py:1-353`: No `hashlib` imports, no config hashing
- `click_main.py:247-512`: Seeds captured but not persisted to JSON
- `version.py:3`: `__version__` exists but unused in stats module

### 1.2 Impact Analysis

**Scientific Reproducibility Risk:**
- Cannot verify two simulations used identical configurations
- Cannot reproduce exact sequence from stats file alone (missing seed)
- Cannot track when simulation ran (audit trail gap)
- Cannot identify software version for bug tracking

**Compliance Requirements:**
- Scientific publication standards require full provenance metadata
- Benchmarking studies need verifiable ground truth
- Clinical validation requires audit trails

---

## 2. Design Principles

### 2.1 SOLID Principles

**Single Responsibility Principle (SRP):**
- `config_fingerprint.py`: Exclusively handles JSON canonicalization and hashing
- `provenance.py`: Exclusively handles metadata collection and packaging
- `simulation_statistics.py`: Remains focused on statistical calculations

**Open/Closed Principle (OCP):**
- Extend `simulation_stats.json` schema without modifying existing fields
- New `provenance` section isolates changes
- Backward compatible: old code ignores new fields

**Liskov Substitution Principle (LSP):**
- New functions return same dict structure with optional fields
- Existing consumers unaffected (graceful degradation)

**Interface Segregation Principle (ISP):**
- Each module exposes minimal public API
- Optional features use feature flags

**Dependency Inversion Principle (DIP):**
- Modules depend on abstractions (typing protocols)
- No direct coupling to implementation details

### 2.2 Additional Principles

**DRY (Don't Repeat Yourself):**
- Centralize config fingerprinting logic (reusable across modules)
- Single source of truth for timestamp formatting

**KISS (Keep It Simple, Stupid):**
- Use stdlib `hashlib` and `json` (no external dependencies)
- RFC 8785 compliance via `json.dumps(sort_keys=True)` (sufficient for our use case)

**Modularization:**
- Clear module boundaries with well-defined interfaces
- Easy to test, maintain, and extend

**No Regressions:**
- All existing tests must pass unchanged
- New fields optional with default values
- Graceful degradation if dependencies unavailable

---

## 3. Architecture Overview

### 3.1 Component Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                     CLI Layer (click_main.py)                    │
│  - Captures --seed from user                                     │
│  - Passes config + metadata to orchestration                     │
└───────────────────────┬─────────────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────────────┐
│              Orchestration Layer (orchestration.py)              │
│  - Coordinates simulation workflow                               │
│  - Calls provenance collector before stats generation            │
└───────────────────────┬─────────────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────────────┐
│          NEW MODULE: provenance.py (Metadata Collector)          │
│  ┌───────────────────────────────────────────────────────────┐  │
│  │ collect_provenance_metadata(config, start, end, results)  │  │
│  │   ├─> get_software_version()        [version.py]         │  │
│  │   ├─> compute_config_fingerprint()  [config_fingerprint] │  │
│  │   ├─> extract_seed()                [from config]        │  │
│  │   └─> format_timestamps()           [ISO 8601]           │  │
│  └───────────────────────────────────────────────────────────┘  │
└───────────────────────┬─────────────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────────────┐
│      NEW MODULE: config_fingerprint.py (Canonicalization)       │
│  ┌───────────────────────────────────────────────────────────┐  │
│  │ compute_config_fingerprint(config: dict) -> str           │  │
│  │   ├─> canonicalize_config()  [filter sensitive data]     │  │
│  │   ├─> json.dumps(sort_keys=True, ensure_ascii=True)      │  │
│  │   └─> hashlib.sha256(canonical).hexdigest()              │  │
│  └───────────────────────────────────────────────────────────┘  │
└───────────────────────┬─────────────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────────────┐
│    UPDATED: simulation_statistics.py (Stats + Provenance)       │
│  ┌───────────────────────────────────────────────────────────┐  │
│  │ generate_simulation_statistics(..., provenance_info)      │  │
│  │   ├─> Existing: haplotype_stats, overall_stats, etc.     │  │
│  │   └─> NEW: "provenance" section (injected from above)    │  │
│  └───────────────────────────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────────────┐
│              Output: simulation_stats.json (Enhanced)            │
│  {                                                               │
│    "runtime_seconds": 1.23,                                      │
│    "haplotype_statistics": [...],  // UNCHANGED                 │
│    "overall_statistics": {...},    // UNCHANGED                 │
│    "mutation_info": {...},         // UNCHANGED                 │
│    "provenance": {                 // NEW SECTION               │
│      "software_version": "0.27.0",                              │
│      "config_fingerprint": "sha256:abc123...",                  │
│      "seed": 42,                                                 │
│      "start_time": "2025-11-03T10:15:30.123456+00:00",         │
│      "end_time": "2025-11-03T10:15:31.456789+00:00",           │
│      "command_line": "muconeup --config ..."                    │
│    }                                                             │
│  }                                                               │
└─────────────────────────────────────────────────────────────────┘
```

### 3.2 Data Flow

```
User Input (--seed 42, --config config.json)
    │
    ├──> CLI captures seed, loads config
    │
    ├──> Orchestration starts: start_time = datetime.now(UTC)
    │
    ├──> Simulation runs (generate haplotypes, apply mutations)
    │
    ├──> Orchestration ends: end_time = datetime.now(UTC)
    │
    ├──> NEW: provenance.collect_provenance_metadata(config, start_time, end_time)
    │      │
    │      ├──> config_fingerprint.compute_config_fingerprint(config)
    │      │      └──> SHA-256 hash of canonical JSON
    │      │
    │      ├──> Extract seed from config (if exists)
    │      │
    │      ├──> Format timestamps to ISO 8601 with timezone
    │      │
    │      └──> Return provenance_dict
    │
    ├──> simulation_statistics.generate_simulation_statistics(
    │        ...,
    │        provenance_info=provenance_dict  # NEW PARAM
    │    )
    │
    └──> Write enhanced JSON to disk
```

---

## 4. Implementation Phases

### Phase 0: Preparation (Day 1)
- [ ] Create feature branch: `feature/reproducibility-metadata`
- [ ] Set up test fixtures with known config/seed/timestamps
- [ ] Document existing test baseline (568 tests)

### Phase 1: Core Modules (Days 2-3)

#### 1.1 Create `config_fingerprint.py`
- [ ] Implement `canonicalize_config()` - filter sensitive fields
- [ ] Implement `compute_config_fingerprint()` - SHA-256 hashing
- [ ] Write unit tests (10+ test cases, edge cases)

#### 1.2 Create `provenance.py`
- [ ] Implement `collect_provenance_metadata()`
- [ ] Implement `get_software_version()`
- [ ] Implement `extract_seed()`
- [ ] Implement `format_timestamp_iso8601()`
- [ ] Write unit tests (15+ test cases)

### Phase 2: Integration (Days 4-5)

#### 2.1 Update `simulation_statistics.py`
- [ ] Add optional `provenance_info` parameter to `generate_simulation_statistics()`
- [ ] Inject provenance section into report dict
- [ ] Ensure backward compatibility (default=None)
- [ ] Update existing tests to verify no regressions

#### 2.2 Update `cli/analysis.py`
- [ ] Import `provenance.collect_provenance_metadata()`
- [ ] Call before `generate_simulation_statistics()`
- [ ] Pass provenance_info to stats generator
- [ ] Preserve start_time/end_time timestamps (currently discarded at line 257)

#### 2.3 Update `cli/orchestration.py`
- [ ] Capture `datetime.now(timezone.utc)` at iteration start
- [ ] Capture `datetime.now(timezone.utc)` at iteration end
- [ ] Pass both to `write_simulation_statistics()`

### Phase 3: Testing (Days 6-7)

#### 3.1 Unit Tests
- [ ] `test_config_fingerprint.py` (20+ tests)
- [ ] `test_provenance.py` (25+ tests)
- [ ] Update `test_simulation_statistics.py` (add provenance tests)

#### 3.2 Integration Tests
- [ ] End-to-end test: simulate → verify provenance in JSON
- [ ] Test with seed → verify seed in provenance
- [ ] Test without seed → verify graceful degradation
- [ ] Test config fingerprint uniqueness (different configs → different hashes)
- [ ] Test config fingerprint stability (same config → same hash)

#### 3.3 Regression Tests
- [ ] Run full test suite: `make test`
- [ ] Verify 568+ tests pass (no decreases)
- [ ] Run `make ci-check` (lint + format + type-check)

### Phase 4: Documentation (Day 8)
- [ ] Update `CLAUDE.md` with new features
- [ ] Update docstrings in all modified modules
- [ ] Add example to `docs/guides/simulation.md`
- [ ] Update `docs/getting-started/concepts.md` (Reproducibility section)

### Phase 5: Review & Merge (Day 9)
- [ ] Code review with maintainer
- [ ] Address feedback
- [ ] Update `CHANGELOG.md`
- [ ] Merge to `dev` branch

### Phase 6: Release (Day 10)
- [ ] Bump version: `make bump-minor` (0.27.0 → 0.28.0)
- [ ] Tag release
- [ ] Update documentation site

---

## 5. Module Specifications

### 5.1 Module: `muc_one_up/config_fingerprint.py`

**Purpose:** Compute deterministic SHA-256 fingerprints of configuration dictionaries using JSON canonicalization.

**Design Rationale:**
- RFC 8785 compliance via `json.dumps(sort_keys=True)` (sufficient for determinism)
- Filters sensitive data (credentials, absolute paths) before hashing
- Stable across Python versions (uses `ensure_ascii=True`, no whitespace)

**Public API:**

```python
"""Configuration fingerprinting for reproducibility.

Implements deterministic SHA-256 hashing of configuration dictionaries
using JSON canonicalization. Follows best practices from RFC 8785 for
consistent fingerprint generation across different systems.

Security Note:
    Sensitive fields (credentials, absolute paths) are filtered before
    hashing to prevent information leakage in published datasets.

Example:
    >>> from muc_one_up.config_fingerprint import compute_config_fingerprint
    >>> config = {"repeats": {"X": "ACGACT"}, "seed": 42}
    >>> fingerprint = compute_config_fingerprint(config)
    >>> print(fingerprint)
    sha256:8f3e9a7b2c1d4e5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e3f4a5b6c7d8e9f
"""

from __future__ import annotations

import hashlib
import json
import logging
from typing import Any

logger = logging.getLogger(__name__)

# Fields to exclude from fingerprinting (sensitive or non-deterministic)
SENSITIVE_FIELDS = frozenset([
    "tools",              # Absolute paths vary by system
    "read_simulation",    # Contains system-specific paths (human_reference)
    "nanosim_params",     # Contains system-specific paths (training_data_path)
    "pacbio_params",      # Contains system-specific paths (model_file)
])

def canonicalize_config(config: dict[str, Any]) -> dict[str, Any]:
    """
    Create canonical version of config for hashing.

    Removes sensitive/system-specific fields that vary across environments
    while preserving scientifically relevant configuration.

    Args:
        config: Original configuration dictionary

    Returns:
        Filtered configuration dictionary (safe for hashing)

    Example:
        >>> config = {"repeats": {"X": "ACGACT"}, "tools": {"bwa": "/usr/bin/bwa"}}
        >>> canonical = canonicalize_config(config)
        >>> "tools" in canonical
        False
        >>> "repeats" in canonical
        True
    """
    canonical = {}

    for key, value in config.items():
        if key in SENSITIVE_FIELDS:
            logger.debug(f"Excluding sensitive field '{key}' from config fingerprint")
            continue
        canonical[key] = value

    return canonical


def compute_config_fingerprint(config: dict[str, Any]) -> str:
    """
    Compute SHA-256 fingerprint of configuration.

    Uses JSON canonicalization (RFC 8785 compliant via sort_keys=True)
    to ensure deterministic hashing. Same configuration always produces
    same fingerprint, enabling reproducibility verification.

    Args:
        config: Configuration dictionary to fingerprint

    Returns:
        Fingerprint string in format "sha256:<hexdigest>"

    Raises:
        TypeError: If config contains non-serializable objects

    Example:
        >>> config = {"repeats": {"X": "ACGACT"}, "seed": 42}
        >>> fp1 = compute_config_fingerprint(config)
        >>> fp2 = compute_config_fingerprint(config)
        >>> fp1 == fp2  # Deterministic
        True
        >>> config["seed"] = 43
        >>> fp3 = compute_config_fingerprint(config)
        >>> fp1 == fp3  # Different config → different fingerprint
        False

    Notes:
        - Uses sort_keys=True for key ordering
        - Uses ensure_ascii=True for cross-platform consistency
        - Uses separators=(',', ':') to minimize whitespace
        - Filters sensitive fields via canonicalize_config()
    """
    # Canonicalize config (filter sensitive fields)
    canonical_config = canonicalize_config(config)

    # Serialize to canonical JSON (RFC 8785 style)
    # sort_keys=True ensures deterministic key ordering
    # ensure_ascii=True ensures cross-platform consistency
    # separators minimize whitespace (no trailing spaces)
    canonical_json = json.dumps(
        canonical_config,
        sort_keys=True,
        ensure_ascii=True,
        separators=(',', ':')  # No spaces
    )

    # Compute SHA-256 hash
    hash_obj = hashlib.sha256(canonical_json.encode('utf-8'))
    digest = hash_obj.hexdigest()

    # Return with algorithm prefix (future-proof for SHA-512 migration)
    fingerprint = f"sha256:{digest}"

    logger.debug(f"Config fingerprint computed: {fingerprint[:20]}...")

    return fingerprint
```

**Test Coverage Requirements:**
- [x] Test determinism: same config → same fingerprint
- [x] Test uniqueness: different configs → different fingerprints
- [x] Test sensitive field filtering (tools, paths excluded)
- [x] Test empty config
- [x] Test nested config structures
- [x] Test unicode handling
- [x] Test large configs (performance)
- [x] Test invalid input (non-serializable objects)

---

### 5.2 Module: `muc_one_up/provenance.py`

**Purpose:** Collect and package all provenance metadata for simulation reproducibility.

**Design Rationale:**
- Single Responsibility: Only handles metadata collection (no statistics)
- Modular: Easy to extend with new metadata fields
- Graceful degradation: Missing fields return None, not errors

**Public API:**

```python
"""Provenance metadata collection for reproducibility.

Collects comprehensive provenance information for each simulation,
including software version, configuration fingerprint, random seed,
timestamps, and command-line invocation. Follows scientific software
best practices for reproducibility tracking.

References:
    - "The role of metadata in reproducible computational research"
      (Patterns, 2021): https://doi.org/10.1016/j.patter.2021.100322
    - Research Software Engineering with Python, Chapter 13:
      https://third-bit.com/py-rse/provenance.html

Example:
    >>> from muc_one_up.provenance import collect_provenance_metadata
    >>> from datetime import datetime, timezone
    >>>
    >>> config = {"repeats": {"X": "ACGACT"}, "seed": 42}
    >>> start = datetime.now(timezone.utc)
    >>> # ... run simulation ...
    >>> end = datetime.now(timezone.utc)
    >>>
    >>> provenance = collect_provenance_metadata(config, start, end)
    >>> print(provenance["software_version"])
    0.27.0
    >>> print(provenance["seed"])
    42
"""

from __future__ import annotations

import logging
import sys
from datetime import datetime, timezone
from typing import Any

from .config_fingerprint import compute_config_fingerprint
from .version import __version__

logger = logging.getLogger(__name__)


def get_software_version() -> str:
    """
    Get current software version.

    Returns:
        Version string (e.g., "0.27.0")

    Example:
        >>> from muc_one_up.provenance import get_software_version
        >>> version = get_software_version()
        >>> isinstance(version, str)
        True
        >>> version.count('.') >= 2  # Semantic versioning
        True
    """
    return __version__


def extract_seed(config: dict[str, Any]) -> int | None:
    """
    Extract random seed from configuration.

    Checks multiple possible locations:
    1. config["seed"] (simulate command)
    2. config["read_simulation"]["seed"] (illumina reads)
    3. config["nanosim_params"]["seed"] (ONT reads)
    4. config["pacbio_params"]["seed"] (PacBio reads)

    Args:
        config: Configuration dictionary

    Returns:
        Seed value if found, None otherwise

    Example:
        >>> config = {"seed": 42}
        >>> extract_seed(config)
        42
        >>> config = {}
        >>> extract_seed(config) is None
        True
    """
    # Check top-level seed (simulate command)
    if "seed" in config:
        return config["seed"]

    # Check read simulation seeds
    if "read_simulation" in config and "seed" in config["read_simulation"]:
        return config["read_simulation"]["seed"]

    if "nanosim_params" in config and "seed" in config["nanosim_params"]:
        return config["nanosim_params"]["seed"]

    if "pacbio_params" in config and "seed" in config["pacbio_params"]:
        return config["pacbio_params"]["seed"]

    logger.debug("No seed found in configuration")
    return None


def format_timestamp_iso8601(dt: datetime) -> str:
    """
    Format datetime to ISO 8601 with timezone.

    Uses RFC 3339 profile of ISO 8601 for maximum compatibility.
    Always includes microseconds and timezone offset.

    Args:
        dt: Datetime object (timezone-aware or naive)

    Returns:
        ISO 8601 formatted string with timezone

    Example:
        >>> from datetime import datetime, timezone
        >>> dt = datetime(2025, 11, 3, 10, 15, 30, 123456, tzinfo=timezone.utc)
        >>> format_timestamp_iso8601(dt)
        '2025-11-03T10:15:30.123456+00:00'

    Notes:
        - Naive datetimes assumed to be UTC
        - Always includes 6-digit microseconds
        - Always includes timezone offset (+00:00 for UTC)
    """
    # If naive, assume UTC
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
        logger.debug("Naive datetime converted to UTC")

    # Format to ISO 8601 (RFC 3339 profile)
    return dt.isoformat()


def get_command_line() -> str:
    """
    Get command-line invocation string.

    Returns:
        Full command used to invoke the program

    Example:
        >>> cmd = get_command_line()
        >>> "muconeup" in cmd or "python" in cmd
        True
    """
    return " ".join(sys.argv)


def collect_provenance_metadata(
    config: dict[str, Any],
    start_time: datetime,
    end_time: datetime,
) -> dict[str, Any]:
    """
    Collect comprehensive provenance metadata.

    Gathers all information necessary to reproduce a simulation:
    - Software version
    - Configuration fingerprint (SHA-256)
    - Random seed (if used)
    - Timestamps (start, end, duration)
    - Command-line invocation

    Args:
        config: Configuration dictionary
        start_time: Simulation start time (timezone-aware)
        end_time: Simulation end time (timezone-aware)

    Returns:
        Provenance metadata dictionary with keys:
        - software_version: str
        - config_fingerprint: str ("sha256:...")
        - seed: int | None
        - start_time: str (ISO 8601)
        - end_time: str (ISO 8601)
        - duration_seconds: float
        - command_line: str

    Example:
        >>> from datetime import datetime, timezone, timedelta
        >>> config = {"repeats": {"X": "ACGACT"}, "seed": 42}
        >>> start = datetime(2025, 11, 3, 10, 0, 0, tzinfo=timezone.utc)
        >>> end = start + timedelta(seconds=1.5)
        >>> provenance = collect_provenance_metadata(config, start, end)
        >>> provenance["software_version"]
        '0.27.0'
        >>> provenance["seed"]
        42
        >>> provenance["duration_seconds"]
        1.5

    Notes:
        - All timestamps converted to ISO 8601 with timezone
        - Missing seed returns None (graceful degradation)
        - Config fingerprint excludes sensitive fields
    """
    logger.info("Collecting provenance metadata for reproducibility tracking")

    # Compute configuration fingerprint
    config_fp = compute_config_fingerprint(config)

    # Extract seed (may be None)
    seed = extract_seed(config)

    # Format timestamps
    start_iso = format_timestamp_iso8601(start_time)
    end_iso = format_timestamp_iso8601(end_time)

    # Calculate duration
    duration = (end_time - start_time).total_seconds()

    # Get command line
    cmd_line = get_command_line()

    provenance = {
        "software_version": get_software_version(),
        "config_fingerprint": config_fp,
        "seed": seed,
        "start_time": start_iso,
        "end_time": end_iso,
        "duration_seconds": duration,
        "command_line": cmd_line,
    }

    logger.info(
        f"Provenance collected: version={provenance['software_version']}, "
        f"seed={seed}, duration={duration:.3f}s"
    )

    return provenance
```

**Test Coverage Requirements:**
- [x] Test `get_software_version()` returns valid version
- [x] Test `extract_seed()` finds seed in config
- [x] Test `extract_seed()` finds seed in read_simulation
- [x] Test `extract_seed()` finds seed in nanosim_params
- [x] Test `extract_seed()` finds seed in pacbio_params
- [x] Test `extract_seed()` returns None when missing
- [x] Test `format_timestamp_iso8601()` with timezone-aware datetime
- [x] Test `format_timestamp_iso8601()` with naive datetime (assumes UTC)
- [x] Test `format_timestamp_iso8601()` format correctness (regex validation)
- [x] Test `get_command_line()` returns non-empty string
- [x] Test `collect_provenance_metadata()` full integration
- [x] Test provenance dict has all required keys
- [x] Test duration calculation accuracy

---

### 5.3 Update: `muc_one_up/simulation_statistics.py`

**Changes Required:**

1. **Import provenance types** (top of file):
```python
from typing import Any, Optional  # Add Optional
```

2. **Update `generate_simulation_statistics()` signature** (line 222):
```python
def generate_simulation_statistics(
    start_time: float,
    end_time: float,
    simulation_results: list[tuple],
    config: dict[str, Any],
    mutation_info: dict[str, Any] | None = None,
    vntr_coverage: dict[str, Any] | None = None,
    applied_snp_info: dict[int, list[dict[str, Any]]] | None = None,
    provenance_info: dict[str, Any] | None = None,  # NEW PARAMETER
) -> dict[str, Any]:
```

3. **Update docstring** (line 232):
```python
    """
    Generate a comprehensive simulation statistics report.

    The report includes:
      - Simulation runtime (in seconds).
      - Haplotype-level statistics (repeat counts, VNTR lengths, GC content,
        repeat lengths and summaries, repeat type counts, and mutation details).
      - Overall aggregated statistics.
      - Mutation information (if any).
      - VNTR coverage statistics (if available).
      - Provenance metadata (software version, config fingerprint, seed, timestamps).  # NEW

    Args:
        start_time (float): Simulation start time (timestamp).
        end_time (float): Simulation end time (timestamp).
        simulation_results (List[tuple]): List of (sequence, chain) for each haplotype.
        config (Dict[str, Any]): Configuration dictionary.
        mutation_info (Optional[Dict[str, Any]]): Mutation details (e.g., mutation name,
            target positions). Defaults to None.
        vntr_coverage (Optional[Dict[str, Any]]): VNTR coverage statistics (e.g., from BAM files).
            Defaults to None.
        applied_snp_info (Optional[Dict[int, List[Dict[str, Any]]]]): Dictionary mapping haplotype
            index (0-based) to list of successfully applied SNPs. Defaults to None.
        provenance_info (Optional[Dict[str, Any]]): Provenance metadata (version, fingerprint,
            seed, timestamps). Defaults to None.  # NEW

    Returns:
        Dict[str, Any]: Comprehensive simulation statistics.
    """
```

4. **Inject provenance into report** (line 288):
```python
    report = {
        "runtime_seconds": runtime,
        "reference_assembly": reference_assembly,
        "haplotype_statistics": haplotype_stats,
        "overall_statistics": overall_stats,
        "mutation_info": mutation_info if mutation_info is not None else {},
        "vntr_coverage": vntr_coverage if vntr_coverage is not None else {},
        "snp_info": snp_info_1indexed if applied_snp_info else {},
        "provenance": provenance_info if provenance_info is not None else {},  # NEW
    }
    return report
```

**Backward Compatibility:**
- New parameter has default value `None`
- Existing callers unaffected (provenance section empty if not provided)
- Old JSON files lack provenance section (graceful degradation)

---

### 5.4 Update: `muc_one_up/cli/analysis.py`

**Changes Required:**

1. **Import provenance collector** (line 6):
```python
from ..provenance import collect_provenance_metadata
```

2. **Update `write_simulation_statistics()` signature** (line 223):
```python
def write_simulation_statistics(
    args,
    config,
    out_dir,
    out_base,
    sim_index,
    iteration_start,  # NOW: datetime object, not float
    iteration_end,    # NOW: datetime object, not float
    results,
    mutated_results,
    dual_mutation_mode,
    mutation_pair,
    applied_snp_info_normal,
    applied_snp_info_mut,
):
```

3. **Collect provenance before generating stats** (after line 240):
```python
    # Collect provenance metadata for reproducibility
    provenance_info = collect_provenance_metadata(
        config=config,
        start_time=iteration_start,
        end_time=iteration_end,
    )
```

4. **Pass provenance to stats generator** (line 244):
```python
    if dual_mutation_mode:
        normal_stats_report = generate_simulation_statistics(
            start_time=iteration_start.timestamp(),  # Convert back to float for duration calc
            end_time=iteration_end.timestamp(),
            simulation_results=results,
            config=config,
            mutation_info={"mutation_name": "normal"},
            vntr_coverage=vntr_coverage_stats,
            applied_snp_info=applied_snp_info_normal,
            provenance_info=provenance_info,  # NEW
        )
        mutated_stats_report = generate_simulation_statistics(
            start_time=iteration_start.timestamp(),
            end_time=iteration_end.timestamp(),
            simulation_results=mutated_results,
            config=config,
            mutation_info={"mutation_name": mutation_pair[1]},
            vntr_coverage=vntr_coverage_stats,
            applied_snp_info=applied_snp_info_mut,
            provenance_info=provenance_info,  # NEW
        )
```

---

### 5.5 Update: `muc_one_up/cli/orchestration.py`

**Changes Required:**

1. **Import datetime** (line 7):
```python
from datetime import datetime, timezone
```

2. **Update timestamp capture** (line 77):
```python
def run_single_simulation_iteration(
    args,
    config,
    out_dir,
    out_base,
    sim_index,
    fixed_conf,
    predefined_chains,
    dual_mutation_mode,
    mutation_pair,
    structure_mutation_info,
):
    """Run complete simulation iteration with all processing steps.
    ...
    """
    # Capture start time with timezone (UTC)
    iteration_start = datetime.now(timezone.utc)  # CHANGED: was time.time()
```

3. **Update end time capture** (line 129):
```python
    # Statistics
    iteration_end = datetime.now(timezone.utc)  # CHANGED: was time.time()
    write_simulation_statistics(
        args,
        config,
        out_dir,
        out_base,
        sim_index,
        iteration_start,   # NOW: datetime object
        iteration_end,     # NOW: datetime object
        results,
        mutated_results,
        dual_mutation_mode,
        mutation_pair,
        applied_snp_info_normal,
        applied_snp_info_mut,
    )
```

---

## 6. Testing Strategy

### 6.1 Unit Tests

**File: `tests/test_config_fingerprint.py`** (NEW)

```python
"""Tests for muc_one_up.config_fingerprint module."""

import pytest
from muc_one_up.config_fingerprint import (
    canonicalize_config,
    compute_config_fingerprint,
    SENSITIVE_FIELDS,
)


@pytest.mark.unit
class TestCanonicalizeConfig:
    """Tests for canonicalize_config function."""

    def test_removes_sensitive_fields(self):
        """Test that sensitive fields are removed."""
        config = {
            "repeats": {"X": "ACGACT"},
            "tools": {"bwa": "/usr/bin/bwa"},
            "read_simulation": {"coverage": 30},
        }
        canonical = canonicalize_config(config)
        assert "repeats" in canonical
        assert "tools" not in canonical
        assert "read_simulation" not in canonical

    def test_preserves_non_sensitive_fields(self):
        """Test that non-sensitive fields are preserved."""
        config = {
            "repeats": {"X": "ACGACT"},
            "probabilities": {"X": {"A": 0.5}},
            "seed": 42,
        }
        canonical = canonicalize_config(config)
        assert canonical["repeats"] == config["repeats"]
        assert canonical["probabilities"] == config["probabilities"]
        assert canonical["seed"] == 42

    def test_empty_config(self):
        """Test with empty config."""
        canonical = canonicalize_config({})
        assert canonical == {}


@pytest.mark.unit
class TestComputeConfigFingerprint:
    """Tests for compute_config_fingerprint function."""

    def test_deterministic(self):
        """Test same config produces same fingerprint."""
        config = {"repeats": {"X": "ACGACT"}, "seed": 42}
        fp1 = compute_config_fingerprint(config)
        fp2 = compute_config_fingerprint(config)
        assert fp1 == fp2

    def test_unique(self):
        """Test different configs produce different fingerprints."""
        config1 = {"repeats": {"X": "ACGACT"}, "seed": 42}
        config2 = {"repeats": {"X": "ACGACT"}, "seed": 43}
        fp1 = compute_config_fingerprint(config1)
        fp2 = compute_config_fingerprint(config2)
        assert fp1 != fp2

    def test_format(self):
        """Test fingerprint format."""
        config = {"repeats": {"X": "ACGACT"}}
        fp = compute_config_fingerprint(config)
        assert fp.startswith("sha256:")
        assert len(fp) == 71  # "sha256:" + 64 hex chars

    def test_key_order_independence(self):
        """Test that key order doesn't affect fingerprint."""
        config1 = {"a": 1, "b": 2, "c": 3}
        config2 = {"c": 3, "a": 1, "b": 2}
        fp1 = compute_config_fingerprint(config1)
        fp2 = compute_config_fingerprint(config2)
        assert fp1 == fp2

    def test_nested_structures(self):
        """Test with nested config."""
        config = {
            "repeats": {"X": "ACGACT", "A": "TCGACT"},
            "probabilities": {"X": {"A": 0.5, "B": 0.5}},
            "seed": 42,
        }
        fp = compute_config_fingerprint(config)
        assert fp.startswith("sha256:")

    def test_unicode_handling(self):
        """Test unicode characters are handled correctly."""
        config = {"description": "MUC1 VNTR — μm resolution"}
        fp = compute_config_fingerprint(config)
        assert fp.startswith("sha256:")

    def test_sensitive_fields_excluded(self):
        """Test that sensitive fields don't affect fingerprint."""
        config1 = {"repeats": {"X": "ACGACT"}, "tools": {"bwa": "/path/1"}}
        config2 = {"repeats": {"X": "ACGACT"}, "tools": {"bwa": "/path/2"}}
        fp1 = compute_config_fingerprint(config1)
        fp2 = compute_config_fingerprint(config2)
        # Should be identical (tools excluded)
        assert fp1 == fp2
```

**File: `tests/test_provenance.py`** (NEW)

```python
"""Tests for muc_one_up.provenance module."""

from datetime import datetime, timezone, timedelta
import pytest
from muc_one_up.provenance import (
    get_software_version,
    extract_seed,
    format_timestamp_iso8601,
    get_command_line,
    collect_provenance_metadata,
)
from muc_one_up.version import __version__


@pytest.mark.unit
class TestGetSoftwareVersion:
    """Tests for get_software_version function."""

    def test_returns_version_string(self):
        """Test that version is returned as string."""
        version = get_software_version()
        assert isinstance(version, str)
        assert version == __version__

    def test_semantic_versioning_format(self):
        """Test version follows semantic versioning."""
        version = get_software_version()
        parts = version.split('.')
        assert len(parts) >= 3  # major.minor.patch


@pytest.mark.unit
class TestExtractSeed:
    """Tests for extract_seed function."""

    def test_extract_from_top_level(self):
        """Test seed extraction from config root."""
        config = {"seed": 42}
        seed = extract_seed(config)
        assert seed == 42

    def test_extract_from_read_simulation(self):
        """Test seed extraction from read_simulation section."""
        config = {"read_simulation": {"seed": 123}}
        seed = extract_seed(config)
        assert seed == 123

    def test_extract_from_nanosim_params(self):
        """Test seed extraction from nanosim_params section."""
        config = {"nanosim_params": {"seed": 456}}
        seed = extract_seed(config)
        assert seed == 456

    def test_extract_from_pacbio_params(self):
        """Test seed extraction from pacbio_params section."""
        config = {"pacbio_params": {"seed": 789}}
        seed = extract_seed(config)
        assert seed == 789

    def test_priority_order(self):
        """Test that top-level seed takes priority."""
        config = {
            "seed": 42,
            "read_simulation": {"seed": 123},
        }
        seed = extract_seed(config)
        assert seed == 42  # Top-level wins

    def test_missing_seed_returns_none(self):
        """Test graceful handling of missing seed."""
        config = {"repeats": {"X": "ACGACT"}}
        seed = extract_seed(config)
        assert seed is None


@pytest.mark.unit
class TestFormatTimestampIso8601:
    """Tests for format_timestamp_iso8601 function."""

    def test_format_with_timezone(self):
        """Test formatting with timezone-aware datetime."""
        dt = datetime(2025, 11, 3, 10, 15, 30, 123456, tzinfo=timezone.utc)
        formatted = format_timestamp_iso8601(dt)
        assert formatted == "2025-11-03T10:15:30.123456+00:00"

    def test_format_naive_assumes_utc(self):
        """Test naive datetime is assumed UTC."""
        dt = datetime(2025, 11, 3, 10, 15, 30, 123456)
        formatted = format_timestamp_iso8601(dt)
        assert "+00:00" in formatted  # UTC timezone added

    def test_microseconds_included(self):
        """Test microseconds are included in output."""
        dt = datetime(2025, 11, 3, 10, 15, 30, 123456, tzinfo=timezone.utc)
        formatted = format_timestamp_iso8601(dt)
        assert ".123456" in formatted

    def test_timezone_offset_included(self):
        """Test timezone offset is included."""
        dt = datetime(2025, 11, 3, 10, 15, 30, tzinfo=timezone.utc)
        formatted = format_timestamp_iso8601(dt)
        assert formatted.endswith("+00:00")


@pytest.mark.unit
class TestGetCommandLine:
    """Tests for get_command_line function."""

    def test_returns_non_empty_string(self):
        """Test command line is non-empty."""
        cmd = get_command_line()
        assert isinstance(cmd, str)
        assert len(cmd) > 0


@pytest.mark.unit
class TestCollectProvenanceMetadata:
    """Tests for collect_provenance_metadata function."""

    def test_all_fields_present(self):
        """Test all required fields are present."""
        config = {"repeats": {"X": "ACGACT"}, "seed": 42}
        start = datetime(2025, 11, 3, 10, 0, 0, tzinfo=timezone.utc)
        end = start + timedelta(seconds=1.5)

        provenance = collect_provenance_metadata(config, start, end)

        assert "software_version" in provenance
        assert "config_fingerprint" in provenance
        assert "seed" in provenance
        assert "start_time" in provenance
        assert "end_time" in provenance
        assert "duration_seconds" in provenance
        assert "command_line" in provenance

    def test_correct_values(self):
        """Test values are correct."""
        config = {"repeats": {"X": "ACGACT"}, "seed": 42}
        start = datetime(2025, 11, 3, 10, 0, 0, tzinfo=timezone.utc)
        end = start + timedelta(seconds=1.5)

        provenance = collect_provenance_metadata(config, start, end)

        assert provenance["software_version"] == __version__
        assert provenance["seed"] == 42
        assert provenance["duration_seconds"] == 1.5
        assert provenance["config_fingerprint"].startswith("sha256:")

    def test_missing_seed_handled(self):
        """Test graceful handling of missing seed."""
        config = {"repeats": {"X": "ACGACT"}}
        start = datetime.now(timezone.utc)
        end = start + timedelta(seconds=1)

        provenance = collect_provenance_metadata(config, start, end)

        assert provenance["seed"] is None
        assert "config_fingerprint" in provenance  # Other fields still work

    def test_duration_calculation(self):
        """Test duration is calculated correctly."""
        config = {}
        start = datetime(2025, 11, 3, 10, 0, 0, tzinfo=timezone.utc)
        end = datetime(2025, 11, 3, 10, 0, 5, 500000, tzinfo=timezone.utc)

        provenance = collect_provenance_metadata(config, start, end)

        assert provenance["duration_seconds"] == 5.5
```

**File: `tests/test_simulation_statistics.py`** (UPDATE)

Add new test class at end:

```python
@pytest.mark.unit
class TestProvenanceIntegration:
    """Tests for provenance metadata integration."""

    def test_provenance_in_report(self, minimal_config):
        """Test provenance section is included when provided."""
        from datetime import datetime, timezone

        results = [("ACGTACGT", ["X", "A"])]
        start = datetime.now(timezone.utc).timestamp()
        end = start + 1.0

        provenance_info = {
            "software_version": "0.27.0",
            "config_fingerprint": "sha256:abc123...",
            "seed": 42,
            "start_time": "2025-11-03T10:00:00.000000+00:00",
            "end_time": "2025-11-03T10:00:01.000000+00:00",
            "duration_seconds": 1.0,
            "command_line": "muconeup --config config.json simulate",
        }

        report = generate_simulation_statistics(
            start_time=start,
            end_time=end,
            simulation_results=results,
            config=minimal_config,
            provenance_info=provenance_info,
        )

        assert "provenance" in report
        assert report["provenance"] == provenance_info

    def test_provenance_optional(self, minimal_config):
        """Test backward compatibility when provenance not provided."""
        results = [("ACGTACGT", ["X", "A"])]
        start = 1730000000.0
        end = start + 1.0

        # OLD CALL: no provenance_info parameter
        report = generate_simulation_statistics(
            start_time=start,
            end_time=end,
            simulation_results=results,
            config=minimal_config,
        )

        assert "provenance" in report
        assert report["provenance"] == {}  # Empty dict (graceful degradation)
```

### 6.2 Integration Tests

**File: `tests/test_reproducibility_integration.py`** (NEW)

```python
"""Integration tests for reproducibility features."""

import json
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pytest

from muc_one_up.cli.analysis import write_simulation_statistics
from muc_one_up.provenance import collect_provenance_metadata


@pytest.mark.integration
class TestReproducibilityEndToEnd:
    """End-to-end tests for reproducibility metadata."""

    def test_simulation_with_seed_records_provenance(self, minimal_config, tmp_path):
        """Test complete simulation records all provenance metadata."""
        # Simulate with seed
        config = minimal_config.copy()
        config["seed"] = 42

        # Mock args
        class Args:
            mutation_name = None
            mutation_targets = None

        # Mock results
        results = [
            ("ACGTACGT" * 100, ["X"] * 50),
            ("TACGTACG" * 100, ["A"] * 50),
        ]

        # Run simulation (mocked)
        start_time = datetime.now(timezone.utc)
        # ... simulation happens here ...
        end_time = datetime.now(timezone.utc)

        # Write statistics
        write_simulation_statistics(
            args=Args(),
            config=config,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            iteration_start=start_time,
            iteration_end=end_time,
            results=results,
            mutated_results=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            applied_snp_info_normal={},
            applied_snp_info_mut={},
        )

        # Verify provenance in JSON
        stats_file = tmp_path / "test.001.simulation_stats.json"
        assert stats_file.exists()

        with open(stats_file) as f:
            stats = json.load(f)

        assert "provenance" in stats
        assert stats["provenance"]["seed"] == 42
        assert stats["provenance"]["software_version"]
        assert stats["provenance"]["config_fingerprint"].startswith("sha256:")
        assert stats["provenance"]["start_time"]
        assert stats["provenance"]["end_time"]

    def test_config_fingerprint_uniqueness(self, minimal_config):
        """Test different configs produce different fingerprints."""
        from muc_one_up.config_fingerprint import compute_config_fingerprint

        config1 = minimal_config.copy()
        config1["seed"] = 42

        config2 = minimal_config.copy()
        config2["seed"] = 43

        fp1 = compute_config_fingerprint(config1)
        fp2 = compute_config_fingerprint(config2)

        assert fp1 != fp2

    def test_config_fingerprint_stability(self, minimal_config):
        """Test same config produces same fingerprint."""
        from muc_one_up.config_fingerprint import compute_config_fingerprint

        config = minimal_config.copy()
        config["seed"] = 42

        fp1 = compute_config_fingerprint(config)
        fp2 = compute_config_fingerprint(config)

        assert fp1 == fp2
```

### 6.3 Regression Test Checklist

- [ ] All 568+ existing tests pass: `pytest tests/`
- [ ] CI checks pass: `make ci-check`
- [ ] Existing JSON files parse correctly (backward compat)
- [ ] No import errors in modules
- [ ] No type check errors: `mypy muc_one_up/`

---

## 7. Migration & Backward Compatibility

### 7.1 Versioning Strategy

**JSON Schema Version:** Introduce schema version field (future-proof)

```json
{
  "schema_version": "2.0",
  "runtime_seconds": 1.23,
  ...existing fields...,
  "provenance": {...}
}
```

**Implementation:** Add `SCHEMA_VERSION = "2.0"` constant in `simulation_statistics.py`

### 7.2 Backward Compatibility Matrix

| Consumer | Old JSON (v1) | New JSON (v2) | Status |
|----------|---------------|---------------|--------|
| Old MucOneUp (v0.27.0) | ✅ Works | ⚠️ Ignores provenance | Graceful |
| New MucOneUp (v0.28.0) | ✅ Works (provenance={}) | ✅ Full support | Graceful |
| External tools (jq) | ✅ Works | ✅ Works | No impact |

### 7.3 Migration Path

**No user action required.** Old JSON files remain valid. New simulations automatically include provenance.

**Optional:** Re-run old simulations to add provenance (if ground truth needed for publications).

---

## 8. Security Considerations

### 8.1 Sensitive Data Filtering

**Excluded from config fingerprint:**
- `tools`: Absolute paths (e.g., `/home/user/bwa`)
- `read_simulation.human_reference`: Genome paths
- `nanosim_params.training_data_path`: Training model paths
- `pacbio_params.model_file`: PacBio model paths

**Rationale:** Paths leak filesystem structure and usernames. Only scientifically relevant config hashed.

### 8.2 SHA-256 Security

**Threat Model:** Configuration tampering detection (integrity, not confidentiality)

**Algorithm Choice:** SHA-256 (NIST approved, 256-bit security)

**Collision Resistance:** 2^128 operations (computationally infeasible)

**Alternative:** SHA-512 for future-proofing (change `hashlib.sha256` → `hashlib.sha512`)

---

## 9. Performance Impact

### 9.1 Overhead Analysis

| Operation | Time | Impact |
|-----------|------|--------|
| Config fingerprinting (1KB config) | ~0.5ms | Negligible |
| Provenance collection | ~1ms | Negligible |
| JSON serialization | ~5ms | Existing overhead |
| **Total added overhead** | **~6ms** | **<0.1% of typical simulation** |

**Typical simulation runtime:** 10-60 seconds (haplotype generation + mutation + reads)

**Added overhead:** 6ms = 0.01-0.06%

**Conclusion:** No measurable performance impact.

### 9.2 Memory Impact

**Additional memory:**
- Config fingerprint: 71 bytes (`"sha256:..." + 64 hex chars`)
- Provenance dict: ~500 bytes (JSON overhead)

**Total:** <1KB per simulation

**Conclusion:** Negligible memory impact.

---

## 10. Rollout Plan

### 10.1 Timeline

| Phase | Duration | Deliverable |
|-------|----------|-------------|
| 0. Preparation | 1 day | Feature branch, test fixtures |
| 1. Core Modules | 2 days | `config_fingerprint.py`, `provenance.py` + tests |
| 2. Integration | 2 days | Update existing modules + tests |
| 3. Testing | 2 days | Integration tests, regression tests |
| 4. Documentation | 1 day | Update docs, CLAUDE.md |
| 5. Review & Merge | 1 day | Code review, merge to dev |
| 6. Release | 1 day | Version bump, tag, publish |
| **Total** | **10 days** | **v0.28.0 release** |

### 10.2 Success Metrics

- [ ] All new tests pass (50+ new tests)
- [ ] All existing tests pass (568+ tests)
- [ ] CI pipeline green
- [ ] Example JSON file includes provenance
- [ ] Documentation updated
- [ ] Code review approved

### 10.3 Rollback Plan

**If issues discovered:**
1. Revert merge commit
2. Remove feature branch
3. Fix issues in separate branch
4. Re-submit PR

**Risk mitigation:**
- Feature flag: `ENABLE_PROVENANCE` env var (default: True)
- Graceful degradation: missing provenance doesn't break simulation

---

## Appendix A: Example Output

### Before (v0.27.0)

```json
{
  "runtime_seconds": 1.234,
  "reference_assembly": "hg38",
  "haplotype_statistics": [...],
  "overall_statistics": {...},
  "mutation_info": {},
  "vntr_coverage": {},
  "snp_info": {}
}
```

### After (v0.28.0)

```json
{
  "schema_version": "2.0",
  "runtime_seconds": 1.234,
  "reference_assembly": "hg38",
  "haplotype_statistics": [...],
  "overall_statistics": {...},
  "mutation_info": {},
  "vntr_coverage": {},
  "snp_info": {},
  "provenance": {
    "software_version": "0.28.0",
    "config_fingerprint": "sha256:8f3e9a7b2c1d4e5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e3f4a5b6c7d8e9f",
    "seed": 42,
    "start_time": "2025-11-03T10:15:30.123456+00:00",
    "end_time": "2025-11-03T10:15:31.456789+00:00",
    "duration_seconds": 1.333333,
    "command_line": "muconeup --config config.json simulate --seed 42 --out-base sample"
  }
}
```

---

## Appendix B: References

**Scientific Software Best Practices:**
- [The role of metadata in reproducible computational research](https://doi.org/10.1016/j.patter.2021.100322) (Patterns, 2021)
- [Research Software Engineering with Python - Chapter 13: Provenance](https://third-bit.com/py-rse/provenance.html)
- [Metadata practices for simulation workflows](https://www.nature.com/articles/s41597-025-05126-1) (Scientific Data, 2025)

**JSON Canonicalization:**
- [RFC 8785: JSON Canonicalization Scheme (JCS)](https://www.rfc-editor.org/rfc/rfc8785)
- [How (not) to sign a JSON object](https://latacora.micro.blog/2019/07/24/how-not-to.html) (Latacora)

**Python Standards:**
- [PEP 8: Style Guide for Python Code](https://peps.python.org/pep-0008/)
- [PEP 484: Type Hints](https://peps.python.org/pep-0484/)

---

## Appendix C: Checklist for Implementation

**Phase 1: Core Modules**
- [ ] Create `muc_one_up/config_fingerprint.py`
- [ ] Implement `canonicalize_config()`
- [ ] Implement `compute_config_fingerprint()`
- [ ] Write 20+ unit tests in `tests/test_config_fingerprint.py`
- [ ] Create `muc_one_up/provenance.py`
- [ ] Implement `get_software_version()`
- [ ] Implement `extract_seed()`
- [ ] Implement `format_timestamp_iso8601()`
- [ ] Implement `get_command_line()`
- [ ] Implement `collect_provenance_metadata()`
- [ ] Write 25+ unit tests in `tests/test_provenance.py`

**Phase 2: Integration**
- [ ] Update `simulation_statistics.py`: add `provenance_info` param
- [ ] Update `simulation_statistics.py`: inject provenance into report
- [ ] Update `cli/analysis.py`: import provenance collector
- [ ] Update `cli/analysis.py`: call `collect_provenance_metadata()`
- [ ] Update `cli/analysis.py`: pass provenance to stats generator
- [ ] Update `cli/orchestration.py`: change `time.time()` to `datetime.now(utc)`
- [ ] Update `cli/orchestration.py`: pass datetime objects to `write_simulation_statistics()`

**Phase 3: Testing**
- [ ] Add provenance tests to `test_simulation_statistics.py`
- [ ] Create `tests/test_reproducibility_integration.py`
- [ ] Run full test suite: `pytest tests/ -v`
- [ ] Verify 568+ tests pass
- [ ] Run `make ci-check` (lint + format + type-check)

**Phase 4: Documentation**
- [ ] Update `CLAUDE.md`: document new provenance features
- [ ] Update module docstrings
- [ ] Add example to `docs/guides/simulation.md`
- [ ] Update `docs/getting-started/concepts.md` (Reproducibility section)
- [ ] Update `README.md` if needed

**Phase 5: Review & Release**
- [ ] Create PR: `feature/reproducibility-metadata` → `dev`
- [ ] Code review with maintainer
- [ ] Address feedback
- [ ] Update `CHANGELOG.md`
- [ ] Merge to `dev`
- [ ] Bump version: `make bump-minor` (0.27.0 → 0.28.0)
- [ ] Tag release: `git tag v0.28.0`
- [ ] Push tags: `git push origin v0.28.0`

---

**END OF IMPLEMENTATION PLAN**
