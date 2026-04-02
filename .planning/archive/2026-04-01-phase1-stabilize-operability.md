# Phase 1: Stabilize Basic Operability — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `import muc_one_up.cli`, `muconeup --help`, `mypy`, and `pytest` collection work in all environments, and enforce a single config-loading boundary.

**Architecture:** Lazy-load the `rfc8785` dependency so it's only imported at runtime when computing fingerprints. Guard the import in `provenance.py`. Replace all raw `json.load()` config reads in CLI commands with calls to `load_config()`. Add smoke tests to prevent regression.

**Tech Stack:** Python 3.10+, Click, pytest, mypy, ruff

---

### Task 1: Lazy-load `rfc8785` in `config_fingerprint.py`

**Files:**
- Modify: `muc_one_up/config_fingerprint.py:37`
- Test: `tests/test_config_fingerprint.py` (existing, verify still passes)

- [ ] **Step 1: Write test for lazy import behavior**

Create a new test file that verifies the module can be imported without `rfc8785` at the module level:

```python
# tests/test_config_fingerprint_import.py
"""Test that config_fingerprint module imports without rfc8785 at module level."""


def test_config_fingerprint_module_imports_without_rfc8785_at_top_level():
    """Importing config_fingerprint should not fail even if rfc8785 has issues.

    The actual rfc8785 import should only happen inside compute_config_fingerprint().
    We verify this by checking that the module-level code does not reference rfc8785.
    """
    import ast
    import inspect

    import muc_one_up.config_fingerprint as mod

    source = inspect.getsource(mod)
    tree = ast.parse(source)

    # Check that no top-level import statement imports rfc8785
    for node in ast.iter_child_nodes(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                assert alias.name != "rfc8785", (
                    "rfc8785 should not be imported at module level"
                )
        elif isinstance(node, ast.ImportFrom):
            assert node.module != "rfc8785", (
                "rfc8785 should not be imported at module level"
            )
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_config_fingerprint_import.py -v`
Expected: FAIL — `rfc8785` is currently imported at module level (line 37).

- [ ] **Step 3: Move `import rfc8785` into `compute_config_fingerprint()`**

Two changes in `muc_one_up/config_fingerprint.py`:

**Change 1:** Remove the module-level import at line 37:

```python
# REMOVE this line:
import rfc8785
```

**Change 2:** Add an `import rfc8785` with `ImportError` guard at the top of `compute_config_fingerprint()` (line 201). The existing function body stays the same — just wrap it with the import:

```python
    try:
        import rfc8785
    except ImportError:
        logger.warning(
            "rfc8785 package not available; config fingerprinting disabled"
        )
        return "error:fingerprint_failed:ImportError"

    try:
        # Step 1: Canonicalize config (filter sensitive fields)
        canonical_config = canonicalize_config(config)
        # ... rest of existing function body unchanged ...
```

The existing `except (TypeError, ValueError, RecursionError)` and `except Exception` handlers remain as-is.

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_config_fingerprint_import.py tests/test_config_fingerprint.py -v`
Expected: All PASS.

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/config_fingerprint.py tests/test_config_fingerprint_import.py
git commit -m "fix: lazy-load rfc8785 in config_fingerprint to prevent import-time failure"
```

---

### Task 2: Guard `config_fingerprint` import in `provenance.py`

**Files:**
- Modify: `muc_one_up/provenance.py:45`
- Test: `tests/test_provenance.py` (existing)

- [ ] **Step 1: Write test for guarded import**

```python
# tests/test_provenance_import.py
"""Test that provenance module handles missing config_fingerprint gracefully."""


def test_provenance_module_imports_successfully():
    """Provenance module should always import without errors."""
    import muc_one_up.provenance  # noqa: F401


def test_collect_provenance_metadata_works_without_rfc8785(monkeypatch):
    """Provenance collection should degrade gracefully if fingerprinting fails."""
    import time

    from muc_one_up.provenance import collect_provenance_metadata

    config = {"repeats": {"X": "ACGACT"}, "seed": 42}
    start = time.time()
    end = start + 1.0

    # Even if fingerprinting fails, provenance should return a dict
    result = collect_provenance_metadata(config, start, end)
    assert isinstance(result, dict)
    assert "software_version" in result
```

- [ ] **Step 2: Run test to verify it passes (baseline)**

Run: `pytest tests/test_provenance_import.py -v`
Expected: PASS (provenance already works when rfc8785 is installed). This establishes the baseline.

- [ ] **Step 3: Change provenance.py to use a guarded import**

In `muc_one_up/provenance.py`, replace the direct import with a guarded one:

```python
# REPLACE this line (line 45):
# from .config_fingerprint import compute_config_fingerprint

# WITH:
try:
    from .config_fingerprint import compute_config_fingerprint
except ImportError:
    def compute_config_fingerprint(config: dict) -> str:  # type: ignore[misc]
        """Fallback when config_fingerprint module cannot be loaded."""
        return "error:fingerprint_failed:ImportError"
```

- [ ] **Step 4: Run all provenance tests**

Run: `pytest tests/test_provenance_import.py tests/test_provenance.py -v`
Expected: All PASS.

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/provenance.py tests/test_provenance_import.py
git commit -m "fix: guard config_fingerprint import in provenance for graceful degradation"
```

---

### Task 3: Enforce single config-loading boundary in CLI commands

**Files:**
- Modify: `muc_one_up/cli/click_main.py:504-506`, `muc_one_up/cli/click_main.py:654-656`, `muc_one_up/cli/click_main.py:876-878`, `muc_one_up/cli/click_main.py:1046-1048`
- Test: `tests/test_config_boundary.py` (new)

- [ ] **Step 1: Write test that raw json.load is not used in CLI commands**

```python
# tests/test_config_boundary.py
"""Test that all CLI config reads go through the normalized load_config() loader."""

import ast
import inspect

import muc_one_up.cli.click_main as cli_module


def test_no_raw_json_load_in_cli_commands():
    """CLI commands must not use json.load() to read config files.

    All config loading must go through config.load_config() which normalizes
    flat constants into the nested assembly-aware shape.
    """
    source = inspect.getsource(cli_module)
    tree = ast.parse(source)

    violations = []
    for node in ast.walk(tree):
        # Look for calls like json.load(f)
        if isinstance(node, ast.Call):
            func = node.func
            if (
                isinstance(func, ast.Attribute)
                and func.attr == "load"
                and isinstance(func.value, ast.Name)
                and func.value.id == "json"
            ):
                violations.append(node.lineno)

    assert violations == [], (
        f"Found json.load() calls at lines {violations} in click_main.py. "
        f"Use config.load_config() instead to ensure constants normalization."
    )
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_config_boundary.py -v`
Expected: FAIL — finds `json.load()` at lines ~505, 655, 877, 1047.

- [ ] **Step 3: Replace raw json.load with load_config in Illumina reads command**

In `muc_one_up/cli/click_main.py`, find the Illumina reads command (~line 504-506):

```python
# REPLACE:
        config_path = Path(ctx.obj["config_path"])
        with config_path.open() as f:
            config = json.load(f)

# WITH:
        from ..config import load_config

        config = load_config(str(ctx.obj["config_path"]))
```

- [ ] **Step 4: Replace raw json.load with load_config in ONT reads command**

In `muc_one_up/cli/click_main.py`, find the ONT reads command (~line 654-656):

```python
# REPLACE:
        config_path = Path(ctx.obj["config_path"])
        with config_path.open() as f:
            config = json.load(f)

# WITH:
        from ..config import load_config

        config = load_config(str(ctx.obj["config_path"]))
```

- [ ] **Step 5: Replace raw json.load with load_config in PacBio reads command**

In `muc_one_up/cli/click_main.py`, find the PacBio reads command (~line 876-878):

```python
# REPLACE:
        config_path = Path(ctx.obj["config_path"])
        with config_path.open() as f:
            config = json.load(f)

# WITH:
        from ..config import load_config

        config = load_config(str(ctx.obj["config_path"]))
```

- [ ] **Step 6: Replace raw json.load with load_config in ORF analysis command**

In `muc_one_up/cli/click_main.py`, find the ORF analysis command (~line 1046-1048):

```python
# REPLACE:
        config_path = Path(ctx.obj["config_path"])
        with config_path.open() as f:
            config = json.load(f)

# WITH:
        from ..config import load_config

        config = load_config(str(ctx.obj["config_path"]))
```

- [ ] **Step 7: Fix ORF command's flat constant access**

After loading through `load_config()`, constants are nested by assembly. Fix the flat access at ~line 1049:

```python
# REPLACE:
        left_const = config.get("constants", {}).get("left")
        right_const = config.get("constants", {}).get("right")

# WITH:
        ref_assembly = config.get("reference_assembly", "hg38")
        assembly_constants = config.get("constants", {}).get(ref_assembly, {})
        left_const = assembly_constants.get("left")
        right_const = assembly_constants.get("right")
```

- [ ] **Step 8: Remove unused json import if no longer needed**

Check if `json` is still used elsewhere in `click_main.py`. If the only uses were for config loading, remove `import json` from the top. If `json` is still used for other purposes (e.g., `json.dump` for output), keep it.

- [ ] **Step 9: Run the boundary test and existing tests**

Run: `pytest tests/test_config_boundary.py tests/cli/ tests/test_click_cli.py tests/test_config.py -v`
Expected: All PASS.

- [ ] **Step 10: Commit**

```bash
git add muc_one_up/cli/click_main.py tests/test_config_boundary.py
git commit -m "fix: enforce single config-loading boundary in all CLI commands

All CLI config reads now go through load_config() which normalizes
flat constants into the nested assembly-aware shape. Fixes flat-key
access bug in ORF analysis command."
```

---

### Task 4: Fix flat constant access in `cli/analysis.py`

**Files:**
- Modify: `muc_one_up/cli/analysis.py:60-61`
- Test: `tests/test_config_boundary.py` (extend)

- [ ] **Step 1: Write test for analysis.py constant access**

Add to `tests/test_config_boundary.py`:

```python
import muc_one_up.cli.analysis as analysis_module


def test_no_flat_constant_access_in_analysis():
    """analysis.py must not use config['constants']['left'] (flat format).

    After load_config() normalization, constants are nested:
    config['constants']['hg38']['left'], not config['constants']['left'].
    """
    source = inspect.getsource(analysis_module)
    # Check for the flat access pattern: .get("constants", {}).get("left")
    assert '.get("constants", {}).get("left")' not in source, (
        "analysis.py uses flat constant access. Use assembly-keyed access instead."
    )
    assert '.get("constants", {}).get("right")' not in source, (
        "analysis.py uses flat constant access. Use assembly-keyed access instead."
    )
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_config_boundary.py::test_no_flat_constant_access_in_analysis -v`
Expected: FAIL — `analysis.py:60` uses `.get("constants", {}).get("left")`.

- [ ] **Step 3: Fix the flat constant access in analysis.py**

In `muc_one_up/cli/analysis.py`, find lines 60-61 inside `run_orf_prediction()`:

```python
# REPLACE:
            left_const_val = config.get("constants", {}).get("left")
            right_const_val = config.get("constants", {}).get("right")

# WITH:
            ref_assembly = config.get("reference_assembly", "hg38")
            assembly_constants = config.get("constants", {}).get(ref_assembly, {})
            left_const_val = assembly_constants.get("left")
            right_const_val = assembly_constants.get("right")
```

- [ ] **Step 4: Check for other flat constant access in analysis.py**

Search for any other `.get("constants", {}).get(` patterns in the file and fix them the same way.

- [ ] **Step 5: Run tests**

Run: `pytest tests/test_config_boundary.py -v`
Expected: All PASS.

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/cli/analysis.py tests/test_config_boundary.py
git commit -m "fix: use assembly-keyed constant access in analysis.py"
```

---

### Task 5: Add CLI smoke tests

**Files:**
- Create: `tests/test_cli_smoke.py`

- [ ] **Step 1: Write smoke tests**

```python
# tests/test_cli_smoke.py
"""Smoke tests for CLI import and basic operability.

These tests verify that the CLI can be imported and basic commands
work without a full environment (no config file, no external tools).
"""

from click.testing import CliRunner


def test_cli_module_imports():
    """The CLI module should import without errors."""
    from muc_one_up.cli import main  # noqa: F401


def test_cli_help_exits_zero():
    """muconeup --help should exit 0 without requiring config."""
    from muc_one_up.cli.click_main import cli

    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "MucOneUp" in result.output


def test_cli_version_exits_zero():
    """muconeup --version should exit 0 without requiring config."""
    from muc_one_up.cli.click_main import cli
    from muc_one_up.version import __version__

    runner = CliRunner()
    result = runner.invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.output


def test_simulate_help_exits_zero():
    """muconeup simulate --help should exit 0 without requiring config."""
    from muc_one_up.cli.click_main import cli

    runner = CliRunner()
    result = runner.invoke(cli, ["simulate", "--help"])
    assert result.exit_code == 0


def test_reads_help_exits_zero():
    """muconeup reads --help should exit 0 without requiring config."""
    from muc_one_up.cli.click_main import cli

    runner = CliRunner()
    result = runner.invoke(cli, ["reads", "--help"])
    assert result.exit_code == 0


def test_analyze_help_exits_zero():
    """muconeup analyze --help should exit 0 without requiring config."""
    from muc_one_up.cli.click_main import cli

    runner = CliRunner()
    result = runner.invoke(cli, ["analyze", "--help"])
    assert result.exit_code == 0
```

- [ ] **Step 2: Run smoke tests**

Run: `pytest tests/test_cli_smoke.py -v`
Expected: All PASS.

- [ ] **Step 3: Commit**

```bash
git add tests/test_cli_smoke.py
git commit -m "test: add CLI smoke tests for import and help accessibility"
```

---

### Task 6: Verify full test suite and CI readiness

**Files:** None (verification only)

- [ ] **Step 1: Run the full test suite**

Run: `pytest --tb=short -q`
Expected: All tests pass. Note any failures.

- [ ] **Step 2: Run ruff linter**

Run: `ruff check muc_one_up/ tests/`
Expected: No errors.

- [ ] **Step 3: Run ruff formatter**

Run: `ruff format --check muc_one_up/ tests/`
Expected: No formatting issues. If there are, run `ruff format muc_one_up/ tests/` to fix.

- [ ] **Step 4: Run mypy**

Run: `mypy muc_one_up/`
Expected: Passes (or only pre-existing errors unrelated to our changes).

- [ ] **Step 5: Verify the grep-based success criteria**

Run: `grep -rn "json.load" muc_one_up/cli/`
Expected: Zero hits outside of `config.py` (which legitimately uses `json.load` inside `load_config()`).

Run: `grep -rn 'get("constants", {}).get("left")' muc_one_up/`
Expected: Zero hits (all flat constant access has been replaced).

- [ ] **Step 6: Commit any formatting fixes if needed**

```bash
git add -u
git commit -m "style: fix formatting after Phase 1 changes"
```

- [ ] **Step 7: Run full test suite one final time**

Run: `pytest --tb=short -q`
Expected: All pass. Phase 1 is complete.
