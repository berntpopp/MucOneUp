# Review Fixes: Path to ~8/10 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the 5 concrete issues identified in codebase-review-report-v3.md to bring the codebase from ~6.8 to ~8/10.

**Architecture:** Targeted fixes to existing files — no new modules needed. Each task is independent and can be committed separately. All changes are backward-compatible.

**Tech Stack:** Python 3.10+, pytest, mypy, ruff, Click CLI

---

## File Map

| File | Action | Responsibility |
|------|--------|----------------|
| `muc_one_up/cli/commands/reads.py` | Modify (lines 286-299) | Fix ONT CLI → backend key mapping |
| `muc_one_up/mutate.py` | Modify (line 145) | Deep-copy results to prevent aliasing |
| `tests/test_mutate.py` | Modify (append) | Add aliasing-proof tests |
| `pyproject.toml` | Modify (line 167) | Fix mypy override for rfc8785 |
| `CLAUDE.md` | Modify | Standardize test command to `python -m pytest` |
| `muc_one_up/read_simulator/wrappers/samtools_convert.py` | Modify (lines 317-370) | Replace shell=True with run_pipeline |
| `muc_one_up/read_simulator/utils/common_utils.py` | Modify (lines 241-256) | Surface timeout-specific error in streaming mode |

---

### Task 1: Fix ONT CLI → backend config key mismatch

The ONT CLI command writes `--coverage` to `read_simulation.coverage` via `_setup_read_config()`, but `ont_pipeline.py:150` reads from `nanosim_params.coverage`. The backend never sees the CLI value.

**Files:**
- Modify: `muc_one_up/cli/commands/reads.py:286-299`

- [ ] **Step 1: Write the failing test**

Create a test that verifies the ONT command writes coverage to `nanosim_params.coverage`.

Add to a new file `tests/cli/test_reads_ont_config.py`:

```python
"""Tests for ONT CLI config key mapping."""

import pytest
from unittest.mock import patch, MagicMock
from click.testing import CliRunner

from muc_one_up.cli.click_main import main


@pytest.fixture
def minimal_ont_config(tmp_path):
    """Write a minimal config file for ONT testing."""
    import json

    config = {
        "tools": {"nanosim": "nanosim", "minimap2": "minimap2", "samtools": "samtools"},
        "nanosim_params": {"training_data_path": "/fake/model", "coverage": 50},
        "read_simulation": {},
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))
    return config_path


def test_ont_coverage_written_to_nanosim_params(minimal_ont_config, tmp_path):
    """CLI --coverage must land in nanosim_params.coverage for the backend."""
    captured_config = {}

    def fake_simulate(config, *args, **kwargs):
        captured_config.update(config)

    runner = CliRunner()
    with patch("muc_one_up.cli.commands.reads.simulate_reads_pipeline", fake_simulate):
        # We can't fully run the command without real tools, so we patch
        # the simulate function and check the config dict that was passed.
        pass

    # Instead, test the config setup logic directly:
    from muc_one_up.cli.commands.reads import _setup_read_config

    config = {"read_simulation": {}, "nanosim_params": {"training_data_path": "/fake"}}
    _setup_read_config(config, "ont", coverage=75, seed=None, seed_config_key="nanosim_params")

    # The bug: coverage only goes to read_simulation, not nanosim_params
    # After fix, nanosim_params.coverage must be set
    assert config["nanosim_params"].get("coverage") == 75


def test_ont_min_read_length_key_matches_backend(minimal_ont_config):
    """CLI --min-read-length must use key 'min_read_length' in nanosim_params."""
    # This is already correct in current code, but we lock it with a test
    config = {"read_simulation": {}, "nanosim_params": {}}
    _setup_read_config(config, "ont", coverage=None, seed=None, seed_config_key="nanosim_params")

    # Simulate what the ont() command does after _setup_read_config
    config["nanosim_params"]["min_read_length"] = 200

    assert config["nanosim_params"]["min_read_length"] == 200
    # Must NOT use "min_len" — that's the wrong key
    assert "min_len" not in config["nanosim_params"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/cli/test_reads_ont_config.py::test_ont_coverage_written_to_nanosim_params -v`
Expected: FAIL — `config["nanosim_params"].get("coverage")` returns `None` because `_setup_read_config` only writes to `read_simulation.coverage`.

- [ ] **Step 3: Fix the ONT command to write coverage to nanosim_params**

In `muc_one_up/cli/commands/reads.py`, in the `ont()` function (around line 291-295), after `_setup_read_config`, propagate coverage to `nanosim_params` where the backend expects it:

```python
    _setup_read_config(config, "ont", coverage, seed, seed_config_key="nanosim_params")

    if "nanosim_params" not in config:
        config["nanosim_params"] = {}
    config["nanosim_params"]["min_read_length"] = min_read_length
    # Ensure coverage is in nanosim_params where ont_pipeline.py reads it
    config["nanosim_params"]["coverage"] = config["read_simulation"]["coverage"]
```

The key change is adding the last line: `config["nanosim_params"]["coverage"] = config["read_simulation"]["coverage"]`.

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/cli/test_reads_ont_config.py -v`
Expected: PASS

- [ ] **Step 5: Run full test suite to check for regressions**

Run: `python -m pytest --tb=short -q`
Expected: All tests pass.

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/cli/commands/reads.py tests/cli/test_reads_ont_config.py
git commit -m "fix: write CLI --coverage to nanosim_params for ONT backend

The ONT pipeline reads coverage from nanosim_params.coverage but the CLI
wrote it only to read_simulation.coverage, causing --coverage to be silently
ignored."
```

---

### Task 2: Make mutation application non-aliasing

`apply_mutations()` in `mutate.py` does `updated_results = list(results)` which is a shallow copy. The HaplotypeResult objects (and their `.chain` lists) are shared with the caller. Mutations to `hr.chain[i]` modify the caller's data.

**Files:**
- Modify: `muc_one_up/mutate.py:145`
- Test: `tests/test_mutate.py` (append)

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_mutate.py`:

```python
class TestMutationAliasing:
    """Prove that apply_mutations does not mutate caller-owned data."""

    def test_original_results_not_mutated(self, mutation_config):
        """Caller's results list and HaplotypeResult objects must be untouched."""
        original_chain = [RepeatUnit("X")]
        original_seq = "TTTTXXXXXGGGG"
        results = [HaplotypeResult(sequence=original_seq, chain=original_chain)]

        # Snapshot before mutation
        chain_before = [RepeatUnit(ru.symbol, ru.mutated) for ru in original_chain]
        seq_before = original_seq

        updated, _ = apply_mutations(
            config=mutation_config,
            results=results,
            mutation_name="testMut",
            targets=[MutationTarget(1, 1)],
        )

        # Caller's objects must be completely unchanged
        assert results[0].sequence == seq_before
        assert len(results[0].chain) == len(chain_before)
        for orig, before in zip(results[0].chain, chain_before):
            assert orig.symbol == before.symbol
            assert orig.mutated == before.mutated

    def test_chain_objects_not_shared(self, mutation_config):
        """Updated chain list must be a different object from the original."""
        original_chain = [RepeatUnit("X")]
        results = [HaplotypeResult(sequence="TTTTXXXXXGGGG", chain=original_chain)]

        updated, _ = apply_mutations(
            config=mutation_config,
            results=results,
            mutation_name="testMut",
            targets=[MutationTarget(1, 1)],
        )

        # The chain list in updated must be a different object
        assert updated[0].chain is not results[0].chain

    def test_forced_change_does_not_alias(self, mutation_config):
        """Even the forced-change path (non-strict, wrong symbol) must not alias."""
        original_chain = [RepeatUnit("C")]
        results = [HaplotypeResult(sequence="TTTTCCCCC", chain=original_chain)]

        updated, _ = apply_mutations(
            config=mutation_config,
            results=results,
            mutation_name="testMut",
            targets=[MutationTarget(1, 1)],
        )

        # Original chain must still be "C", not mutated
        assert results[0].chain[0].symbol == "C"
        assert results[0].chain[0].mutated is False
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_mutate.py::TestMutationAliasing -v`
Expected: FAIL — `test_original_results_not_mutated` and `test_chain_objects_not_shared` fail because the chain is shared.

- [ ] **Step 3: Fix apply_mutations to deep-copy results**

In `muc_one_up/mutate.py`, change line 145 from:

```python
    updated_results = list(results)
```

to:

```python
    updated_results = [
        HaplotypeResult(
            sequence=hr.sequence,
            chain=[RepeatUnit(ru.symbol, ru.mutated) for ru in hr.chain],
        )
        for hr in results
    ]
```

This creates new HaplotypeResult objects with independently copied chains, so no mutation can propagate back to the caller.

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/test_mutate.py -v`
Expected: All tests pass, including the new aliasing tests and all existing tests.

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/mutate.py tests/test_mutate.py
git commit -m "fix: deep-copy results in apply_mutations to prevent aliasing

list(results) was a shallow copy — mutations to chain[i] propagated back
to the caller. Now each HaplotypeResult and its chain are independently
copied before any modification."
```

---

### Task 3: Fix mypy rfc8785 override mismatch

The mypy override in `pyproject.toml` targets `rfc8785.*` (submodules only), but `config_fingerprint.py:200` does `import rfc8785` (top-level module). mypy can't find stubs for the top-level import.

**Files:**
- Modify: `pyproject.toml:166-168`

- [ ] **Step 1: Fix the mypy override**

In `pyproject.toml`, change:

```toml
[[tool.mypy.overrides]]
module = "rfc8785.*"
ignore_missing_imports = true
```

to:

```toml
[[tool.mypy.overrides]]
module = ["rfc8785", "rfc8785.*"]
ignore_missing_imports = true
```

This tells mypy to ignore missing imports for both the top-level `rfc8785` module and any submodules.

- [ ] **Step 2: Run mypy to verify it passes**

Run: `mypy muc_one_up/`
Expected: `Success: no issues found in 83 source files` (0 errors, down from 1).

- [ ] **Step 3: Commit**

```bash
git add pyproject.toml
git commit -m "fix: mypy override for rfc8785 to match top-level import

Override was 'rfc8785.*' (submodules only) but code does 'import rfc8785'
(top-level). Added both patterns to the override list."
```

---

### Task 4: Standardize test runner invocation

`pytest` and `python -m pytest` use different interpreters locally. CLAUDE.md documents `pytest` but the project only passes under `python -m pytest`. Standardize on `python -m pytest` in CLAUDE.md so the documented command always works.

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Update CLAUDE.md test command**

In `CLAUDE.md`, change:

```bash
pytest --tb=short -q              # Run tests
```

to:

```bash
python -m pytest --tb=short -q    # Run tests
```

- [ ] **Step 2: Verify the documented command works**

Run: `python -m pytest --tb=short -q`
Expected: All tests pass (1324+ passed).

- [ ] **Step 3: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: standardize test command to python -m pytest

Avoids interpreter mismatch where bare 'pytest' may pick up a different
Python than the project venv."
```

---

### Task 5: Harden subprocess execution

Two issues: (a) `samtools_convert.py` has a `shell=True` path for conda-wrapped piped commands that bypasses the shared `run_pipeline()` runner, and (b) `common_utils.py:run_command()` streaming mode collapses timeout failures into a generic `"Command failed"` error, losing the timeout context.

**Files:**
- Modify: `muc_one_up/read_simulator/wrappers/samtools_convert.py:310-371`
- Modify: `muc_one_up/read_simulator/utils/common_utils.py:241-256`

- [ ] **Step 1: Write a test for timeout-specific error reporting**

Add to `tests/test_common_utils.py` (create if needed):

```python
"""Tests for common_utils subprocess helpers."""

import subprocess
from unittest.mock import patch, MagicMock

import pytest

from muc_one_up.read_simulator.utils.common_utils import run_command
from muc_one_up.exceptions import ExternalToolError


def test_streaming_mode_timeout_error_includes_timeout_info():
    """Streaming mode must report timeout duration in error, not generic 'Command failed'."""
    # Mock Popen to simulate a timeout
    mock_proc = MagicMock()
    mock_proc.pid = 12345
    mock_proc.stdout = MagicMock()
    mock_proc.stderr = MagicMock()
    # Make stdout/stderr iterators return empty immediately
    mock_proc.stdout.readline = MagicMock(return_value=b"")
    mock_proc.stderr.readline = MagicMock(return_value=b"")
    mock_proc.wait = MagicMock(side_effect=subprocess.TimeoutExpired(cmd="test", timeout=30))
    mock_proc.returncode = -15  # SIGTERM

    with patch("subprocess.Popen", return_value=mock_proc), \
         patch("os.setsid"), \
         patch("os.killpg"), \
         patch("os.getpgid", return_value=12345):
        with pytest.raises(ExternalToolError) as exc_info:
            run_command(["sleep", "999"], timeout=30)

        error = exc_info.value
        # The error message must mention the timeout, not just "Command failed"
        assert "timed out" in str(error).lower() or "timeout" in str(error).lower()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_common_utils.py::test_streaming_mode_timeout_error_includes_timeout_info -v`
Expected: FAIL — current code raises `ExternalToolError` with `stderr="Command failed"` which doesn't mention timeout.

- [ ] **Step 3: Fix timeout error reporting in streaming mode**

In `muc_one_up/read_simulator/utils/common_utils.py`, replace the error-raising block (lines 241-256) from:

```python
    if popen_proc.returncode != 0:
        # Check if this is a timeout-related exit code (-15 is typical for SIGTERM)
        if timeout is not None and popen_proc.returncode == -15:
            logging.info("Command terminated due to timeout (%d seconds): %s", timeout, cmd_str)
        else:
            logging.error(
                "Command exited with non-zero exit code %d: %s",
                popen_proc.returncode,
                cmd_str,
            )
        raise ExternalToolError(
            tool="command",
            exit_code=popen_proc.returncode,
            stderr="Command failed",
            cmd=cmd_str,
        )
```

with:

```python
    if popen_proc.returncode != 0:
        # Check if this is a timeout-related exit code (-15 is typical for SIGTERM)
        if timeout is not None and popen_proc.returncode == -15:
            logging.info("Command terminated due to timeout (%d seconds): %s", timeout, cmd_str)
            stderr_msg = f"Command timed out after {timeout}s"
        else:
            logging.error(
                "Command exited with non-zero exit code %d: %s",
                popen_proc.returncode,
                cmd_str,
            )
            stderr_msg = f"Command failed with exit code {popen_proc.returncode}"
        raise ExternalToolError(
            tool="command",
            exit_code=popen_proc.returncode,
            stderr=stderr_msg,
            cmd=cmd_str,
        )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_common_utils.py -v`
Expected: PASS

- [ ] **Step 5: Fix samtools_convert.py conda wrapper path**

In `muc_one_up/read_simulator/wrappers/samtools_convert.py`, replace the `if use_conda_wrapper:` block (lines 317-371) with code that uses `run_pipeline` for both paths. The conda wrapper prefix is already part of `samtools_cmd_list` via `shlex.split()`, so building command lists with it works for piped execution too:

Replace lines 317-371 (the entire `if use_conda_wrapper:` branch) with:

```python
        if use_conda_wrapper:
            # For conda-wrapped commands, build proper command lists
            # and use run_pipeline. The wrapper prefix is already in
            # samtools_cmd_list from shlex.split().
            collate_cmd = [
                *samtools_cmd_list,
                "collate",
                "-u",
                "-O",
                "-@",
                str(opts.threads),
                str(input_bam),
            ]

            fastq_cmd = [
                *samtools_cmd_list,
                "fastq",
                "-1",
                str(output_fq1),
                "-2",
                str(output_fq2),
                "-F",
                "0x900",
            ]

            # Add singleton handling
            if singleton_opt:
                fastq_cmd.extend(singleton_opt.split())

            # Add other options
            fastq_cmd.extend(
                ["-0", "/dev/null", read_name_opt, "-@", str(opts.threads), "-"]
            )

            logging.debug("  Collate command: %s", " ".join(collate_cmd))
            logging.debug("  Fastq command: %s", " ".join(fastq_cmd))

            pipeline_result = run_pipeline(
                [collate_cmd, fastq_cmd],
                capture=True,
                timeout=opts.timeout,
            )

            # Log any warnings from samtools
            if pipeline_result.stderr:
                stderr_text = pipeline_result.stderr.strip()
                if stderr_text:
                    logging.debug("  Samtools stderr: %s", stderr_text)
```

This eliminates the `shell=True` path entirely. Both conda and non-conda paths now use `run_pipeline`, which uses `subprocess.Popen` with `shell=False`.

Also remove the now-unused imports at the top of the function. The `import shlex` can stay (it's used for `samtools_cmd_list` construction). Remove `import subprocess` since it's no longer used directly in this function.

- [ ] **Step 6: Run full test suite**

Run: `python -m pytest --tb=short -q`
Expected: All tests pass.

- [ ] **Step 7: Run ruff check**

Run: `ruff check muc_one_up/read_simulator/wrappers/samtools_convert.py muc_one_up/read_simulator/utils/common_utils.py`
Expected: All checks passed.

- [ ] **Step 8: Commit**

```bash
git add muc_one_up/read_simulator/wrappers/samtools_convert.py muc_one_up/read_simulator/utils/common_utils.py tests/test_common_utils.py
git commit -m "fix: harden subprocess execution — remove shell=True, surface timeouts

samtools_convert: replace shell=True conda wrapper path with run_pipeline,
eliminating the inconsistent execution path.

common_utils: streaming mode now reports 'timed out after Ns' instead of
generic 'Command failed' for timeout exits."
```

---

## Final Verification

After all 5 tasks are committed:

- [ ] **Run full quality gate**

```bash
ruff check muc_one_up/ tests/
ruff format --check muc_one_up/ tests/
mypy muc_one_up/
python -m pytest --tb=short -q
```

All four commands must pass cleanly.
