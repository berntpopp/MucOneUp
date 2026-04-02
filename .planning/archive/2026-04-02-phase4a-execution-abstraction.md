# Phase 4A: Execution Abstraction & Subprocess Migration — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace all direct `subprocess.run`/`Popen` calls across the read-simulator codebase with a unified execution abstraction (`RunResult`, extended `run_command`, `run_pipeline`), so that subprocess usage is centralized in a single module.

**Architecture:** Extend `common_utils.py` with a `RunResult` dataclass and two new modes: `capture=True` for output capture, and `run_pipeline()` for pipe chains. Migrate the 33 direct subprocess calls across 5 files (bwa_wrapper, samtools_wrapper, nanosim_wrapper, utils/samtools, utils/bed, utils/tool_version) to use the new abstractions. The conda `shell=True` exception in samtools_wrapper stays as-is (it's properly justified and secured).

**Tech Stack:** Python 3.10+, subprocess, dataclasses, pytest

**Spec:** `.planning/specs/2026-04-01-codebase-refactoring-design.md` (sections 4.2, 4.3, 4.6)

**Note:** Phase 4B (AssemblyContext, temp-file fix, alignment dedup) is a separate plan.

---

## File Structure

### Modified Files
| File | Changes |
|------|---------|
| `muc_one_up/read_simulator/utils/common_utils.py` | Add `RunResult` dataclass; extend `run_command()` with `capture` mode; add `run_pipeline()` |
| `muc_one_up/read_simulator/wrappers/bwa_wrapper.py` | Replace 7 direct subprocess calls with `run_command`/`run_pipeline` |
| `muc_one_up/read_simulator/wrappers/samtools_wrapper.py` | Replace 10 direct subprocess calls (keep 1 conda shell=True exception) |
| `muc_one_up/read_simulator/wrappers/nanosim_wrapper.py` | Replace 1 direct subprocess call |
| `muc_one_up/read_simulator/utils/samtools.py` | Replace 9 direct subprocess calls |
| `muc_one_up/read_simulator/utils/bed.py` | Replace 2 direct subprocess calls |
| `muc_one_up/read_simulator/utils/tool_version.py` | Replace 1 direct subprocess call |

### New Test Files
| File | Responsibility |
|------|---------------|
| `tests/read_simulator/test_run_result.py` | Tests for RunResult, capture mode, run_pipeline |

---

## Task 1: Add RunResult Dataclass and Capture Mode to run_command

**Files:**
- Modify: `muc_one_up/read_simulator/utils/common_utils.py`
- Create: `tests/read_simulator/test_run_result.py`

- [ ] **Step 1: Write failing tests for RunResult and capture mode**

```python
# tests/read_simulator/test_run_result.py
"""Tests for RunResult and extended run_command."""

from __future__ import annotations

import subprocess
import sys

import pytest

from muc_one_up.read_simulator.utils.common_utils import RunResult, run_command


class TestRunResult:
    """Tests for RunResult dataclass."""

    def test_create(self):
        r = RunResult(returncode=0, stdout="hello", stderr="", command="echo hello")
        assert r.returncode == 0
        assert r.stdout == "hello"
        assert r.stderr == ""
        assert r.command == "echo hello"

    def test_none_stdout_in_streaming_mode(self):
        r = RunResult(returncode=0, stdout=None, stderr=None, command="cmd")
        assert r.stdout is None


class TestRunCommandCapture:
    """Tests for run_command with capture=True."""

    def test_capture_stdout(self):
        result = run_command(
            [sys.executable, "-c", "print('hello world')"],
            capture=True,
        )
        assert isinstance(result, RunResult)
        assert result.returncode == 0
        assert "hello world" in result.stdout

    def test_capture_stderr(self):
        result = run_command(
            [sys.executable, "-c", "import sys; sys.stderr.write('err msg\\n')"],
            capture=True,
        )
        assert result.returncode == 0
        assert "err msg" in result.stderr

    def test_capture_nonzero_raises(self):
        with pytest.raises(Exception):
            run_command(
                [sys.executable, "-c", "import sys; sys.exit(1)"],
                capture=True,
            )

    def test_streaming_mode_returns_run_result(self):
        """Default streaming mode now also returns RunResult."""
        result = run_command(
            [sys.executable, "-c", "print('streamed')"],
        )
        assert isinstance(result, RunResult)
        assert result.returncode == 0
        assert result.stdout is None  # streaming mode doesn't capture

    def test_capture_timeout(self):
        with pytest.raises(Exception):
            run_command(
                [sys.executable, "-c", "import time; time.sleep(30)"],
                capture=True,
                timeout=1,
            )
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/read_simulator/test_run_result.py -v --no-cov`
Expected: FAIL with ImportError (RunResult doesn't exist yet)

- [ ] **Step 3: Implement RunResult and extend run_command**

In `muc_one_up/read_simulator/utils/common_utils.py`:

1. Add `RunResult` dataclass after imports:

```python
from dataclasses import dataclass

@dataclass(frozen=True, slots=True)
class RunResult:
    """Result of running an external command.

    Attributes:
        returncode: Process exit code.
        stdout: Captured stdout (None when streaming mode).
        stderr: Captured stderr (None when streaming mode).
        command: The command string (for error messages).
    """

    returncode: int
    stdout: str | None
    stderr: str | None
    command: str
```

2. Extend `run_command()` signature to add `capture` parameter and return `RunResult`:

```python
def run_command(
    cmd: list[str],
    timeout: int | None = None,
    stderr_log_level: int = logging.ERROR,
    stderr_prefix: str = "",
    capture: bool = False,
    cwd: Path | None = None,
) -> RunResult:
```

3. Add capture mode implementation before the existing streaming code:

```python
    if capture:
        try:
            proc = subprocess.run(
                cmd,
                shell=False,
                capture_output=True,
                text=True,
                timeout=timeout,
                cwd=cwd,
            )
        except subprocess.TimeoutExpired as e:
            raise ExternalToolError(
                tool="command", exit_code=-1,
                stderr=f"Command timed out after {timeout}s",
                cmd=cmd_str,
            ) from e
        except Exception as e:
            raise ExternalToolError(
                tool="command", exit_code=1, stderr=str(e), cmd=cmd_str
            ) from e

        result = RunResult(
            returncode=proc.returncode,
            stdout=proc.stdout,
            stderr=proc.stderr,
            command=cmd_str,
        )
        if proc.returncode != 0:
            logging.error("Command exited with code %d: %s", proc.returncode, cmd_str)
            raise ExternalToolError(
                tool="command", exit_code=proc.returncode,
                stderr=proc.stderr or "Command failed", cmd=cmd_str
            )
        return result
```

4. Change the existing streaming path to also return `RunResult`:

Replace `return proc.returncode` with:
```python
    return RunResult(
        returncode=proc.returncode,
        stdout=None,
        stderr=None,
        command=cmd_str,
    )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/read_simulator/test_run_result.py -v --no-cov`
Expected: All PASS

- [ ] **Step 5: Run existing tests to verify backward compat**

Run: `uv run pytest tests/read_simulator/test_command_utils.py tests/read_simulator/ --tb=short -q --no-cov`
Expected: All PASS (existing callers compare `run_command() == 0` which now gets `RunResult == 0` — need to verify; if this breaks, callers that check `== 0` will need updating)

**IMPORTANT:** If existing callers break because they compare the return value to `int`, fix those callers in the same commit. The wrapper files that use `run_command()` (minimap2, pbsim3, reseq, ccs, ucsc_tools) may need updating if they check `result == 0`.

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/read_simulator/utils/common_utils.py tests/read_simulator/test_run_result.py
git commit -m "feat: add RunResult dataclass and capture mode to run_command"
```

---

## Task 2: Add run_pipeline for Pipe Chains

**Files:**
- Modify: `muc_one_up/read_simulator/utils/common_utils.py`
- Modify: `tests/read_simulator/test_run_result.py`

- [ ] **Step 1: Write failing tests for run_pipeline**

Append to `tests/read_simulator/test_run_result.py`:

```python
from muc_one_up.read_simulator.utils.common_utils import run_pipeline


class TestRunPipeline:
    """Tests for run_pipeline (pipe chains)."""

    def test_two_command_pipe(self):
        """Echo piped to grep-like filter."""
        result = run_pipeline(
            [
                [sys.executable, "-c", "print('hello\\nworld\\nfoo')"],
                [sys.executable, "-c", "import sys; [print(l, end='') for l in sys.stdin if 'world' in l]"],
            ],
            capture=True,
        )
        assert result.returncode == 0
        assert "world" in result.stdout
        assert "hello" not in result.stdout

    def test_single_command_pipeline(self):
        """Single command in pipeline works like run_command."""
        result = run_pipeline(
            [[sys.executable, "-c", "print('solo')"]],
            capture=True,
        )
        assert result.returncode == 0
        assert "solo" in result.stdout

    def test_pipeline_failure_raises(self):
        """Failure in any pipeline stage raises."""
        with pytest.raises(Exception):
            run_pipeline(
                [
                    [sys.executable, "-c", "import sys; sys.exit(1)"],
                    [sys.executable, "-c", "pass"],
                ],
                capture=True,
            )

    def test_pipeline_timeout(self):
        with pytest.raises(Exception):
            run_pipeline(
                [[sys.executable, "-c", "import time; time.sleep(30)"]],
                capture=True,
                timeout=1,
            )
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/read_simulator/test_run_result.py::TestRunPipeline -v --no-cov`
Expected: FAIL with ImportError

- [ ] **Step 3: Implement run_pipeline**

Add to `muc_one_up/read_simulator/utils/common_utils.py`:

```python
def run_pipeline(
    cmds: list[list[str]],
    *,
    capture: bool = True,
    timeout: int | None = None,
    cwd: Path | None = None,
) -> RunResult:
    """Run a pipeline of commands connected by pipes.

    Connects stdout of each process to stdin of the next.
    Returns the result of the final process.

    Args:
        cmds: List of commands, each as a list of strings.
        capture: If True, capture final stdout/stderr.
        timeout: Timeout in seconds for entire pipeline.
        cwd: Working directory for all processes.

    Returns:
        RunResult from the final process in the pipeline.

    Raises:
        ExternalToolError: If any process in the pipeline fails.
    """
    if not cmds:
        raise ValueError("Pipeline must have at least one command")

    for cmd in cmds:
        if not isinstance(cmd, list):
            raise TypeError("Each command must be a list of strings")

    pipeline_str = " | ".join(" ".join(cmd) for cmd in cmds)
    logging.info("Running pipeline: %s (timeout=%s)", pipeline_str, timeout)

    if len(cmds) == 1:
        return run_command(cmds[0], timeout=timeout, capture=capture, cwd=cwd)

    processes: list[subprocess.Popen] = []
    try:
        # Start all processes in the pipeline
        for i, cmd in enumerate(cmds):
            stdin = processes[-1].stdout if i > 0 else None
            stdout = subprocess.PIPE
            stderr = subprocess.PIPE if (i == len(cmds) - 1) else subprocess.DEVNULL

            proc = subprocess.Popen(
                cmd,
                shell=False,
                stdin=stdin,
                stdout=stdout,
                stderr=stderr,
                cwd=cwd,
            )
            processes.append(proc)

            # Allow previous process to receive SIGPIPE
            if i > 0 and processes[-2].stdout:
                processes[-2].stdout.close()

        # Wait for final process
        last = processes[-1]
        try:
            stdout_bytes, stderr_bytes = last.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            for p in processes:
                p.kill()
            last.communicate()
            raise ExternalToolError(
                tool="pipeline", exit_code=-1,
                stderr=f"Pipeline timed out after {timeout}s",
                cmd=pipeline_str,
            )

        # Wait for all upstream processes
        for p in processes[:-1]:
            p.wait()

        # Check for failures in any stage
        for i, p in enumerate(processes):
            if p.returncode != 0:
                raise ExternalToolError(
                    tool="pipeline", exit_code=p.returncode,
                    stderr=f"Pipeline stage {i} failed",
                    cmd=pipeline_str,
                )

        stdout_str = stdout_bytes.decode("utf-8", errors="replace") if capture else None
        stderr_str = stderr_bytes.decode("utf-8", errors="replace") if capture else None

        return RunResult(
            returncode=0,
            stdout=stdout_str,
            stderr=stderr_str,
            command=pipeline_str,
        )
    except ExternalToolError:
        raise
    except Exception as e:
        raise ExternalToolError(
            tool="pipeline", exit_code=1, stderr=str(e), cmd=pipeline_str
        ) from e
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/read_simulator/test_run_result.py -v --no-cov`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/read_simulator/utils/common_utils.py tests/read_simulator/test_run_result.py
git commit -m "feat: add run_pipeline for pipe chain execution"
```

---

## Task 3: Migrate utils/samtools.py to run_command

**Files:**
- Modify: `muc_one_up/read_simulator/utils/samtools.py`

- [ ] **Step 1: Read the file and identify all 9 subprocess.run calls**

Read the file completely. Each `subprocess.run()` call needs to be replaced with `run_command(cmd, capture=True)` and use `result.stdout` where output was captured.

- [ ] **Step 2: Replace all subprocess.run calls**

For each call:
- `subprocess.run(cmd, capture_output=True, text=True, check=True)` → `run_command(cmd, capture=True)`
- Access `result.stdout` instead of `proc.stdout`
- `subprocess.run(cmd, check=True, stdout=fh)` → use `run_command(cmd, capture=True)` and write `result.stdout` to file, OR use the existing streaming mode with shell redirection handled separately

For calls that redirect stdout to a file (like `stdout=fh`), use capture mode and write the captured output to file:
```python
result = run_command(cmd, capture=True)
Path(output_file).write_text(result.stdout)
```

- [ ] **Step 3: Remove subprocess import**

Remove `import subprocess` from the file.

- [ ] **Step 4: Run tests**

Run: `uv run pytest tests/read_simulator/ tests/test_samtools.py --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/read_simulator/utils/samtools.py
git commit -m "refactor: migrate utils/samtools.py to use run_command"
```

---

## Task 4: Migrate utils/bed.py and utils/tool_version.py

**Files:**
- Modify: `muc_one_up/read_simulator/utils/bed.py`
- Modify: `muc_one_up/read_simulator/utils/tool_version.py`

- [ ] **Step 1: Read both files**

- [ ] **Step 2: Replace subprocess calls in bed.py (2 calls)**

Replace `subprocess.run()` calls with `run_command(cmd, capture=True)`. For bedtools subtract which writes to a file, capture and write.

- [ ] **Step 3: Replace subprocess calls in tool_version.py (1 call)**

Replace the version query `subprocess.run()` call with `run_command(cmd, capture=True)`.

- [ ] **Step 4: Remove subprocess imports from both files**

- [ ] **Step 5: Run tests**

Run: `uv run pytest tests/read_simulator/ --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/read_simulator/utils/bed.py muc_one_up/read_simulator/utils/tool_version.py
git commit -m "refactor: migrate bed.py and tool_version.py to use run_command"
```

---

## Task 5: Migrate bwa_wrapper.py to run_command/run_pipeline

**Files:**
- Modify: `muc_one_up/read_simulator/wrappers/bwa_wrapper.py`

- [ ] **Step 1: Read the file and identify all 7 subprocess calls**

The key pattern is a pipe chain: BWA mem → samtools view (Popen piping) followed by samtools sort and samtools index (subprocess.run).

- [ ] **Step 2: Replace pipe chain with run_pipeline**

The `bwa mem | samtools view` pipe chain becomes:
```python
result = run_pipeline(
    [bwa_cmd, samtools_view_cmd],
    capture=True,
)
# Write BAM output to file
Path(output_bam).write_bytes(result.stdout.encode())
```

**Note:** BAM output is binary. If `run_pipeline` captures as text, this won't work for binary data. You may need to add a `binary=True` option to `run_pipeline` or handle the BAM pipe chain differently. Read the current implementation carefully.

If binary output is needed, keep using `Popen` directly for the BAM pipe chain but wrap it in a helper function, or add `text=False` support to `run_pipeline`.

- [ ] **Step 3: Replace samtools sort/index with run_command**

```python
run_command(sort_cmd, capture=True)
run_command(index_cmd, capture=True)
```

- [ ] **Step 4: Remove direct subprocess import**

- [ ] **Step 5: Run tests**

Run: `uv run pytest tests/read_simulator/test_bwa_wrapper.py --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/read_simulator/wrappers/bwa_wrapper.py
git commit -m "refactor: migrate bwa_wrapper.py to run_command/run_pipeline"
```

---

## Task 6: Migrate samtools_wrapper.py to run_command/run_pipeline

**Files:**
- Modify: `muc_one_up/read_simulator/wrappers/samtools_wrapper.py`

- [ ] **Step 1: Read the file and identify all 11 subprocess calls**

10 can be migrated. 1 (the conda `shell=True` exception at line ~806) stays as-is.

- [ ] **Step 2: Replace subprocess.run calls with run_command**

For calls that write stdout to file:
```python
result = run_command(cmd, capture=True)
Path(output_file).write_text(result.stdout)
```

- [ ] **Step 3: Replace Popen pipe chains with run_pipeline**

The `samtools collate | samtools fastq` pipe chain:
```python
result = run_pipeline(
    [collate_cmd, fastq_cmd],
    capture=True,
)
```

- [ ] **Step 4: Keep conda shell=True exception**

The conda/mamba workaround at ~line 806 stays unchanged (properly secured with shlex.quote, marked with nosec B602).

- [ ] **Step 5: Remove subprocess import (except for the conda exception)**

Keep `import subprocess` since the conda exception still uses it directly.

- [ ] **Step 6: Run tests**

Run: `uv run pytest tests/read_simulator/test_samtools_wrapper.py --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 7: Commit**

```bash
git add muc_one_up/read_simulator/wrappers/samtools_wrapper.py
git commit -m "refactor: migrate samtools_wrapper.py to run_command/run_pipeline (except conda exception)"
```

---

## Task 7: Migrate nanosim_wrapper.py

**Files:**
- Modify: `muc_one_up/read_simulator/wrappers/nanosim_wrapper.py`

- [ ] **Step 1: Read the file and identify the 1 subprocess call**

- [ ] **Step 2: Replace with run_command**

The SAM→BAM conversion subprocess.run becomes:
```python
result = run_command(cmd, capture=True)
```
If the call writes stdout to a file, capture and write.

- [ ] **Step 3: Remove subprocess import**

- [ ] **Step 4: Run tests**

Run: `uv run pytest tests/read_simulator/test_nanosim_wrapper.py --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/read_simulator/wrappers/nanosim_wrapper.py
git commit -m "refactor: migrate nanosim_wrapper.py to run_command"
```

---

## Task 8: Verify Full Suite and Success Criteria

**Files:** None (verification only)

- [ ] **Step 1: Run full test suite**

Run: `uv run pytest --tb=short -q`
Expected: All tests pass.

- [ ] **Step 2: Run ruff linter and formatter**

Run: `uv run ruff check muc_one_up/ tests/ && uv run ruff format --check muc_one_up/ tests/`
Expected: Clean.

- [ ] **Step 3: Run mypy**

Run: `uv run mypy muc_one_up/`
Expected: No new errors.

- [ ] **Step 4: Verify subprocess centralization**

```bash
# Only common_utils.py and samtools_wrapper.py (conda exception) should import subprocess
grep -rn "import subprocess" muc_one_up/read_simulator/ | grep -v __pycache__
# Expected: only common_utils.py and samtools_wrapper.py (conda shell=True exception)
```

- [ ] **Step 5: Commit any fixes**

```bash
git add -u
git commit -m "style: fix formatting after Phase 4A subprocess migration"
```

- [ ] **Step 6: Final test run**

Run: `uv run pytest --tb=short -q`
Expected: All pass. Phase 4A is complete.
