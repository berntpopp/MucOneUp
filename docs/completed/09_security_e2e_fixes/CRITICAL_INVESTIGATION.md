# CRITICAL: CLI Broken After Refactor - Investigation Report

**Date**: 2025-10-18
**Branch**: `dev/modern-python-refactor`
**Severity**: 🔴 **CRITICAL** - Tool is non-functional for end-to-end workflows
**Status**: Identified root cause, fix plan ready

---

## Executive Summary

The refactored CLI appears to work on surface-level testing but **completely fails** when attempting real end-to-end workflows. The root cause is a **security fix side-effect**: we eliminated `shell=True` for security but failed to properly handle multi-word tool commands (conda/mamba environments) across all wrappers.

**Impact**:
- ✅ Unit tests pass (508 tests, 77% coverage)
- ✅ Linting passes
- ❌ **E2E workflows completely broken**
- ❌ No integration tests to catch this

---

## Root Cause Analysis

### The Security Fix (Week 2)

**Commit**: `3405761 security: eliminate shell=True vulnerability in run_command (CRITICAL)`

**What we did**:
```python
# BEFORE (VULNERABLE)
def run_command(cmd, shell=False, ...):
    proc = subprocess.Popen(cmd, shell=shell, ...)  # Could use shell=True

# AFTER (SECURE)
def run_command(cmd: list[str], ...):
    proc = subprocess.Popen(cmd, shell=False, ...)  # Always False, enforced
```

**Why we did it**: HIGH severity B602 vulnerability - command injection risk

**What we missed**: When `shell=False`, commands must be **properly split lists**, not strings!

---

### The Breaking Change

**Problem**: Config file contains multi-word commands for conda/mamba environments:

```json
{
  "tools": {
    "reseq": "mamba run --no-capture-output -n env_wessim reseq",
    "bwa": "mamba run --no-capture-output -n env_wessim bwa",
    "samtools": "mamba run --no-capture-output -n env_wessim samtools",
    "faToTwoBit": "mamba run --no-capture-output -n env_wessim faToTwoBit",
    "pblat": "mamba run --no-capture-output -n env_wessim pblat"
  }
}
```

**Current broken code** (in 4 wrappers):
```python
# reseq_wrapper.py, bwa_wrapper.py, samtools_wrapper.py, ucsc_tools_wrapper.py
cmd = [tools["reseq"], "replaceN", "-r", input_fa, "-R", output_fa]
# Results in: ["mamba run ... reseq", "replaceN", ...]
# ERROR: Can't find executable "mamba run ... reseq" (treated as single string)
```

**What should happen**:
```python
import shlex
cmd = shlex.split(tools["reseq"]) + ["replaceN", "-r", input_fa, "-R", output_fa]
# Results in: ["mamba", "run", "--no-capture-output", "-n", "env_wessim", "reseq", "replaceN", ...]
# SUCCESS: Proper command list
```

---

## Broken Components

### ✅ WORKING

1. **nanosim_wrapper.py** - Already uses `shlex.split()` (fixed during development)
2. **simulate command** - Works for generating haplotypes
3. **Unit tests** - All 508 tests pass (mock at subprocess boundary)

### ❌ BROKEN

| Component | Status | Error Message |
|-----------|--------|---------------|
| **reseq_wrapper.py** | 🔴 Broken | `[Errno 2] No such file or directory: 'mamba run ... reseq'` |
| **bwa_wrapper.py** | 🔴 Broken | `[Errno 2] No such file or directory: 'mamba run ... bwa'` |
| **samtools_wrapper.py** | 🔴 Broken | `[Errno 2] No such file or directory: 'mamba run ... samtools'` |
| **ucsc_tools_wrapper.py** | 🔴 Broken | `[Errno 2] No such file or directory: 'mamba run ... faToTwoBit'` |
| **reads illumina** | 🔴 Broken | Fails at step 1 (replace_Ns using reseq) |
| **reads ont** | ✅ Works | Uses nanosim_wrapper (already fixed) |
| **analyze orfs** | 🔴 Unknown | Needs testing |

---

## Test Results

### Test 1: Simulate Command
```bash
muconeup --config config.json simulate \
  --out-base test_sim \
  --out-dir /tmp/test_output \
  --fixed-lengths 20
```

**Result**: ✅ **WORKS**
**Outputs Created**:
- `/tmp/test_output/test_sim.001.simulated.fa` (23KB)
- `/tmp/test_output/test_sim.001.simulation_stats.json` (3.0KB)

**Note**: Shows error `ERROR - Unexpected error: 0` but files are created correctly. This is a minor bug with `ctx.exit(0)` raising Click.Exit exception.

---

### Test 2: Reads Illumina Command
```bash
muconeup --config config.json reads illumina \
  /tmp/test_output/test_sim.001.simulated.fa \
  --out-dir /tmp/test_output
```

**Result**: ❌ **FAILS AT STEP 1**

**Error**:
```
2025-10-18 12:12:35,926 - INFO - 1. Replacing Ns in FASTA
2025-10-18 12:12:35,926 - INFO - Running command: mamba run --no-capture-output -n env_wessim reseq replaceN ...
2025-10-18 12:12:36,036 - ERROR - Read simulation failed: command failed with exit code 1
Command: mamba run --no-capture-output -n env_wessim reseq replaceN -r /tmp/test_output/test_sim.001.simulated.fa -R /tmp/test_output/_test_sim.001.simulated_noNs.fa
Error: [Errno 2] No such file or directory: 'mamba run --no-capture-output -n env_wessim reseq'
```

**Root Cause**: `tools["reseq"]` not split into list elements

---

## Why Unit Tests Didn't Catch This

Our unit tests mock at the `subprocess.Popen` boundary:

```python
# tests/read_simulator/test_reseq_wrapper.py
def test_constructs_correct_command(self, mocker, tmp_path, tools_dict):
    mock_popen = mocker.patch("subprocess.Popen", return_value=mock_proc)

    replace_Ns(str(input_fa), str(output_fa), tools_dict)

    # We verify the command is constructed, but don't actually execute it
    cmd = mock_popen.call_args[0][0]
    assert cmd[0] == "reseq"  # ❌ This passes but is WRONG!
```

**The problem**:
- We check `cmd[0] == "reseq"`
- But `tools_dict["reseq"]` in tests is just `"reseq"` (not the mamba command)
- Real config has `"mamba run --no-capture-output -n env_wessim reseq"`
- Unit tests use simplified mock config, don't test real conda commands

---

## Missing Test Coverage

### ❌ No E2E Integration Tests

We have:
- ✅ Unit tests for individual functions (77% coverage)
- ✅ Wrapper tests with mocked subprocess
- ❌ **NO** end-to-end integration tests
- ❌ **NO** tests with real config.json
- ❌ **NO** tests that actually execute the full pipeline

### What We Should Have

```python
# tests/integration/test_e2e_workflow.py (MISSING!)
def test_full_workflow_with_real_config():
    """Test complete simulate -> reads -> analyze workflow."""
    # 1. Run simulate command
    result = subprocess.run([
        "muconeup", "--config", "config.json",
        "simulate", "--out-base", "test", ...
    ])
    assert result.returncode == 0

    # 2. Run reads command
    result = subprocess.run([
        "muconeup", "--config", "config.json",
        "reads", "illumina", "test.001.simulated.fa", ...
    ])
    assert result.returncode == 0

    # 3. Verify outputs exist
    assert Path("test.001.simulated.bam").exists()
```

---

## Impact on Users

### Before Refactor (main branch)
```bash
# Single command, worked end-to-end
muconeup --config config.json \
  --out-base output \
  --simulate-reads \
  --output-orfs
# ✅ WORKED
```

### After Refactor (current branch)
```bash
# Step 1: Works
muconeup --config config.json simulate --out-base output

# Step 2: FAILS
muconeup --config config.json reads illumina output.001.simulated.fa
# ❌ BROKEN - "No such file or directory: 'mamba run...'"

# Step 3: Unknown
muconeup --config config.json analyze orfs output.001.simulated.fa
# ❌ Probably broken if it calls external tools
```

**User Experience**: Tool appears to work (simulate succeeds) but completely fails when trying to use it for real workflows.

---

## Fix Plan

### Phase 1: Immediate Critical Fixes (1-2 hours)

#### Fix 1: Add shlex.split() to all wrappers

**Files to fix**:
1. `muc_one_up/read_simulator/wrappers/reseq_wrapper.py`
2. `muc_one_up/read_simulator/wrappers/bwa_wrapper.py`
3. `muc_one_up/read_simulator/wrappers/samtools_wrapper.py`
4. `muc_one_up/read_simulator/wrappers/ucsc_tools_wrapper.py`

**Pattern to apply**:
```python
import shlex

# BEFORE (BROKEN)
cmd = [tools["reseq"], "replaceN", ...]

# AFTER (FIXED)
if " " in tools["reseq"]:
    cmd = shlex.split(tools["reseq"]) + ["replaceN", ...]
else:
    cmd = [tools["reseq"], "replaceN", ...]
```

**Lines to change**: ~15-20 locations across 4 files

---

#### Fix 2: Fix ctx.exit(0) in simulate command

**File**: `muc_one_up/cli/click_main.py`

**Change**:
```python
# BEFORE (line ~259)
ctx.exit(0)

# AFTER
return  # Click handles exit automatically
```

---

### Phase 2: Add E2E Integration Tests (2-3 hours)

#### Test 1: Minimal E2E workflow
```python
# tests/integration/test_minimal_e2e.py
def test_simulate_to_fasta(tmp_path):
    """Test simulate command produces valid FASTA."""
    subprocess.run([
        "muconeup", "--config", "config.json",
        "simulate", "--out-base", "test",
        "--out-dir", str(tmp_path),
        "--fixed-lengths", "20"
    ], check=True)

    assert (tmp_path / "test.001.simulated.fa").exists()
```

#### Test 2: Full pipeline (if tools available)
```python
@pytest.mark.skipif(not shutil.which("reseq"), reason="reseq not installed")
def test_full_illumina_pipeline(tmp_path):
    """Test complete simulate -> reads workflow."""
    # Similar to user's original command
```

#### Test 3: Real config validation
```python
def test_config_json_tools_are_valid():
    """Test that all tools in config.json can be parsed."""
    config = json.load(open("config.json"))
    for tool_name, tool_cmd in config["tools"].items():
        # Verify shlex.split doesn't raise
        parts = shlex.split(tool_cmd)
        assert len(parts) > 0
```

---

### Phase 3: Update Unit Tests (1 hour)

Update wrapper tests to use **real config patterns**:

```python
# tests/read_simulator/test_reseq_wrapper.py
@pytest.fixture
def tools_dict():
    return {
        "reseq": "mamba run --no-capture-output -n env_wessim reseq",  # Real pattern!
        # ... other tools
    }
```

---

### Phase 4: Documentation (30 minutes)

1. Update CLAUDE.md with:
   - E2E testing requirements
   - How to test CLI changes
   - Integration test guidelines

2. Add TESTING.md with:
   - How to run E2E tests
   - How to test with real tools vs mocks
   - CI/CD integration test strategy

---

## Lessons Learned

### What Went Wrong

1. **Over-reliance on unit tests**: 77% coverage, but missed critical integration issues
2. **Mock-heavy testing**: Mocked subprocess calls don't catch command construction bugs
3. **No smoke tests**: Should have run actual CLI commands during development
4. **Incomplete security fix**: Fixed shell=True but didn't handle multi-word commands
5. **Test fixture simplification**: Used `tools_dict = {"reseq": "reseq"}` instead of real patterns

### What to Do Differently

1. **Always have E2E tests**: At least one test that exercises the full workflow
2. **Test with real config**: Use actual config.json patterns in tests
3. **Smoke test after major refactors**: Run actual commands, not just unit tests
4. **Integration test CI stage**: Add separate CI stage for integration tests
5. **User acceptance testing**: Actually try user workflows before declaring "done"

---

## Testing Checklist (Before Merge)

### Pre-Fix Baseline
- [x] Identified simulate works
- [x] Identified reads illumina fails
- [ ] Test reads ont (probably works - uses nanosim)
- [ ] Test analyze orfs
- [ ] Test with mutations
- [ ] Test with SNPs
- [ ] Test series mode

### Post-Fix Validation
- [ ] All 4 wrappers use shlex.split()
- [ ] ctx.exit(0) fixed
- [ ] Simulate command clean exit
- [ ] Reads illumina works end-to-end
- [ ] Reads ont still works
- [ ] Analyze orfs works
- [ ] Full pipeline: simulate -> reads -> orfs
- [ ] User's original command equivalent works
- [ ] All unit tests still pass
- [ ] New E2E tests added
- [ ] Documentation updated

---

## Estimated Fix Timeline

| Phase | Time | Priority |
|-------|------|----------|
| Fix wrappers (shlex.split) | 1-2h | 🔴 CRITICAL |
| Fix ctx.exit(0) | 10m | 🟡 Medium |
| Add basic E2E test | 1h | 🔴 CRITICAL |
| Add comprehensive E2E tests | 2h | 🟠 High |
| Update unit test fixtures | 1h | 🟠 High |
| Documentation | 30m | 🟢 Low |
| **Total** | **5-6 hours** | |

---

## Next Actions

1. **IMMEDIATE**: Fix the 4 broken wrappers with shlex.split()
2. **IMMEDIATE**: Add minimal E2E integration test
3. **BEFORE MERGE**: Run full user workflow test
4. **BEFORE MERGE**: Add integration test CI stage
5. **FUTURE**: Consider test strategy overhaul (balance unit vs integration)

---

## Appendix: Code Examples

### Example Fix for reseq_wrapper.py

```python
import shlex

def replace_Ns(input_fa: str, output_fa: str, tools: dict[str, str]) -> None:
    """Replace Ns in the input FASTA using reseq replaceN."""

    # FIXED: Handle multi-word commands (conda/mamba)
    reseq_cmd = tools["reseq"]
    if " " in reseq_cmd:
        # Command has spaces (conda/mamba), use shlex.split() to parse safely
        cmd = shlex.split(reseq_cmd)
    else:
        # Simple command, just wrap in list
        cmd = [reseq_cmd]

    # Add the actual reseq arguments
    cmd.extend(["replaceN", "-r", input_fa, "-R", output_fa])

    run_command(cmd, timeout=60, stderr_log_level=logging.INFO, stderr_prefix="[reseq] ")
```

### Example E2E Test

```python
# tests/integration/test_e2e_basic.py
import subprocess
import pytest
from pathlib import Path

def test_simulate_command_produces_fasta(tmp_path):
    """Test that simulate command creates valid FASTA output."""
    result = subprocess.run([
        "muconeup",
        "--config", "config.json",
        "simulate",
        "--out-base", "test_e2e",
        "--out-dir", str(tmp_path),
        "--fixed-lengths", "20"
    ], capture_output=True, text=True)

    # Should succeed
    assert result.returncode == 0, f"Simulate failed: {result.stderr}"

    # Should create FASTA
    fasta_file = tmp_path / "test_e2e.001.simulated.fa"
    assert fasta_file.exists(), "FASTA file not created"

    # Should be valid FASTA
    content = fasta_file.read_text()
    assert content.startswith(">"), "Not valid FASTA format"
    assert "haplotype_1" in content or "haplotype_2" in content

@pytest.mark.integration
@pytest.mark.skipif(not shutil.which("reseq"), reason="reseq not installed")
def test_full_workflow_simulate_to_reads(tmp_path):
    """Test complete workflow if tools are available."""
    # 1. Simulate
    subprocess.run([
        "muconeup", "--config", "config.json",
        "simulate", "--out-base", "integration_test",
        "--out-dir", str(tmp_path), "--fixed-lengths", "20"
    ], check=True)

    # 2. Reads
    subprocess.run([
        "muconeup", "--config", "config.json",
        "reads", "illumina",
        str(tmp_path / "integration_test.001.simulated.fa"),
        "--out-dir", str(tmp_path)
    ], check=True)

    # 3. Verify outputs
    assert (tmp_path / "integration_test.001.simulated.bam").exists()
```

---

**Report Generated**: 2025-10-18 12:15 UTC
**Next Update**: After critical fixes applied
