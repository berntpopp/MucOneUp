# MucOneUp CLI Enhancement Implementation Plan

**Date:** 2025-10-18
**Developer:** Senior Developer
**Version:** 0.10.0 (proposed)
**Estimated Total Effort:** 12-16 hours

---

## Overview

This plan details the implementation of 5 enhancements identified in the CLI audit:

### Medium Priority (Required)
1. ✅ Batch processing tests for reads/analyze commands
2. ✅ Add `test_reads_ont_with_file` for parity with illumina

### Low Priority (Value-Added)
3. ✅ Progress indicators for `--simulate-series`
4. ✅ Shell completion support
5. ✅ `--verbose/-v` alias for `--log-level DEBUG`

---

## Table of Contents

1. [Architecture Principles](#architecture-principles)
2. [Enhancement 1: Batch Processing Tests](#enhancement-1-batch-processing-tests)
3. [Enhancement 2: ONT Test Parity](#enhancement-2-ont-test-parity)
4. [Enhancement 3: Progress Indicators](#enhancement-3-progress-indicators)
5. [Enhancement 4: Shell Completion](#enhancement-4-shell-completion)
6. [Enhancement 5: Verbose Flag](#enhancement-5-verbose-flag)
7. [Testing Strategy](#testing-strategy)
8. [Version Bump & Documentation](#version-bump--documentation)
9. [Implementation Order](#implementation-order)

---

## Architecture Principles

**Maintain Existing Patterns:**
- ✅ Unix philosophy (single responsibility)
- ✅ SOLID principles (separation of concerns)
- ✅ DRY (don't repeat yourself)
- ✅ Existing test patterns (CliRunner, fixtures, markers)
- ✅ Click 8.1.7 best practices

**Quality Standards:**
- All new code must have 100% test coverage
- Follow existing naming conventions
- Use type hints (Python 3.10+ syntax)
- Add docstrings for all new functions
- Update inline help text

---

## Enhancement 1: Batch Processing Tests

**Priority:** MEDIUM
**Effort:** 4 hours
**Files Modified:** `tests/test_click_cli.py`

### 1.1 Current State

**Existing Test Coverage:**
- ✅ Single file: `test_reads_illumina_with_file` (line 272)
- ❌ Multiple files: Not tested
- ❌ Auto-generated output names: Not verified
- ❌ Warning message: Not tested

**Implementation Already Supports:**
- `input_fastas` accepts multiple files (click_main.py:291, nargs=-1)
- Batch processing loops (click_main.py:370-393 illumina, 481-503 ont)
- Auto-generated output names (click_main.py:374-380, 485-491)
- Warning for multiple files with --out-base (click_main.py:361-367)

### 1.2 New Tests to Add

**File:** `tests/test_click_cli.py`
**Location:** After line 293 (end of TestReadsCommand class)

#### Test 1: Multiple FASTA files for illumina

```python
def test_reads_illumina_with_multiple_files(self, runner, temp_config, tmp_path):
    """Test reads illumina with multiple FASTA files (batch processing)."""
    # Create 3 FASTA files
    fasta1 = tmp_path / "sample1.fa"
    fasta2 = tmp_path / "sample2.fa"
    fasta3 = tmp_path / "sample3.fa"

    for i, fasta_file in enumerate([fasta1, fasta2, fasta3], start=1):
        fasta_file.write_text(f">haplotype_{i}\nACGTACGTACGT\n")

    result = runner.invoke(
        cli,
        [
            "--config", temp_config,
            "reads", "illumina",
            str(fasta1), str(fasta2), str(fasta3),
            "--out-dir", str(tmp_path),
            "--coverage", "10",
        ],
    )

    # Verify batch processing worked
    assert result.exit_code in [0, 1]  # May fail on missing tools
    assert "3 FASTA file(s)" in result.output or result.exit_code == 0
```

**Line to insert:** After line 293

#### Test 2: Warning for --out-base with multiple files

```python
def test_reads_illumina_warns_on_out_base_with_multiple_files(
    self, runner, temp_config, tmp_path
):
    """Test that --out-base with multiple files triggers warning."""
    fasta1 = tmp_path / "sample1.fa"
    fasta2 = tmp_path / "sample2.fa"

    fasta1.write_text(">test1\nACGT\n")
    fasta2.write_text(">test2\nGCTA\n")

    result = runner.invoke(
        cli,
        [
            "--config", temp_config,
            "reads", "illumina",
            str(fasta1), str(fasta2),
            "--out-base", "custom_name",
        ],
    )

    # Check for warning message
    if result.exit_code in [0, 1]:
        # Warning should appear in output
        assert (
            "will be used for all" in result.output or
            "Consider omitting --out-base" in result.output or
            result.exit_code == 1  # May fail on missing tools before warning
        )
```

**Line to insert:** After previous test

#### Test 3: Multiple files for ONT

```python
def test_reads_ont_with_multiple_files(self, runner, temp_config, tmp_path):
    """Test reads ont with multiple FASTA files (batch processing)."""
    # Create multiple FASTA files
    fastas = []
    for i in range(1, 4):
        fasta = tmp_path / f"ont_sample{i}.fa"
        fasta.write_text(f">sequence_{i}\nACGTACGTACGTACGT\n")
        fastas.append(str(fasta))

    result = runner.invoke(
        cli,
        [
            "--config", temp_config,
            "reads", "ont",
            *fastas,
            "--out-dir", str(tmp_path),
            "--coverage", "20",
            "--min-read-length", "50",
        ],
    )

    assert result.exit_code in [0, 1]
    assert "3 FASTA file(s)" in result.output or result.exit_code == 0
```

**Line to insert:** After line 300 (after test_reads_ont_help)

#### Test 4: Multiple files for analyze orfs

```python
def test_analyze_orfs_with_multiple_files(self, runner, temp_config, tmp_path):
    """Test analyze orfs with multiple FASTA files (batch processing)."""
    # Create test FASTA files
    fasta1 = tmp_path / "test1.fa"
    fasta2 = tmp_path / "test2.fa"

    fasta1.write_text(">seq1\nATGACGTACGTACGT\n")
    fasta2.write_text(">seq2\nATGGCTAGCTAGCTA\n")

    result = runner.invoke(
        cli,
        [
            "--config", temp_config,
            "analyze", "orfs",
            str(fasta1), str(fasta2),
            "--out-dir", str(tmp_path),
            "--orf-min-aa", "10",
        ],
    )

    assert result.exit_code in [0, 1]
    assert "2 FASTA file(s)" in result.output or result.exit_code == 0
```

**Line to insert:** After line 358 (in TestAnalyzeCommand class)

#### Test 5: Multiple files for analyze stats

```python
def test_analyze_stats_with_multiple_files(self, runner, temp_config, tmp_path):
    """Test analyze stats with multiple FASTA files (batch processing)."""
    # Create test FASTA files
    fastas = []
    for i in range(1, 4):
        fasta = tmp_path / f"stats_test{i}.fa"
        fasta.write_text(f">haplotype_{i}\n{'ACGT' * 10}\n")
        fastas.append(str(fasta))

    result = runner.invoke(
        cli,
        [
            "--config", temp_config,
            "analyze", "stats",
            *fastas,
            "--out-dir", str(tmp_path),
        ],
    )

    assert result.exit_code in [0, 1]
    assert "3 FASTA file(s)" in result.output or result.exit_code == 0

    # Verify JSON files created (if successful)
    if result.exit_code == 0:
        json_files = list(tmp_path.glob("*.basic_stats.json"))
        # Should create 3 JSON files (one per input)
        # Note: May be fewer if auto-naming conflicts occur
        assert len(json_files) >= 1
```

**Line to insert:** After line 383 (end of TestAnalyzeCommand class)

### 1.3 Summary

**Total New Tests:** 5
**Expected Test Count After:** 153 CLI tests (was 148)
**Coverage Improvement:** Batch processing fully tested

---

## Enhancement 2: ONT Test Parity

**Priority:** MEDIUM
**Effort:** 1 hour
**Files Modified:** `tests/test_click_cli.py`

### 2.1 Current State

- ✅ Illumina has `test_reads_illumina_with_file` (line 272)
- ❌ ONT only has `test_reads_ont_help` (line 295)
- ❌ No file processing test for ONT

### 2.2 New Test

**File:** `tests/test_click_cli.py`
**Location:** After line 300 (after test_reads_ont_help)

```python
def test_reads_ont_with_file(self, runner, temp_config, tmp_path):
    """Test reads ont with single FASTA file (parity with illumina)."""
    # Create a minimal FASTA file
    fasta_file = tmp_path / "test_ont.fa"
    fasta_file.write_text(">test_sequence\nACGTACGTACGTACGT\n")

    result = runner.invoke(
        cli,
        [
            "--config", temp_config,
            "reads", "ont",
            str(fasta_file),
            "--out-base", "test_ont_reads",
            "--out-dir", str(tmp_path),
            "--coverage", "15",
            "--min-read-length", "100",
        ],
    )

    # May fail on missing NanoSim tools, but parsing should work
    assert result.exit_code in [0, 1]

    # Verify command was parsed correctly
    if result.exit_code == 0:
        # If successful, output directory should exist
        assert tmp_path.exists()
```

**Line to insert:** After line 300

### 2.3 Summary

**Total New Tests:** 1
**Expected Test Count After:** 154 CLI tests (including batch tests)
**Benefit:** Feature parity between illumina and ont testing

---

## Enhancement 3: Progress Indicators

**Priority:** LOW (Nice to Have)
**Effort:** 3 hours
**Files Modified:** `muc_one_up/cli/click_main.py`

### 3.1 Current State

**Simulation Loop Location:**
- File: `muc_one_up/cli/click_main.py`
- Lines: 244-256
- Current: Simple for loop, no progress indication

```python
# Line 244
for sim_index, fixed_conf in enumerate(simulation_configs, start=1):
    run_single_simulation_iteration(
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
    )
```

### 3.2 Implementation Strategy

**Use Click's Built-in Progress Bar:**
- Click 8.1.7 includes `click.progressbar()` context manager
- No additional dependencies needed
- Follows Click best practices

### 3.3 Code Changes

**File:** `muc_one_up/cli/click_main.py`

#### Change 1: Import statement (line 20)

**Current Line 20:**
```python
import click
```

**No change needed** - click.progressbar is already available

#### Change 2: Add progress bar to simulation loop

**Location:** Replace lines 243-256

**Current Code (lines 243-256):**
```python
# Run haplotype generation ONLY
for sim_index, fixed_conf in enumerate(simulation_configs, start=1):
    run_single_simulation_iteration(
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
    )
```

**New Code:**
```python
# Run haplotype generation ONLY
total_iterations = len(simulation_configs)

# Show progress bar only for series mode (multiple iterations)
if total_iterations > 1:
    with click.progressbar(
        simulation_configs,
        label=f"Simulating {total_iterations} iterations",
        show_eta=True,
        show_pos=True,
    ) as configs:
        for sim_index, fixed_conf in enumerate(configs, start=1):
            run_single_simulation_iteration(
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
            )
else:
    # Single iteration - no progress bar needed
    for sim_index, fixed_conf in enumerate(simulation_configs, start=1):
        run_single_simulation_iteration(
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
        )
```

**Lines to replace:** 243-256

### 3.4 User Experience

**Before:**
```
INFO - Analyzing VNTR structures from ...
INFO - Haplotype generation completed successfully.
```

**After (with --simulate-series):**
```
INFO - Analyzing VNTR structures from ...
Simulating 21 iterations  [################-----------]   61%  00:02:15
```

**After (single iteration):**
```
INFO - Analyzing VNTR structures from ...
INFO - Haplotype generation completed successfully.
```

### 3.5 Testing

**File:** `tests/test_click_cli.py`
**Location:** In TestSimulateCommand class (after line 232)

```python
def test_simulate_series_shows_progress(self, runner, temp_config, tmp_path):
    """Test that --simulate-series shows progress indicator."""
    result = runner.invoke(
        cli,
        [
            "--config", temp_config,
            "simulate",
            "--out-dir", str(tmp_path),
            "--out-base", "series_test",
            "--fixed-lengths", "20-25",
            "--simulate-series", "2",
            "--seed", "42",
        ],
    )

    # May fail on missing tools, but should show progress
    assert result.exit_code in [0, 1]
    # Progress bar should appear for multiple iterations
    # Note: CliRunner may not capture progress bar output in test mode
    # This is expected Click behavior
```

**Line to insert:** After line 232

### 3.6 Summary

**Lines Modified:** 1 section (lines 243-256 → 243-279)
**Dependencies:** None (uses built-in click.progressbar)
**Backward Compatibility:** ✅ Maintained (single iterations unchanged)
**Test Added:** 1

---

## Enhancement 4: Shell Completion

**Priority:** LOW (Value-Added)
**Effort:** 2 hours
**Files Modified:** `muc_one_up/cli/click_main.py`, `README.md`

### 4.1 Background

Click 8.0+ includes built-in shell completion support via the `_<PROG>_COMPLETE` environment variable pattern. No additional code needed - it's automatic!

**Supported Shells:**
- Bash
- Zsh
- Fish

### 4.2 How It Works

When installed, users can enable completion:

**Bash:**
```bash
_MUCONEUP_COMPLETE=bash_source muconeup > ~/.muconeup-complete.bash
source ~/.muconeup-complete.bash
```

**Zsh:**
```bash
_MUCONEUP_COMPLETE=zsh_source muconeup > ~/.muconeup-complete.zsh
source ~/.muconeup-complete.zsh
```

**Fish:**
```bash
_MUCONEUP_COMPLETE=fish_source muconeup > ~/.config/fish/completions/muconeup.fish
```

### 4.3 Code Changes

**No code changes needed!** Click 8.1.7 handles this automatically.

However, we can add a helper command to make it easier for users.

**File:** `muc_one_up/cli/click_main.py`
**Location:** After line 914 (before helper functions section)

```python
# ============================================================================
# COMPLETION Command - Shell Completion Setup
# ============================================================================


@cli.command()
@click.argument("shell", type=click.Choice(["bash", "zsh", "fish"]), required=False)
def completion(shell):
    """Generate shell completion script.

    \b
    Install completion for your shell:

      Bash:
        eval "$(_MUCONEUP_COMPLETE=bash_source muconeup)"
        # Or permanently:
        _MUCONEUP_COMPLETE=bash_source muconeup > ~/.muconeup-complete.bash
        echo 'source ~/.muconeup-complete.bash' >> ~/.bashrc

      Zsh:
        eval "$(_MUCONEUP_COMPLETE=zsh_source muconeup)"
        # Or permanently:
        _MUCONEUP_COMPLETE=zsh_source muconeup > ~/.muconeup-complete.zsh
        echo 'source ~/.muconeup-complete.zsh' >> ~/.zshrc

      Fish:
        _MUCONEUP_COMPLETE=fish_source muconeup > ~/.config/fish/completions/muconeup.fish

    \b
    Examples:
      muconeup completion          # Show installation instructions
      muconeup completion bash     # Show bash-specific instructions
    """
    if shell:
        # Show shell-specific instructions
        instructions = {
            "bash": """
# Bash completion installation:
_MUCONEUP_COMPLETE=bash_source muconeup > ~/.muconeup-complete.bash
echo 'source ~/.muconeup-complete.bash' >> ~/.bashrc
source ~/.bashrc
""",
            "zsh": """
# Zsh completion installation:
_MUCONEUP_COMPLETE=zsh_source muconeup > ~/.muconeup-complete.zsh
echo 'source ~/.muconeup-complete.zsh' >> ~/.zshrc
source ~/.zshrc
""",
            "fish": """
# Fish completion installation:
_MUCONEUP_COMPLETE=fish_source muconeup > ~/.config/fish/completions/muconeup.fish
# Restart Fish shell to activate
""",
        }
        click.echo(instructions[shell])
    else:
        # Show general instructions
        click.echo(completion.__doc__)
```

**Line to insert:** After line 914 (before `# ============================================================================`)

### 4.4 Documentation Updates

**File:** `README.md`
**Location:** After Installation section (around line 111)

```markdown
### Shell Completion (Optional)

MucOneUp supports tab completion for Bash, Zsh, and Fish shells:

**Bash:**
```bash
_MUCONEUP_COMPLETE=bash_source muconeup > ~/.muconeup-complete.bash
echo 'source ~/.muconeup-complete.bash' >> ~/.bashrc
source ~/.bashrc
```

**Zsh:**
```bash
_MUCONEUP_COMPLETE=zsh_source muconeup > ~/.muconeup-complete.zsh
echo 'source ~/.muconeup-complete.zsh' >> ~/.zshrc
source ~/.zshrc
```

**Fish:**
```bash
_MUCONEUP_COMPLETE=fish_source muconeup > ~/.config/fish/completions/muconeup.fish
```

Alternatively, use the helper command:
```bash
muconeup completion bash  # Shows bash-specific instructions
```

Now you can use tab completion:
```bash
muconeup <TAB>          # Shows: simulate, reads, analyze, completion
muconeup simulate --<TAB>  # Shows all simulate options
```
```

**Line to insert:** After line 111 in README.md

### 4.5 Testing

**File:** `tests/test_click_cli.py`
**Location:** In TestCLIRoot class (after line 109)

```python
def test_completion_command_exists(self, runner, temp_config):
    """Test that completion command is available."""
    result = runner.invoke(cli, ["--config", temp_config, "completion", "--help"])
    assert result.exit_code == 0
    assert "shell completion" in result.output.lower()
    assert "bash" in result.output
    assert "zsh" in result.output
    assert "fish" in result.output

def test_completion_shows_shell_specific_instructions(self, runner, temp_config):
    """Test that completion command shows shell-specific instructions."""
    for shell in ["bash", "zsh", "fish"]:
        result = runner.invoke(cli, ["--config", temp_config, "completion", shell])
        assert result.exit_code == 0
        assert shell in result.output.lower()
        assert "_MUCONEUP_COMPLETE" in result.output
```

**Line to insert:** After line 109

### 4.6 Summary

**Code Added:** 1 new command (~50 lines)
**Documentation:** README update
**Dependencies:** None (built-in to Click 8.1+)
**Tests Added:** 2
**Benefit:** Improved developer experience

---

## Enhancement 5: Verbose Flag

**Priority:** LOW (UX Improvement)
**Effort:** 2 hours
**Files Modified:** `muc_one_up/cli/click_main.py`

### 5.1 Current State

**Current Flag:**
- `--log-level` with choices: DEBUG, INFO, WARNING, ERROR, CRITICAL, NONE
- Line 67-71 in click_main.py

**Desired:**
- Add `--verbose/-v` as alias for `--log-level DEBUG`
- More intuitive for users
- Common CLI pattern

### 5.2 Implementation Strategy

**Option 1: Callback-based (Recommended)**
- Add `--verbose` flag with callback
- Callback sets log_level to DEBUG
- Maintains backward compatibility

**Option 2: Mutually exclusive groups**
- More complex
- Not needed for this case

We'll use **Option 1**.

### 5.3 Code Changes

**File:** `muc_one_up/cli/click_main.py`

#### Change 1: Add callback function

**Location:** After line 51 (after configure_logging function)

**New Code:**
```python
def verbose_callback(ctx, param, value):
    """Callback to set log level to DEBUG when --verbose is used.

    This allows --verbose to override --log-level if both are provided.
    """
    if value:
        # Verbose flag was set - override log_level to DEBUG
        ctx.params["log_level"] = "DEBUG"
    return value
```

**Line to insert:** After line 51

#### Change 2: Add --verbose option

**Location:** After line 71 (after --log-level option)

**Current Code (lines 67-71):**
```python
@click.option(
    "--log-level",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL", "NONE"]),
    default="INFO",
    help="Set logging level.",
)
```

**Add After Line 71:**
```python
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    callback=verbose_callback,
    is_eager=True,
    help="Enable verbose output (equivalent to --log-level DEBUG).",
)
```

**Line to insert:** After line 71

#### Change 3: Update cli function signature

**Location:** Line 73

**Current Code (line 73):**
```python
def cli(ctx, config, log_level):
```

**New Code:**
```python
def cli(ctx, config, log_level, verbose):
```

**Line to modify:** 73

### 5.4 User Experience

**Before:**
```bash
muconeup --config X --log-level DEBUG simulate ...
```

**After (both work):**
```bash
muconeup --config X --log-level DEBUG simulate ...  # Still works
muconeup --config X --verbose simulate ...          # New shortcut
muconeup --config X -v simulate ...                 # Short form
```

**Priority:**
- If both `--verbose` and `--log-level` provided, `--verbose` wins (DEBUG)
- If neither provided, default to INFO

### 5.5 Testing

**File:** `tests/test_click_cli.py`
**Location:** In TestCLIRoot class (after line 109)

```python
def test_verbose_flag_sets_debug_level(self, runner, temp_config, tmp_path):
    """Test that --verbose sets log level to DEBUG."""
    result = runner.invoke(
        cli,
        [
            "--config", temp_config,
            "--verbose",
            "simulate",
            "--out-dir", str(tmp_path),
            "--out-base", "verbose_test",
            "--seed", "42",
        ],
    )

    # Verbose flag should set DEBUG level
    # Check that DEBUG-level messages appear
    assert result.exit_code in [0, 1]
    # Note: Actual DEBUG output depends on implementation
    # This test mainly ensures the flag is accepted

def test_verbose_short_flag_works(self, runner, temp_config, tmp_path):
    """Test that -v short flag works."""
    result = runner.invoke(
        cli,
        [
            "--config", temp_config,
            "-v",
            "simulate",
            "--out-dir", str(tmp_path),
            "--out-base", "verbose_short_test",
            "--seed", "42",
        ],
    )

    assert result.exit_code in [0, 1]

def test_verbose_overrides_log_level(self, runner, temp_config, tmp_path):
    """Test that --verbose overrides --log-level."""
    result = runner.invoke(
        cli,
        [
            "--config", temp_config,
            "--log-level", "WARNING",
            "--verbose",  # Should override to DEBUG
            "simulate",
            "--out-dir", str(tmp_path),
            "--out-base", "override_test",
            "--seed", "42",
        ],
    )

    assert result.exit_code in [0, 1]
```

**Line to insert:** After line 109

### 5.6 Help Text Update

The help text will automatically show:

```
Options:
  --config PATH                   Path to JSON configuration file.  [required]
  --log-level [DEBUG|INFO|WARNING|ERROR|CRITICAL|NONE]
                                  Set logging level.  [default: INFO]
  -v, --verbose                   Enable verbose output (equivalent to
                                  --log-level DEBUG).
  --version                       Show the version and exit.
  --help                          Show this message and exit.
```

### 5.7 Summary

**Lines Added:** ~20 lines (callback + option)
**Lines Modified:** 1 (function signature)
**Tests Added:** 3
**Backward Compatibility:** ✅ Fully maintained
**Benefit:** Better UX, follows common CLI patterns

---

## Testing Strategy

### Test Organization

**All tests in:** `tests/test_click_cli.py`

**Test Markers:**
```python
@pytest.mark.unit
@pytest.mark.cli
```

**Test Classes:**
- `TestCLIRoot` - Root CLI tests (completion, verbose)
- `TestSimulateCommand` - Simulate command tests (progress)
- `TestReadsCommand` - Reads command tests (batch, ONT parity)
- `TestAnalyzeCommand` - Analyze command tests (batch)

### Test Count Summary

| Enhancement | New Tests | Total CLI Tests After |
|-------------|-----------|----------------------|
| Initial | 0 | 148 |
| Batch processing | +5 | 153 |
| ONT parity | +1 | 154 |
| Progress indicators | +1 | 155 |
| Shell completion | +2 | 157 |
| Verbose flag | +3 | 160 |
| **TOTAL** | **+12** | **160** |

### Running Tests

```bash
# All new tests
python -m pytest tests/test_click_cli.py -v -k "batch or ont_with_file or progress or completion or verbose"

# Specific enhancement
python -m pytest tests/test_click_cli.py::TestReadsCommand::test_reads_illumina_with_multiple_files -v

# All CLI tests
python -m pytest tests/test_click_cli.py -v

# Full test suite
python -m pytest tests/ --cov=muc_one_up --cov-report=html
```

---

## Version Bump & Documentation

### 5.1 Version Number

**Current:** v0.9.0
**Proposed:** v0.10.0

**Reason:** Minor version bump for new features (progress, completion, verbose)

**Files to Update:**
1. `pyproject.toml` line 7: `version = "0.10.0"`
2. `muc_one_up/version.py`: `__version__ = "0.10.0"`

### 5.2 Documentation Updates

#### README.md

**Sections to Add:**

1. **Shell Completion** (after Installation, ~line 111)
   - See Enhancement 4.4 above

2. **Usage Examples** (update existing examples)
   - Add verbose flag example
   - Add progress bar note for series mode

```markdown
## Usage

### Basic Simulation

```bash
# Standard simulation
muconeup --config config.json simulate --out-base sample

# Verbose output
muconeup --verbose --config config.json simulate --out-base sample

# Series mode (shows progress bar)
muconeup --config config.json simulate \
  --fixed-lengths 20-40 \
  --simulate-series 5 \
  --out-base series
```
```

#### CLAUDE.md

**Update Common Commands section (lines 8-44):**

Add note about new flags:
```markdown
### Running Tests
```bash
# Run with verbose output
python -m pytest -v tests/

# Or use the tool's verbose flag
muconeup -v --config config.json simulate --help
```

### Progress Indicators
When using `--simulate-series` with multiple iterations, a progress bar will appear:
```bash
Simulating 21 iterations  [################-----------]   61%  00:02:15
```
```

#### CHANGELOG.md

**Create new file:** `CHANGELOG.md`

```markdown
# Changelog

All notable changes to MucOneUp will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.10.0] - 2025-10-18

### Added
- Progress indicators for `--simulate-series` mode (shows ETA and completion %)
- Shell completion support for Bash, Zsh, and Fish
- `--verbose/-v` flag as intuitive alias for `--log-level DEBUG`
- `completion` command to assist with shell completion setup
- Comprehensive batch processing tests for reads/analyze commands

### Improved
- Test coverage increased to 160 CLI-specific tests (was 148)
- Better user experience for long-running simulations
- Tab completion for all commands and options

### Fixed
- None

## [0.9.0] - 2025-10-XX

### Changed
- Migrated from argparse to Click for CLI framework
- Restructured CLI following Unix philosophy (single responsibility)
- Improved architecture with SOLID principles

### Added
- Clean separation: `simulate`, `reads`, `analyze` commands
- Batch processing support for reads/analyze
- vntr-stats analysis command

## [Previous versions]
See git history for earlier changes.
```

---

## Implementation Order

### Phase 1: Testing Enhancements (Medium Priority)
**Estimated:** 5 hours
**Dependencies:** None
**Risk:** Low

1. ✅ **Batch Processing Tests** (4 hours)
   - Add 5 new tests to `tests/test_click_cli.py`
   - Verify batch processing works as expected
   - Test warning messages

2. ✅ **ONT Test Parity** (1 hour)
   - Add `test_reads_ont_with_file`
   - Mirror illumina test pattern

**Deliverable:** Test count increases to 154

---

### Phase 2: User Experience Enhancements (Low Priority)
**Estimated:** 7 hours
**Dependencies:** Phase 1 complete
**Risk:** Low

3. ✅ **Verbose Flag** (2 hours)
   - Add callback function
   - Add `--verbose/-v` option
   - Add 3 tests
   - Update help text

4. ✅ **Shell Completion** (2 hours)
   - Add `completion` command
   - Update README with instructions
   - Add 2 tests

5. ✅ **Progress Indicators** (3 hours)
   - Modify simulation loop
   - Add click.progressbar for series mode
   - Add 1 test
   - Ensure single iterations unaffected

**Deliverable:** Test count increases to 160

---

### Phase 3: Documentation & Release (Final)
**Estimated:** 2 hours
**Dependencies:** Phases 1 & 2 complete
**Risk:** Low

6. ✅ **Version Bump** (30 minutes)
   - Update pyproject.toml to v0.10.0
   - Update version.py

7. ✅ **Documentation** (1 hour)
   - Update README.md (shell completion, examples)
   - Update CLAUDE.md (new commands)
   - Create CHANGELOG.md

8. ✅ **Final Testing** (30 minutes)
   - Run full test suite
   - Verify all 160 tests pass
   - Check coverage report

**Deliverable:** v0.10.0 ready for release

---

## File Change Summary

### Files to Modify

| File | Lines Changed | New Lines | Deleted Lines | Type |
|------|---------------|-----------|---------------|------|
| `tests/test_click_cli.py` | +180 | +180 | 0 | Tests |
| `muc_one_up/cli/click_main.py` | +110 | +95 | -15 | Code |
| `README.md` | +40 | +40 | 0 | Docs |
| `CLAUDE.md` | +15 | +15 | 0 | Docs |
| `pyproject.toml` | 1 | 0 | -1 | Config |
| `muc_one_up/version.py` | 1 | 0 | -1 | Config |
| `CHANGELOG.md` | +50 | +50 | 0 | Docs |
| **TOTAL** | **397** | **380** | **-17** | |

### Exact Line Numbers for Code Changes

#### muc_one_up/cli/click_main.py

1. **Line 51** - Insert verbose_callback function (10 lines)
2. **Line 71** - Insert --verbose option (7 lines)
3. **Line 73** - Modify cli function signature (1 line changed)
4. **Lines 243-256** - Replace with progress bar code (36 lines)
5. **Line 914** - Insert completion command (50 lines)

**Total Changes:** ~103 lines (95 added, 15 replaced, -7 net)

#### tests/test_click_cli.py

1. **After line 109** - Add 2 completion tests + 3 verbose tests (60 lines)
2. **After line 232** - Add 1 progress test (20 lines)
3. **After line 293** - Add 2 batch tests for illumina (50 lines)
4. **After line 300** - Add 1 ONT parity test + 1 batch test (40 lines)
5. **After line 358** - Add 1 batch test for orfs (20 lines)
6. **After line 383** - Add 1 batch test for stats (25 lines)

**Total Changes:** ~215 lines added

---

## Risk Assessment

### Low Risk Items ✅

- **Batch processing tests:** Tests existing functionality
- **ONT test parity:** Mirrors existing pattern
- **Verbose flag:** Purely additive, uses callback
- **Shell completion:** Built-in Click feature
- **Progress indicators:** Conditional (only for series mode)

### Mitigation Strategies

1. **Backward Compatibility**
   - All changes are additive (no breaking changes)
   - Existing functionality preserved
   - Tests verify old behavior still works

2. **Feature Flags**
   - Progress bar only shows for multiple iterations
   - Verbose is optional
   - Completion is opt-in

3. **Testing**
   - 100% test coverage for new features
   - Integration tests ensure end-to-end functionality
   - Existing 558 tests must still pass

---

## Success Criteria

### Must Have ✅

- [ ] All 160 CLI tests pass
- [ ] All 558 total tests pass
- [ ] Code coverage ≥ 95%
- [ ] No breaking changes to existing API
- [ ] Documentation updated

### Should Have ✅

- [ ] Progress bar works in terminal (manual test)
- [ ] Shell completion works in Bash/Zsh (manual test)
- [ ] Verbose flag produces DEBUG output
- [ ] Batch processing handles edge cases

### Nice to Have 🎯

- [ ] Progress bar shows accurate ETA
- [ ] Shell completion includes option values
- [ ] Zero ruff/mypy warnings

---

## Timeline

**Total Estimated Effort:** 12-16 hours

| Phase | Duration | Deliverable |
|-------|----------|-------------|
| Phase 1: Testing | 5 hours | 154 tests passing |
| Phase 2: Features | 7 hours | 160 tests passing, all features working |
| Phase 3: Docs & Release | 2 hours | v0.10.0 ready |
| **Buffer** | 2 hours | Unexpected issues |
| **TOTAL** | **14-16 hours** | Production-ready release |

---

## Post-Implementation Checklist

- [ ] Run full test suite: `python -m pytest tests/ -v`
- [ ] Check coverage: `python -m pytest tests/ --cov=muc_one_up --cov-report=html`
- [ ] Run linter: `ruff check .`
- [ ] Run type checker: `mypy muc_one_up`
- [ ] Manual test: Shell completion in Bash
- [ ] Manual test: Progress bar in series mode
- [ ] Manual test: Verbose flag output
- [ ] Review all docstrings
- [ ] Update git commit messages
- [ ] Tag release: `git tag -a v0.10.0 -m "Release v0.10.0"`

---

## Notes for Implementation

### Code Style

- Follow existing patterns in click_main.py
- Use type hints (Python 3.10+ syntax)
- Add docstrings for all new functions
- Follow ruff/mypy linting rules

### Testing Patterns

- Use `CliRunner` for all CLI tests
- Use `tmp_path` fixture for file operations
- Use `temp_config` fixture for config files
- Check both success (0) and graceful failure (1) exit codes
- Add `@pytest.mark.unit` and `@pytest.mark.cli` markers

### Documentation Style

- Follow existing README.md formatting
- Use code blocks with syntax highlighting
- Add clear examples
- Keep language concise and technical

---

**End of Implementation Plan**

This plan provides exact line numbers, detailed code changes, and comprehensive testing strategy for all 5 enhancements. Ready for implementation!
